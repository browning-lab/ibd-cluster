/*
 * Copyright 2023 Brian L. Browning
 *
 * This file is part of the ibd-cluster program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package ibdcluster;

import beagleutil.PbwtDivUpdater;
import blbutil.BGZIPOutputStream;
import blbutil.Const;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import ints.IndexArray;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.stream.IntStream;
import java.util.NoSuchElementException;
import vcf.Samples;

/**
 * <p>Class {@code Cluster} performs an ibd-cluster analysis.</p>
 *
 * <p>Instances of class {@code Cluster} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Cluster {

    private static final int MAX_BUFFERED_BYTES_PER_THREAD = 1 << 27;

    private final int nThreads;
    private final int maxBufferedSitesPerThread;
    private final int maxBufferedSites;

    private final float minIbsCm;
    private final float trimCm;

    private final int nHaps;
    private final PbwtDivUpdater pbwt;
    private final int[] prefix;
    private final int[] div;
    private final EndList partitions;

    private int chrom = -1;
    private AggMarker next;

    private int windowStart = 0;
    private int trimStart = 0;
    private int trimEnd = 0;            // exclusive end
    private int prevWindowStart = 0;
    private int prevTrimEnd = 0;

    private int nAggregateMarkers = 0;
    private int nIbdSets = 0;

    /**
     * Constructs a new {@code Clust} class for the specified data.
     * @param minIbsCm the minimum IBS segment length in cM
     * @param trimCm the cM length trimmed from segment ends that do not
     * include the first or last marker on the chromosome
     * @param nHaps the number of haplotypes
     * @param nThreads the number of threads used to compress otput records
     * @throws IllegalArgumentException if
     * {@code length <= 0.0 || Float.isFinite(length) == false}
     * @throws IllegalArgumentException if
     * {@code trim <= 0.0 || Float.isFinite(trim) == false}
     * @throws IllegalArgumentException if {@code nHaps < 0}
     * @throws IllegalArgumentException if {@code nThreads < 1}
     */
    public Cluster(float minIbsCm, float trimCm, int nHaps, int nThreads) {
        if (minIbsCm <= 0.0 || Float.isFinite(minIbsCm)==false) {
            throw new IllegalArgumentException(String.valueOf(minIbsCm));
        }
        if (trimCm <= 0.0 || Float.isFinite(trimCm)==false) {
            throw new IllegalArgumentException(String.valueOf(trimCm));
        }
        if (nThreads < 1) {
            throw new IllegalArgumentException(String.valueOf(nThreads));
        }
        this.nThreads = nThreads;
        // singleton clusters are 8 bytes per hap
        this.maxBufferedSitesPerThread
                = Math.max(((MAX_BUFFERED_BYTES_PER_THREAD>>3) / nHaps), 16);
        this.maxBufferedSites = nThreads*maxBufferedSitesPerThread;
        this.minIbsCm = minIbsCm;
        this.trimCm = trimCm;
        this.nHaps = nHaps;
        this.pbwt = new PbwtDivUpdater(nHaps);
        this.partitions = new EndList(1 << 9);
        this.prefix = new int[nHaps];
        this.div = new int[nHaps];
    }

    /**
     * Runs an ibd-cluster analysis.  If an I/O error occurs during this
     * method, the Java Virtual Machine will exit with an error message.
     * @param it an iterator returns VCF records
     * @param out the PrintWriter to which write partitions will be written

     * @throws IllegalArgumentException if
     * {@code 2*it.samples().size() != this.nHaps()}
     * @throws NoSuchElementException if {@code it.hasNext() == false}
     * @throws NullPointerException if {@code it == null || writer == null}
     */
    public void run(SampleFileIt<AggMarker> it, OutputStream out) {
        printHeader(it.samples(), out);
        next = it.next();
        while (next!=null) {
            analyzeChrom(it, out);
        }
        writeEmptyBGZipBlock(out);
    }

    private void printHeader(Samples samples, OutputStream out) {
        PrintWriter pw = new PrintWriter(new BGZIPOutputStream(out, false));
        String[] ids = samples.ids();
        if ((ids.length<<1)!=nHaps) {
            throw new IllegalArgumentException(String.valueOf(nHaps));
        }
        char tab = Const.tab;
        pw.print("CHROM");
        pw.print(tab);
        pw.print("POS");
        pw.print(tab);
        pw.print("CM");
        for (String id : ids) {
            pw.print(tab);
            pw.print(id);
        }
        pw.println();
        pw.flush();
    }

    private void analyzeChrom(SampleFileIt<AggMarker> it, OutputStream out) {
        initializeFields();
        int outCnt = partitions.end();
        chrom = next.position().chrom();
        double targetGenPos = next.position().genPos() + minIbsCm;
        while (next!=null && next.position().chrom()==chrom
                && next.position().genPos()<targetGenPos) {
            advanceWindow(it);
        }
        while (next!=null && next.position().chrom()==chrom) {
            advanceWindow(it);
            mergeClusters();
            int n = windowStart - outCnt;
            if (n >= maxBufferedSites) {
                writeAndRemovePartitions(n, out);
                outCnt = windowStart;
            }
        }
        int n = partitions.end() - outCnt;
        writeAndRemovePartitions(n, out);
    }

    private void initializeFields() {
        int cnt = partitions.end();
        prevWindowStart = cnt;
        prevTrimEnd = cnt;
        windowStart = cnt;
        trimStart = cnt;
        trimEnd = cnt;
        for (int j=0; j<prefix.length; ++j) {
            prefix[j] = j;
            div[j] = cnt;
        }
    }

   private void advanceWindow(SampleFileIt<AggMarker> it) {
        if (windowStart < partitions.end()) {
            double endGenPos = next.position().genPos();
            prevWindowStart = windowStart;
            prevTrimEnd = trimEnd;
            windowStart = updateIndexLE(windowStart, endGenPos - minIbsCm);
            double startGenPos = partitions.get(windowStart).position().genPos();
            trimStart = updateIndexGE(trimStart, startGenPos + trimCm);
            trimEnd = updateIndexLE(trimEnd-1, (endGenPos-trimCm)) + 1;  // trimEnd is exclusive end
        }
        IndexArray alleles = next.alleles();
        pbwt.fwdUpdate(alleles.intArray(), alleles.valueSize(),
                nAggregateMarkers++, prefix, div);
        partitions.add(new Partition(next.position(), nHaps));
        next = it.hasNext() ? it.next() : null;
    }

    private int updateIndexLE(int index, double maxGenPos) {
        int indexP1 = index+1;
        while(indexP1<partitions.end()
                && partitions.get(indexP1).position().genPos() <= maxGenPos) {
            ++index;
            indexP1 = index + 1;
        }
        return index;
    }

    private int updateIndexGE(int index, double minGenPos) {
        while(index<partitions.end()
                && partitions.get(index).position().genPos() < minGenPos) {
            ++index;
        }
        return index;
    }

    private void mergeClusters() {
        for (int j=1; j<nHaps; ++j) {
            if (div[j]<=windowStart) {
                int start = trimStart;
                if (div[j]<=prevWindowStart) {
                    start = prevTrimEnd;
                }
                for (int m=start; m<trimEnd; ++m) {
                    partitions.get(m).union(prefix[j-1], prefix[j]);
                }
            }
        }
    }

    private void writeAndRemovePartitions(int nPartitions, OutputStream out) {
        int batchSize = (nPartitions + nThreads - 1) / nThreads;
        byte[][] compressed = IntStream.range(0, nThreads)
                .parallel()
                .mapToObj(batch -> bgzipCompress(partitions, nPartitions, batch, batchSize))
                .toArray(byte[][]::new);
        try {
            for (byte[] ba : compressed) {
                out.write(ba);
            }
        }
        catch (IOException ex) {
            Utilities.exit(ex, "Error writing output IBD cluster file");
        }
        for (int j=0; j<nPartitions; ++j) {
            Partition p = partitions.removeStart();
            nIbdSets += p.nSets();
        }
    }

    private static byte[] bgzipCompress(EndList partitions, int nPartitions,
            int batch, int batchSize) {
        int partitionsStart = partitions.start();
        int start = partitionsStart + batch*batchSize;
        int end = Math.min(start+batchSize, (partitionsStart + nPartitions));
        boolean writeEmptyBlock = false;
        ByteArrayOutputStream baos = new ByteArrayOutputStream(batchSize<<8);
        try (BGZIPOutputStream bgzos = new BGZIPOutputStream(baos, writeEmptyBlock);
                PrintWriter out = new PrintWriter(bgzos)) {
            for (int j=start; j<end; ++j) {
                partitions.get(j).write(out);
            }
        } catch (IOException ex) {
            Utilities.exit(ex);
        }
        return baos.toByteArray();
    }

    private static void writeEmptyBGZipBlock(OutputStream out) {
        try {
            BGZIPOutputStream.writeEmptyBlock(out);
        } catch (IOException ex) {
            Utilities.exit(ex, "Error writing output IBD cluster file");
        }
    }

    /**
     * Returns the number of analyzed VCF variants after MAC-filtering.
     * Returns 0 if {@code this.run()} has not been invoked.
     *
     * @return the number of analyzed VCF variants
     */
    public int nAggregateMarkers() {
        return nAggregateMarkers;
    }

    /**
     * Returns the total number of IBD sets for all aggregate markers.
     * Returns 0 if {@code this.run()} has not been invoked.
     * @return the total number of IBD sets for all aggregate markers
     */
    public long nIbdSets() {
        return nIbdSets;
    }
}
