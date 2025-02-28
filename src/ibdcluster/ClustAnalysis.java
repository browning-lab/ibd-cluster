/*
 * Copyright 2023-2025 Brian L. Browning
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

import blbutil.BGZIPOutputStream;
import blbutil.Const;
import blbutil.FileUtil;
import blbutil.Utilities;
import ints.WrappedIntArray;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import vcf.MarkerMap;
import vcf.PlinkGenMap;
import vcf.RefGT;
import vcf.Samples;

/**
 * <p>Class {@code ClustAnalysis} runs an IBD clustering analysis.</p>
 *
 * <p>Instances of class {@code ClustAnalysis} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ClustAnalysis {

    private ClustAnalysis() {
        // private constructor to prevent instantiation
    }

    /**
     * Runs an IBD clustering analysis.
     * @param par the command line parameters
     * @return analysis statistics
     * @throws NullPointerException if {@code par == null}
     */
    public static ClustStats run(ClustPar par) {
        ClustStats stats = new ClustStats();
        File outFile = new File(par.out() + ".ibdclust.gz");
        try (OutputStream os = FileUtil.bufferedOutputStream(outFile);
                ChromIt it = new ChromIt(par)) {
            PlinkGenMap genMap = it.genMap();
            printHeader(it.samples(), os);
            while (it.hasNext()) {
                RefGT refGT = it.next();
                analyzeChrom(par, refGT, genMap, stats, os);
            }
            BGZIPOutputStream.writeEmptyBlock(os);
            stats.setNSamples(it.samples().size());
            stats.addMarkers(it.nMarkers());
            stats.addFilteredMarkers(it.nFilteredMarkers());
        }
        catch (IOException ex) {
            Utilities.exit(ex, "Error writing output IBD cluster file");
        }
        return stats;
    }

    private static void printHeader(Samples samples, OutputStream os) {
        PrintWriter out = new PrintWriter(new BGZIPOutputStream(os, false));
        out.print("CHROM");
        out.print(Const.tab);
        out.print("POS");
        out.print(Const.tab);
        out.print("CM");
        String[] ids = samples.ids();
        for (String id : ids) {
            out.print(Const.tab);
            out.print(id);
        }
        out.println();
        out.flush();
    }

    private static void analyzeChrom(ClustPar par, RefGT refGT,
            PlinkGenMap genMap, ClustStats stats, OutputStream os) {
        ClustData data = new ClustData(par, refGT, genMap);
        HapPairSegment[] ibdSegs = ibdSegs(data, stats);
        clusterHaps(data, stats, ibdSegs, os);
    }

    private static HapPairSegment[] ibdSegs(ClustData data, ClustStats stats) {
        RefGT refGT = data.fwdGT();
        ClustPar par = data.par();
        MarkerMap map = MarkerMap.create(data.genMap(), refGT.markers());
        IbsSegments ibs = new IbsSegments(par, refGT, map.genPos());
        HapPairSegment[] ibsSegs = ibs.hapPairSegments();
        ConcurrentLinkedQueue<IbdEstimator> ibdEstQ = ibdEstQ(data, stats);
        HapPairSegment[] ibdSegs = Arrays.stream(ibsSegs)
                .parallel()
                .map(hps -> ibsToIbd(data, stats, hps, ibdEstQ))
                .filter(hps -> hps!=HapPairSegment.ZERO_LENGTH_SEGMENT)
                .toArray(HapPairSegment[]::new);
        return ibdSegs;
    }

    private static HapPairSegment ibsToIbd(ClustData clustData,
            ClustStats clustStats, HapPairSegment hps,
            ConcurrentLinkedQueue<IbdEstimator> ibdEstQ) {
        IbdEstimator ibdEst = ibdEstQ.poll();
        if (ibdEst==null) {
            ibdEst = new IbdEstimator(clustData, clustStats);
        }
        HapPairSegment ibdSeg = ibdEst.ibdSegment(hps);
        ibdEstQ.add(ibdEst);
        return ibdSeg;
    }

    private static ConcurrentLinkedQueue<IbdEstimator> ibdEstQ(ClustData data,
            ClustStats stats) {
        List<IbdEstimator> list = IntStream.rangeClosed(0, data.par().nthreads())
                .parallel()
                .mapToObj(j -> new IbdEstimator(data, stats))
                .collect(Collectors.toList());
        return new ConcurrentLinkedQueue<>(list);
    }

    private static void clusterHaps(ClustData data, ClustStats stats,
            HapPairSegment[] ibdSegs, OutputStream os) {
        ClustPar par = data.par();
        ibdSegs = sortSegments(data, ibdSegs);
        recordDiscordRate(data, stats, ibdSegs);
        int sitesPerWindow = par.out_window_size();
        double outMorgans = 0.01*par.out_cm();
        double startMorgans = data.morganPos().get(0);
        double endMorgans = data.morganPos().get(data.morganPos().size()-1);
        int fromStepIndex = (int) Math.ceil(startMorgans/outMorgans);
        int toStepIndex = (int) Math.ceil(endMorgans/outMorgans); // exclusive end
        stats.addOutputPositions(toStepIndex - fromStepIndex);
        for (int start=fromStepIndex; start<toStepIndex; start+=sitesPerWindow) {
            int end = Math.min(start + sitesPerWindow, toStepIndex);
            Partition[] partitions  = clusterHaps(data, ibdSegs, outMorgans, start, end);
            recordPartitions(partitions, stats);
            writePartitions(partitions, par.nthreads(), os);
            int minInclEnd = data.morganToBase(end*outMorgans);
            ibdSegs = filterSegments(ibdSegs, minInclEnd);
        }
    }

    private static HapPairSegment[] sortSegments(ClustData data,
            HapPairSegment[] ibdSegs) {
        ibdSegs = Arrays.stream(ibdSegs)
                .parallel()
                .sorted(HapPairSegment.intervalComparator())
                .toArray(HapPairSegment[]::new);
        return ibdSegs;
    }


    private static void recordDiscordRate(ClustData data, ClustStats stats,
            HapPairSegment[] ibdSegs) {
        Arrays.stream(ibdSegs)
                .parallel()
                .forEach(hps -> recordDisordRate(data, stats, hps));
    }

    private static void recordDisordRate(ClustData data, ClustStats stats,
            HapPairSegment hps) {
        int hap1 = hps.hap1();
        int hap2 = hps.hap2();
        int startPos = hps.startPos();
        int inclEndPos = hps.inclEndPos();
        assert startPos <= inclEndPos;
        WrappedIntArray basePos = data.basePos();
        RefGT refGT = data.fwdGT();
        int startMarker = basePos.binarySearch(startPos);
        if (startMarker<0) {
            startMarker = -startMarker-1;
        }
        int endMarker = basePos.binarySearch(startMarker, basePos.size(), inclEndPos);
        if (endMarker<0) {
            endMarker = -endMarker-2; // -2 so that marker is within interval
        }
        if (startMarker<=endMarker ) {
            int discordCnt = 0;
            for (int m=startMarker; m<=endMarker; ++m) {
                if (refGT.allele(m, hap1)!=refGT.allele(m, hap2)) {
                    ++discordCnt;
                }
            }
            int totalChecked = endMarker - startMarker + 1;
            stats.updateDiscordRate(discordCnt, totalChecked);
        }
    }

    private static Partition[] clusterHaps(ClustData data,
            HapPairSegment[] ibdSegs, double outMorgans,
            int startStep, int endStep) {
        return IntStream.range(startStep, endStep)
            .parallel()
            .mapToObj(j -> cluster(data, j*outMorgans, ibdSegs))
            .toArray(Partition[]::new);
    }

    private static Partition cluster(ClustData data, double morganPos,
            HapPairSegment[] ibdSegs) {
        int basePos = data.morganToBase(morganPos);
        double cmPos = 100*morganPos;
        int nHaps = data.fwdGT().nHaps();
        Partition p = new Partition(data.chrom(), basePos, cmPos, nHaps);
        int j=0;
        HapPairSegment hps = j<ibdSegs.length ? ibdSegs[j++] : null;
        while (hps!=null && hps.startPos() <= p.pos()) {
            if (p.pos() <= hps.inclEndPos()) {
                p.union(hps.hap1(), hps.hap2());
            }
            hps = j<ibdSegs.length ? ibdSegs[j++] : null;
        }
        return p;
    }

    private static void recordPartitions(Partition[] partitions, ClustStats stats) {
        long nIbdSets = Arrays.stream(partitions)
                .parallel()
                .mapToLong(p -> p.nSets())
                .sum();
        stats.addIbdSets(nIbdSets);
    }

    private static void writePartitions(Partition[] partitions, int nThreads,
            OutputStream out) {
        int batchSize = (partitions.length + nThreads - 1) / nThreads;
        byte[][] compressed = IntStream.range(0, nThreads)
                .parallel()
                .mapToObj(batch -> bgzipCompress(partitions, batch, batchSize))
                .toArray(byte[][]::new);
        try {
            for (byte[] ba : compressed) {
                out.write(ba);
            }
        }
        catch (IOException ex) {
            Utilities.exit(ex, "Error writing output IBD cluster file");
        }
    }

    private static byte[] bgzipCompress(Partition[] partitions, int batch, int batchSize) {
        int start = batch*batchSize;
        int end = Math.min(start+batchSize, partitions.length);
        boolean writeEmptyBlock = false;
        ByteArrayOutputStream baos = new ByteArrayOutputStream(batchSize<<8);
        try (BGZIPOutputStream bgzos = new BGZIPOutputStream(baos, writeEmptyBlock);
                PrintWriter out = new PrintWriter(bgzos)) {
            for (int j=start; j<end; ++j) {
                partitions[j].write(out);
            }
        } catch (IOException ex) {
            Utilities.exit(ex);
        }
        return baos.toByteArray();
    }

    private static HapPairSegment[] filterSegments(HapPairSegment[] ibdSegs,
            int minInclEnd) {
        return Arrays.stream(ibdSegs)
                .parallel()
                .filter(hps -> (hps.inclEndPos() >= minInclEnd))
                .toArray(HapPairSegment[]::new);
    }
}
