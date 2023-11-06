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

import blbutil.Const;
import blbutil.FileIt;
import blbutil.Filter;
import blbutil.InputIt;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import bref.Bref3It;
import ints.IndexArray;
import ints.IntArray;
import ints.IntList;
import java.io.File;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import vcf.FilterUtil;
import vcf.GTRec;
import vcf.GeneticMap;
import vcf.IntervalVcfIt;
import vcf.Marker;
import vcf.RefGTRec;
import vcf.RefIt;
import vcf.Samples;

/**
 * <p>Class {@code AggregateIt} is an iterator for aggregate markers from a VCF
 * file. The Java Virtual Machine will exit with an error message if the
 * constructor or methods of this class encounter an error while reading
 * or parsing an input file.</p>
 *
 * <p>Instances of class {@code AggregateIt} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AggregateIt implements SampleFileIt<AggMarker> {

    private final SampleFileIt<RefGTRec> it;
    private final int nHaps;
    private final int minMac;
    private final GeneticMap genMap;
    private RefGTRec nextRec;
    private int nMarkers = 0;
    private int nFilteredMarkers = 0;
    private final Deque<AggMarker> buffer;
    private final float aggDist;

    /**
     * Creates a {@code AggregateIt} instance from the specified data.
     * @param par the command line parameters
     * @throws NullPointerException if {@code par == null}
     */
    public AggregateIt(ClustPar par) {
        this.it = refIt(par);
        this.nHaps = 2*it.samples().size();
        this.minMac = (int) Math.ceil(Math.nextDown(par.min_maf())*nHaps);
        this.genMap = GeneticMap.geneticMap(par.map(), par.chromInt());
        this.nextRec = it.hasNext() ? it.next() : null;
        this.buffer = new ArrayDeque<>(1<<12);
        this.aggDist = par.aggregate();
        fillBuffer();
    }

    private static SampleFileIt<RefGTRec> refIt(ClustPar par) {
        String filename = par.gt().toString();
        Filter<Marker> mFilter = FilterUtil.markerFilter(par.excludemarkers());
        SampleFileIt<RefGTRec> refIt;
        if (filename.endsWith(".bref3")) {
            if (par.excludesamples()!=null) {
                StringBuilder sb = new StringBuilder(256);
                sb.append("ERROR: the excludesamples parameter cannot be used if the gt parameter");
                sb.append(Const.nl);
                sb.append("       is a bref3 file");
                sb.append(Const.nl);
                sb.append("Suggestion: replace the bref3 file with a bgzip-compressed VCF file");
                Utilities.exit(ClustPar.usage() + sb.toString());
            }
            refIt = new Bref3It(par.gt(), mFilter);
        }
        else {
            Filter<String> sFilter = FilterUtil.sampleFilter(par.excludesamples());
            int nBufferedBlocks = par.nthreads() << 3;
            FileIt<String> it = InputIt.fromBGZipFile(par.gt(), nBufferedBlocks);
            refIt = RefIt.create(it, sFilter, mFilter);
        }
        if (par.chromInt()!=null) {
            refIt = new IntervalVcfIt<>(refIt, par.chromInt());
        }
        return refIt;
    }

    @Override
    public Samples samples() {
        return it.samples();
    }

    /**
     * Returns the number of markers before MAF-filtering and aggregation.
     * @return the number of markers before MAF-filtering and aggregation
     */
    public int nMarkers() {
        return nMarkers;
    }

    /**
     * Returns the number times the number of markers after MAF-filtering
     * and before aggregation.
     * @return the number times the number of markers after MAF-filtering
     * and before aggregation
     */
    public int nFilteredMarkers() {
        return nFilteredMarkers;
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public void close() {
        nextRec=null;
        buffer.clear();
        it.close();
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * (In other words, returns {@code true} if {@link #next} would
     * return an element rather than throwing an exception.)
     *
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return !buffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public AggMarker next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        AggMarker first = buffer.removeFirst();
        if (buffer.isEmpty()) {
            fillBuffer();
        }
        return first;
    }

    private void fillBuffer() {
        assert buffer.isEmpty();
        while (nextRec!=null) {
            RefGTRec[] recs = readChrom();
            int[] grpEnds = recGrpEnds(recs);
            List<AggMarker> aggMarkerList = IntStream.range(0, grpEnds.length)
                    .parallel()
                    .mapToObj(j -> aggregate(recs, grpEnds, j))
                    .collect(Collectors.toCollection(ArrayList::new));
            buffer.addAll(aggMarkerList);
        }
    }

    private RefGTRec[] readChrom() {
        List<RefGTRec> recList = new ArrayList<>(1<<14);
        int chromIndex = nextRec.marker().chromIndex();
        while (nextRec!=null && nextRec.marker().chromIndex()==chromIndex) {
            recList.add(nextRec);
            nextRec = it.hasNext() ? it.next() : null;
        }
        RefGTRec[] recs = recList.toArray(new RefGTRec[0]);
        nMarkers += recs.length;
        recs = applyMacFilter(recs, minMac);
        nFilteredMarkers += recs.length;
        return recs;
    }

    private static RefGTRec[] applyMacFilter(RefGTRec[] recs, int minMac) {
        return Arrays.stream(recs)
                .parallel()
                .filter(rec -> mac(rec)>=minMac)
                .toArray(RefGTRec[]::new);
    }

    private static int mac(RefGTRec rec) {
        int[] alCnts = GTRec.alleleCounts(rec);
        Arrays.sort(alCnts);
        return alCnts.length==1 ? 0 : alCnts[alCnts.length-2];
    }

    private int[] recGrpEnds(RefGTRec[] recs) {
        double[] pos = Arrays.stream(recs)
                .parallel()
                .mapToDouble(rec -> genMap.genPos(rec.marker()))
                .toArray();
        IntList il = new IntList(recs.length>>2);

        double maxGenPos = pos[0] + aggDist;
        for (int j=0; j<pos.length; ++j) {
            if (pos[j]>maxGenPos) {
                il.add(j);
                maxGenPos = pos[j] + aggDist;
            }
        }
        il.add(pos.length);
        return il.toArray();
    }

    private AggMarker aggregate(RefGTRec[] recs, int[] recGrpEnds, int i) {
        int start = i==0 ? 0 : recGrpEnds[i-1];
        int end = recGrpEnds[i];
        Position medianPos = medianPos(recs, start, end, genMap);
        return new AggMarker(medianPos, codeHaps(recs, start, end));
    }

    private Position medianPos(RefGTRec[] recList, int start, int end,
            GeneticMap genMap) {
        int n = end - start;
        int i1 = start + ((n-1)>>1);
        int i2 = start + (n>>1);
        Marker mkr1 = recList[i1].marker();
        double genPos1 = genMap.genPos(mkr1);

        if (i1==i2) {
            return new Position(mkr1.chromIndex(), mkr1.pos(), genPos1);
        }
        else {
            Marker mkr2 = recList[i2].marker();
            int medianPos = (mkr1.pos() + mkr2.pos()) >>> 1;
            double medianGenPos = 0.5*(genPos1 + genMap.genPos(mkr2));
            return new Position(mkr1.chromIndex(), medianPos, medianGenPos);
        }
    }

    private IndexArray codeHaps(RefGTRec[] recs, int start, int end) {
        int notMapped = -1;
        int[] hapToSeq = new int[nHaps];
        int[] seqMap = new int[(end-start)<<3];
        int seqCnt = 1;
        for (int j=start; j<end; ++j) {
            GTRec rec = recs[j];
            int nAlleles = rec.marker().nAlleles();
            int size = seqCnt*nAlleles;
            if (size>seqMap.length) {
                seqMap = new int[Math.max(size, seqMap.length<<1)];
            }
            Arrays.fill(seqMap, 0, size, notMapped);
            seqCnt = 0;
            for (int h=0; h<nHaps; ++h) {
                int i = nAlleles*hapToSeq[h] + rec.get(h);
                if (seqMap[i] == notMapped) {
                    seqMap[i] = seqCnt++;
                }
                hapToSeq[h] = seqMap[i];
            }
        }
        IntArray intArray = IntArray.packedCreate(hapToSeq, seqCnt);
        return new IndexArray(intArray, seqCnt);
    }
}
