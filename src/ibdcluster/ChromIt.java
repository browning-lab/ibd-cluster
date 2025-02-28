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

import beagleutil.ChromIds;
import blbutil.Const;
import blbutil.FileIt;
import blbutil.InputIt;
import blbutil.SampleFileIt;
import blbutil.Utilities;
import blbutil.VcfFileIt;
import bref.Bref3It;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.Predicate;
import vcf.FilterUtil;
import vcf.GTRec;
import vcf.IntervalVcfIt;
import vcf.Marker;
import vcf.PlinkGenMap;
import vcf.RefGT;
import vcf.RefGTRec;
import vcf.RefIt;
import vcf.Samples;

/**
 * <p>Class {@code ChromIt} is an iterator whose {@code it.next()} method
 * return MAF-filtered, phased genotype data for the next chromosome
 * in an input VCF file.
 * The Java Virtual Machine will exit with an error message if the
 * constructor or methods of this class encounter an error while reading
 * or parsing an input file.</p>
 *
 * <p>Instances of class {@code ChromIt} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ChromIt implements SampleFileIt<RefGT> {

    private final ClustPar par;
    private final SampleFileIt<RefGTRec> it;
    private final int nHaps;
    private final int minMac;
    private final PlinkGenMap genMap;
    private final HashSet<Integer> chroms;
    private RefGTRec next;

    private long nMarkers = 0;
    private long nFilteredMarkers = 0;

    /**
     * Creates a {@code ChromIt} instance from the specified data.
     * @param par the command line parameters
     * @throws NullPointerException if {@code par == null}
     */
    public ChromIt(ClustPar par) {
        this.par = par;
        this.it = refIt(par);
        this.nHaps = 2*it.samples().size();
        this.minMac = (int) Math.ceil(Math.nextDown(par.min_maf())*nHaps);
        this.genMap = par.chrom()==null
                ? PlinkGenMap.fromPlinkMapFile(par.map())
                : PlinkGenMap.fromPlinkMapFile(par.map(), par.chromInt().chrom());
        if (it.hasNext()==false) {
                Utilities.exit(Const.nl
                        + "Error: no VCF records found after filtering. Exiting program.");
        }
        this.chroms = new HashSet<>();
        this.next = it.next();
    }

    private static VcfFileIt<RefGTRec> refIt(ClustPar par) {
        String filename = par.gt().toString();
        Predicate<Marker> mFilter = FilterUtil.markerFilter(par.excludemarkers());
        Predicate<String> sFilter = FilterUtil.excludePredicate(par.excludesamples());
        VcfFileIt<RefGTRec> refIt;
        if (filename.endsWith(".bref3")) {
            refIt = new Bref3It(par.gt(), sFilter, mFilter);
        }
        else {
            int nBufferedBlocks = par.nthreads() << 2;
            FileIt<String> it = InputIt.fromBGZipFile(par.gt(), nBufferedBlocks);
            refIt = RefIt.create(it, sFilter, mFilter);
        }
        if (par.chromInt()!=null) {
            refIt = new IntervalVcfIt<>(refIt, par.chromInt());
        }
        return refIt;
    }

    private RefGT nextChrom() {
        int chromIndex = nextChromIndex();
        int lastMapIndex = genMap.nMapPositions(chromIndex) - 1;
        int firstMapPos = genMap.index2BasePos(chromIndex, 0);
        int lastMapPos = genMap.index2BasePos(chromIndex, lastMapIndex);
        List<RefGTRec> list = new ArrayList<>(8192);
        try {
            while (next!=null && next.marker().chromIndex()==chromIndex) {
                int pos = next.marker().pos();
                if ((firstMapPos <= pos) && (pos <=lastMapPos)) {
                    list.add(next);
                }
                next = it.hasNext() ? it.next() : null;
            }
        } catch (Throwable t) {
            Utilities.exit(t);
        }
        nMarkers += list.size();
        RefGTRec[] recs = applyMacFilter(list, minMac);
        nFilteredMarkers += recs.length;
        if (recs.length==0) {
            exitWithNoRecordsErr(chromIndex);
        }
        return new RefGT(recs);
    }

    private int nextChromIndex() {
        assert next!=null;
        int chromIndex = next.marker().chromIndex();
        boolean newChrom = chroms.add(chromIndex);
        if (newChrom==false) {
            exitWithNoncontiguousChromErr(chromIndex);
        }
        return chromIndex;
    }

    private static RefGTRec[] applyMacFilter(List<RefGTRec> list, int minMac) {
        if (minMac>0) {
            return list.stream()
                    .parallel()
                    .filter(rec -> mac(rec)>=minMac)
                    .toArray(RefGTRec[]::new);
        }
        else {
            return list.toArray(new RefGTRec[0]);
        }
    }

    private static int mac(RefGTRec rec) {
        int[] alCnts = GTRec.alleleCounts(rec);
        Arrays.sort(alCnts);
        return alCnts.length<=1 ? 0 : alCnts[alCnts.length-2];
    }

    private void exitWithNoRecordsErr(int chromIndex) {
        String s = Const.nl
                + "ERROR: there are no VCF records inside the boundaries of the genetic"
                + Const.nl
                + "map for chromosome "
                + ChromIds.instance().id(chromIndex)
                + " after minor allele frequency filtering.";
        Utilities.exit(s);
    }

    private void exitWithNoncontiguousChromErr(int chromIndex) {
        String s = Const.nl
                + "ERROR: the VCF records for chromosome "
                + ChromIds.instance().id(chromIndex)
                + " are not contiguous.";
        Utilities.exit(s);
    }

    /**
     * Returns the genetic map.
     * @return the genetic map
     */
    public PlinkGenMap genMap() {
        return genMap;
    }

    /**
     * Returns the command line parameters.
     * @return the command line parameters
     */
    public ClustPar par() {
        return par;
    }

    /**
     * Returns the cumulative number of markers read by previous
     * invocations of {@code this.next()}.  The returned count
     * includes markers that were excluded by the MAF filter.
     * @return the cumulative number of markers read by previous
     * invocations of {@code this.next()}
     */
    public long nMarkers() {
        return nMarkers;
    }

    /**
     * Returns the cumulative number of markers returned by previous
     * invocations of {@code this.next()}. The returned count
     * excludes markers that were excluded by the MAF filter.
     * @return the cumulative number of markers returned by previous
     * invocations of {@code this.next()}
     */
    public long nFilteredMarkers() {
        return nFilteredMarkers;
    }

    @Override
    public Samples samples() {
        return it.samples();
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public void close() {
        next=null;
        it.close();
    }

    /**
     * Returns {@code true} if the iteration has more elements.
     * In other words, returns {@code true} if {@link #next} would
     * return an element rather than throwing an exception.
     *
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return next!=null;
    }

    /**
     * Returns the next element in the iteration.
     *
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements
     */
    @Override
    public RefGT next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        return nextChrom();
    }
}
