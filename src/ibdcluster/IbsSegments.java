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

import beagleutil.PbwtDivUpdater;
import blbutil.DoubleArray;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;
import vcf.RefGT;
import vcf.RefGTRec;

/**
 * <p>Class {@code IbsSegments} identifies IBS haplotype segments in
 * phased genotype data that have no non-missing alleles.
 * </p>
 * <p>Instances of {@code IbsSegments} are immutable.
 * </p>
 * <p>
 * Reference: Durbin, R (2014) Efficient haplotype matching and storage using
 *            the positional Burrows-Wheeler transform (PBWT).
 *            Bioinformatics 30(9):1266-1272.
 *
 * </p>
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class IbsSegments {

    private final ClustPar par;
    private final RefGT refGT;
    private final DoubleArray cmPos;
    private final double minIbsCm;
    private final int nAnalyses;

    /**
     * Constructs a new {@code IbsSegments} object from the specified data.
     * @param par the analysis parameters
     * @param refGT phased, non-missing genotypes
     * @param cmPos an array of non-decreasing marker positions whose
     * {@code j}-th element is the cM position of the {@code j}-th marker
     * in {@code refGT}
     * @throws IllegalArgumentException if
     * {@code refGT.nMarkers() != cmPos.size()}
     * @throws NullPointerException if {@code (refGT == null || cmPos == null)}
     */
    public IbsSegments(ClustPar par, RefGT refGT, DoubleArray cmPos) {
        if (refGT.nMarkers()!= cmPos.size()) {
            throw new IllegalArgumentException(String.valueOf(refGT.nMarkers()));
        }
        this.par = par;
        this.refGT = refGT;
        this.cmPos = cmPos;
        this.minIbsCm = par.min_ibs_cm();
        this.nAnalyses = par.pbwt();
    }

    /**
     * Performs the specified number of interleaved PBWT analyses, and
     * returns an array of haplotype pair segments.
     * @return an array of haplotype pair segments for IBD-based clustering
     */
    public HapPairSegment[] hapPairSegments() {
        HapPairSegment[] segs = IntStream.range(0, nAnalyses)
                .parallel()
                .mapToObj(j -> ibsSegments(j, refGT.nMarkers(), nAnalyses))
                .flatMap(array -> Arrays.stream(array))
                .sorted(HapPairSegment.hapPairComparator())
                .toArray(HapPairSegment[]::new);
        return mergeSortedSegments(segs);
    }

    /**
     * Returns a list of pairwise IBS haplotype segments that have length
     * greater than or equal to {@code this.minCm()} in the specified marker
     * interval
     *
     * @param start the starting marker index
     * @param end the end marker index (exclusive)
     * @param step the difference between consecutive marker indices
     * @return a list of pairwise IBS haplotype segments that have length
     * greater than or equal to {@code this.minCm()}
     *
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || end > this.refGT().nMarkers())}
     */
    private HapPairSegment[] ibsSegments(int start, int end, int step) {
        ArrayList<HapPairSegment> segList = new ArrayList<>();
        PbwtDivUpdater pbwt = new PbwtDivUpdater(refGT.nHaps());
        int[] a = IntStream.range(0, pbwt.nHaps()).toArray();
        int[] d = IntStream.range(0, pbwt.nHaps()).map(i -> start).toArray();
        int maxIbsStart = start-1;
        int endMinusStep = end-step;
        for (int m=start; m<end; m+=step) {
            RefGTRec rec = refGT.get(m);
            pbwt.fwdUpdate(rec, rec.marker().nAlleles(), m, a, d);
            maxIbsStart = updateMaxIbsStart(m, maxIbsStart);
            if (start<=maxIbsStart) {
                if (m<endMinusStep) {
                    addAdjacentIbsSegs(m, step, a, d, maxIbsStart, segList);
                }
                else {
                    lastAddAdjacentIbsSegs(m, a, d, maxIbsStart, segList);
                }
            }
        }
        return segList.toArray(new HapPairSegment[0]);
    }

    private int updateMaxIbsStart(int marker, int previousMaxIbsStart) {
        assert previousMaxIbsStart < marker;
        assert minIbsCm>0f;
        double maxCmPos = cmPos.get(marker) - minIbsCm;
        int candidate = previousMaxIbsStart+1;
        while (cmPos.get(candidate) <= maxCmPos) {
            ++candidate;
        }
        return candidate - 1;
    }

    private void addAdjacentIbsSegs(int m, int step, int[] a, int[] d,
            int maxIbsStart, ArrayList<HapPairSegment> segList) {
        int inclEndPos = refGT.marker(m).pos();
        RefGTRec rec = refGT.get(m+step);
        int a1 = rec.get(a[0]);
        for (int j=1; j<a.length; ++j) {
            int a2 = rec.get(a[j]);
            if (d[j]<=maxIbsStart && a1!=a2) {
                int startPos = refGT.marker(d[j]).pos();
                if (a[j-1] < a[j]) {
                    segList.add(new HapPairSegment(a[j-1], a[j], startPos, inclEndPos));
                }
                else {
                    segList.add(new HapPairSegment(a[j], a[j-1], startPos, inclEndPos));
                }
            }
            a1 = a2;
        }
    }

    private void lastAddAdjacentIbsSegs(int m, int[] a, int[] d,
            int maxIbsStart, ArrayList<HapPairSegment> segList) {
        int inclEndPos = refGT.marker(m).pos();
        for (int j=1; j<a.length; ++j) {
            if (d[j] <= maxIbsStart) {
                int startPos = refGT.marker(d[j]).pos();
                if (a[j-1] < a[j]) {
                    segList.add(new HapPairSegment(a[j-1], a[j], startPos, inclEndPos));
                }
                else {
                    segList.add(new HapPairSegment(a[j], a[j-1], startPos, inclEndPos));
                }
            }
        }
    }

    private HapPairSegment[] mergeSortedSegments(HapPairSegment[] hpsa) {
        int[] ends = IntStream.rangeClosed(1, hpsa.length)
                .parallel()
                .filter(j -> {
                    HapPairSegment prev = hpsa[j-1];
                    return (j==hpsa.length
                            || prev.hap1()!=hpsa[j].hap1()
                            || prev.hap2()!=hpsa[j].hap2()
                            || prev.inclEndPos() < hpsa[j].startPos());
                })
                .toArray();
        return IntStream.range(0, ends.length)
                .parallel()
                .mapToObj(j -> merge(hpsa, ends, j))
                .toArray(HapPairSegment[]::new);
    }

    private HapPairSegment merge(HapPairSegment[] hpsa, int[] ends, int endsIndex) {
        int from = endsIndex==0 ? 0 : ends[endsIndex-1];
        int to = ends[endsIndex];
        HapPairSegment base = hpsa[from];
        int maxInclEndPos = base.inclEndPos();
        for (int j=(from+1); j<to; ++j) {
            if (hpsa[j].inclEndPos() > maxInclEndPos) {
                maxInclEndPos = hpsa[j].inclEndPos();
            }
        }
        return base.inclEndPos()==maxInclEndPos ? base
                : new HapPairSegment(base.hap1(), base.hap2(), base.startPos(), maxInclEndPos);
    }

    /**
     * Returns the analysis parameters.
     * @return the analysis parameters
     */
    public ClustPar par() {
        return par;
    }

    /**
     * Returns the phased genotype data.
     * @return the phased genotype data
     */
    public RefGT refGT() {
        return refGT;
    }

    /**
     * Returns an array of non-decreasing cM positions whose {@code j}-th
     * element is the genetic position of the {@code j}-th marker
     * in {@code this.refGT()}.
     * @return an array of non-decreasing cM positions
     */
    public DoubleArray cmPos() {
        return cmPos;
    }
}
