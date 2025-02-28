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

import blbutil.Utilities;
import ints.CharArray;
import ints.IntArray;
import ints.IntList;
import ints.UnsignedByteArray;
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;
import vcf.Markers;
import vcf.RefGT;

/**
 * <p>Class {@code IbsCounts} stores the number of ordered haplotype pairs
 * in a subset of haplotypes that are identical by state (IBS) on a subset of
 * marker intervals.</p>
 *
 * <p>Instances of class {@code IbsCounts} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbsCounts {

    private final int nHaps;
    private final int[][] counts;

    /**
     * Constructs a new {@code IbsCounts} instance from the specified data.
     * @param par the analysis parameters
     * @param refGT the phased, non-missing genotype data
     * @throws IllegalArgumentException if {@code refGT.nHaps() < 2}
     * @throws NullPointerException if {@code (par == null || refGT == null)}
     */
    public IbsCounts(ClustPar par, RefGT refGT) {
        if (refGT.nHaps()<2) {
            throw new IllegalArgumentException(String.valueOf(refGT.nHaps()));
        }
        int[] hapList = hapList(par, refGT);
        int n = hapList.length;
        assert (n<=ClustPar.max_local_segments()) && (((long) n)*(n-1)<Integer.MAX_VALUE);
        IntArray[] sampleAlleles = alleles(refGT, hapList);
        Boolean[] isMonomorphic = isMonomorphic(sampleAlleles);
        int maxNAlleles = maxNAlleles(refGT);

        float maxCDF = par.local_max_cdf();
        int minIbsPairs = (int) Math.ceil(((1.0-maxCDF)*n)*(n-1));
        this.nHaps = hapList.length;
        this.counts = IntStream.range(0, refGT.nMarkers())
                .parallel()
                .mapToObj(m -> counts(refGT.markers(), maxNAlleles, m,
                        sampleAlleles, isMonomorphic, minIbsPairs))
                .toArray(int[][]::new);
    }

    /* This constructor reverses markers in the specified {@code IbsCounts) instance */
    private IbsCounts(IbsCounts fwdCnts) {
        int nMarkers = fwdCnts.nMarkers();
        this.nHaps = fwdCnts.nHaps;
        this.counts = IntStream.range(0, nMarkers)
                .parallel()
                .mapToObj(m -> revCounts(fwdCnts, m))
                .toArray(int[][]::new);
    }

    private static int[] revCounts(IbsCounts fwdCnts, int revStart) {
        IntList revCnts = new IntList(1<<8);
        int inclEnd = fwdCnts.nMarkers() - 1 - revStart;
        for (int start=inclEnd; start>=0 && inclEnd<fwdCnts.end(start); --start) {
            revCnts.add(fwdCnts.counts(start, inclEnd));
        }
        return revCnts.toArray();
    }

    private static int[] hapList(ClustPar par, RefGT refGT) {
        int[] allHaps = IntStream.range(0, refGT.nHaps()).parallel().toArray();
        int maxLocalHaps = par.local_segments();
        if (refGT.nHaps()<=maxLocalHaps) {
            return allHaps;
        }
        else {
            Utilities.shuffle(allHaps, maxLocalHaps, new Random(par.seed()));
            Arrays.sort(allHaps, 0, maxLocalHaps);
            return Arrays.copyOf(allHaps, maxLocalHaps);
        }
    }

    private static IntArray[] alleles(RefGT gt, int[] hapList) {
        return IntStream.range(0, gt.nMarkers())
                .parallel()
                .mapToObj(m -> alleles(gt, m, hapList))
                .toArray(IntArray[]::new);
    }

    private static IntArray alleles(RefGT gt, int marker, int[] hapList) {
        int[] alleles = new int[hapList.length];
        for (int j=0; j<alleles.length; ++j) {
            alleles[j] = gt.allele(marker, hapList[j]);
        }
        int nAlleles = gt.marker(marker).nAlleles();
        if (nAlleles<=256) {
            return new UnsignedByteArray(alleles);
        }
        else if (nAlleles<=65536) {
            return new CharArray(alleles);
        }
        else {
            return new WrappedIntArray(alleles);
        }
    }

    private static Boolean[] isMonomorphic(IntArray[] alleles) {
        return Arrays.stream(alleles)
                .parallel()
                .map(ia -> isMonomorphic(ia))
                .toArray(Boolean[]::new);
    }

    private static Boolean isMonomorphic(IntArray alleles) {
        for (int j=1, n=alleles.size(); j<n; ++j) {
            if (alleles.get(j)!=alleles.get(j-1)) {
                return false;
            }
        }
        return true;
    }

    private static int maxNAlleles(RefGT gt) {
        return IntStream.range(0, gt.nMarkers())
                .parallel()
                .map(m -> gt.marker(m).nAlleles())
                .max()
                .orElse(0);
    }

    private static int[] counts(Markers markers, int maxNAlleles, int start,
            IntArray[] alleles, Boolean[] isMonomorphic, int minIbsPairs) {
        int nMarkers = alleles.length;
        int nHaps = alleles[start].size();
        IntList cnts = new IntList(1<<8);
        int[] hap2Seq = new int[nHaps];
        int[] seq2Cnt = new int[nHaps];
        int[] seqAlMap = new int[maxNAlleles*nHaps];
        seq2Cnt[0] = nHaps;
        int nSeq = 1;
        int ibsPairs = nHaps*(nHaps-1);
        for (int m=start; m<nMarkers && ibsPairs>=minIbsPairs; ++m) {
            if (isMonomorphic[m]) {
                cnts.add(ibsPairs);
            }
            else {
                int nAlleles = markers.marker(m).nAlleles();
                Arrays.fill(seqAlMap, 0, nAlleles*nSeq, -1);
                Arrays.fill(seq2Cnt, 0, nSeq, 0);
                nSeq = 0;
                for (int j=0; j<nHaps; ++j) {
                    int seqAlIndex = hap2Seq[j]*nAlleles + alleles[m].get(j);
                    int seqIndex = seqAlMap[seqAlIndex];
                    if (seqIndex<0) {
                        seqIndex = nSeq++;
                        seqAlMap[seqAlIndex] = seqIndex;
                    }
                    hap2Seq[j] = seqIndex;
                    ++seq2Cnt[seqIndex];
                }
                ibsPairs = sumIbsPairs(seq2Cnt, nSeq);
                if (ibsPairs>=minIbsPairs) {
                    cnts.add(ibsPairs);
                }
            }
        }
        return cnts.toArray();
    }

    private static int sumIbsPairs(int[] seqCnts, int nSeq) {
        int sum = 0;
        for (int j=0; j<nSeq; ++j) {
            sum += seqCnts[j]*(seqCnts[j]-1);
        }
        return sum;
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return counts.length;
    }

    /**
     * Returns the number of random haplotypes used to generate IBS segment
     * counts.
     * @return the number of random haplotypes used to generate IBS segment
     * counts
     */
    public int nHaps() {
        return nHaps;
    }

    /**
     * Returns the number of ordered haplotype pairs from a random set of
     * {@code this.nHaps()} haplotypes that are identical by state on the
     * interval from the specified start marker index (inclusive) to the
     * specified end marker index (inclusive).
     * @param start the start marker (inclusive)
     * @param inclEnd the end marker (inclusive)
     * @return the number of ordered haplotype pairs from a random set of
     * {@code this.nHaps()} haplotypes that are IBS on the interval
     * from the specified start marker index (inclusive) to the specified
     * end marker index (inclusive)
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start >= this.nMarkers()}
     * @throws IndexOutOfBoundsException if
     * {@code inclEnd < start || inclEnd >= this.end(start)}
     */
    public int counts(int start, int inclEnd) {
        return counts[start][inclEnd - start];
    }

    /**
     * Returns the exclusive end index for the {@code inclEnd} parameter
     * of the {@code this.counts()} method.
     * @param start a marker index
     * @return the exclusive end index for the {@code inclEnd} parameter
     * of the {@code this.counts()} method
     * @throws IndexOutOfBoundsException if
     * {@code start < 0 || start >= this.nMarkers()}
     */
    public int end(int start) {
        return start + counts[start].length;
    }

    /**
     * Returns an {@code IbsCounts} instance obtained by reversing the marker
     * order of {@code this}.
     * @return an {@code IbsCounts} instance obtained by reversing the marker
     * order of {@code this}
     */
    public IbsCounts reverse() {
        return new IbsCounts(this);
    }
}
