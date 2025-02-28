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

import blbutil.DoubleArray;
import java.util.Arrays;
import java.util.Random;
import java.util.stream.IntStream;
import vcf.RefGT;

/**
 * <p>Class {@code GlobalIbsProbs} stores sampled one-sided global IBS
 * lengths.  A one-sided IBS length for a pair of distinct
 * haplotypes is the distance in Morgans from a focal position to the first
 * position between the focal position and the distal end of the chromosome
 * for which two haplotypes have discordant alleles or to the marker nearest
 * the distal end of the chromosome if there are no discordant alleles.</p>
 *
 * <p>Instances of class {@code GlobalIbsProbs} are immutable</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class GlobalIbsProbs {

    private final double[] lengths;
    private final double reciprocalSize;

    /**
     * Constructs a new {@code GlobalIbsProbs} from the specified data.  The
     * contract for this method is undefined if {@code morganPos} is not
     * sorted in increasing order or if any element of {@code marganPos} is
     * not a finite number.
     * @param par the analysis parameters
     * @param refGT phased, non-missing genotype data
     * @param morganPos an array containing the Morgan position of each marker
     * @throws IllegalArgumentException if {@code refGT.nHaps() < 2}
     * @throws IllegalArgumentException if
     * {@code refGT.nMarkers() != morganPos().size()}
     * @throws NullPointerException if
     * {@code par == null || refGT == null || morganPos == null}
     */
    public GlobalIbsProbs(ClustPar par, RefGT refGT, DoubleArray morganPos) {
        checkArgs(refGT, morganPos);
        int samplesPerLocus = par.global_segments();
        double[][] lengths0 = IntStream.range(0, par.global_loci())
                .parallel()
                .mapToObj(i -> sampleIbsLengths(refGT, morganPos,
                        samplesPerLocus, (par.seed() + i)))
                .toArray(double[][]::new);

        int index = (int) Math.floor(par.global_quantile()*samplesPerLocus);
        double maxValue = maxValue(lengths0, index, par.global_multiple());
        this.lengths = Arrays.stream(lengths0)
                .parallel()
                .filter(da -> da[index]<=maxValue)
                .flatMapToDouble(da -> Arrays.stream(da))
                .sorted()
                .toArray();
        this.reciprocalSize = 1.0/this.lengths.length;
    }

    private static void checkArgs(RefGT refGT, DoubleArray morganPos) {
        if (refGT.nMarkers()!=morganPos.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        if (refGT.nHaps()<2) {
            throw new IllegalArgumentException(String.valueOf(refGT.nHaps()));
        }
    }

    private static double[] sampleIbsLengths(RefGT refGT, DoubleArray morganPos,
            int samplesPerLocus, long seed) {
        Random rand = new Random(seed);
        double pos = randomGenPos(rand, morganPos);
        double midPos = 0.5*(morganPos.get(0) + morganPos.get(morganPos.size()-1));
        return IntStream.range(0, samplesPerLocus)
                .mapToDouble(i -> sampleIbsLength(refGT, morganPos, midPos, pos, rand))
                .sorted()
                .toArray();
    }

    private static double randomGenPos(Random rand, DoubleArray genPos) {
        double startMorgans = genPos.get(0);
        double endMorgans = genPos.get(genPos.size()-1);
        double pos = startMorgans + rand.nextDouble()*(endMorgans - startMorgans);
        if (pos>=endMorgans) {
            pos = Math.nextDown(pos);
        }
        return pos;
    }

    private static double sampleIbsLength(RefGT refGT, DoubleArray morganPos,
            double midPos, double pos, Random rand) {
        int nHaps = refGT.nHaps();
        int h1 = rand.nextInt(nHaps);
        int h2 = rand.nextInt(nHaps);
        assert nHaps>1;
        while (h1==h2) {
            h2 = rand.nextInt(nHaps);
        }
        if (pos<=midPos) {
            return fwdLength(refGT, morganPos, pos, h1, h2);
        }
        else {
            return bwdLength(refGT, morganPos, pos, h1, h2);
        }
    }

    private static double fwdLength(RefGT refGT, DoubleArray genPos, double pos,
            int h1, int h2) {
        int nMarkersM1 = refGT.nMarkers() - 1;
        int m = genPos.binarySearch(pos);
        if (m<0) {
            m = -m - 1;
        }
        while (m<nMarkersM1 && refGT.allele(m, h1)==refGT.allele(m, h2)) {
            ++m;
        }
        return (genPos.get(m) - pos);
    }

    private static double bwdLength(RefGT refGT, DoubleArray genPos, double pos,
            int h1, int h2) {
        int m = genPos.binarySearch(pos);
        if (m<0) {
            m = -m - 2;
        }
        assert genPos.get(m) <= pos;
        while (m>0 && refGT.allele(m, h1)==refGT.allele(m,h2)) {
            --m;
        }
        return (pos - genPos.get(m));
    }

    private double maxValue(double[][] lengths, int index, double maxMultiple) {
        double[] sortedQuantiles = Arrays.stream(lengths)
                .parallel()
                .mapToDouble(da -> da[index])
                .sorted()
                .toArray();
        int n = lengths.length;
        double median = 0.5*(sortedQuantiles[(n-1)>>1] + sortedQuantiles[n>>1]);
        return maxMultiple*median;
    }

    /**
     * Returns the number of filtered, sampled segments lengths.
     * The filtering excludes loci with unusually long one-sided IBS lengths.
     * @return the number of filtered, sampled segments lengths
     */
    public int nLengths() {
        return lengths.length;
    }

    /**
     * Returns the proportion of filtered, sampled one-sided discord distances
     * that are less than or equal to the specified length. The filtering
     * excludes loci with unusually long one-sided IBS lengths.  Returns
     * {@code 1.0/this.nLengths()} if the the value is less than all filtered,
     * sampled discord distances and returns
     * {@code (this.nLengths() - 1.0)/this.nLengths()} if the the value is
     * greater than or equal to all filtered, sampled discord distances.
     * @param morgans a length in Morgans
     * @return the proportion of filtered, sampled one-sided discord distances
     * that are less than or equal to the specified length
     * @throws IllegalArgumentException if {@code Double.isNaN(morgans) == true}
     */
    public double cdf(double morgans) {
        if (Double.isNaN(morgans)) {
            throw new IllegalArgumentException(String.valueOf(morgans));
        }
        int index = Arrays.binarySearch(lengths, morgans);
        if (index>=0) {
            while ((index+1)<lengths.length && lengths[index]==lengths[index+1]) {
                ++index;
            }
            ++index;
        }
        else {
            index = -index - 1;
        }
        if (index==0) {
            ++index;
        }
        if (index==lengths.length) {
            --index;
        }
        return index*reciprocalSize;
    }
}
