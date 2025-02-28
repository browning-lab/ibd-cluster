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
import ints.WrappedIntArray;
import java.util.Arrays;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;
import vcf.GT;
import vcf.GeneticMap;
import vcf.PlinkGenMap;
import vcf.RefGT;
import vcf.ReversedGT;

/**
 * <p>Class {@code ClustData} represents the immutable input data for an
 * ibd-cluster analysis.</p>
 *
 * <p>Instances of class {@code ClustData} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ClustData {

    private static final double MIN_CM_DIST = 1e-6;

    private final ClustPar par;
    private final PlinkGenMap genMap;
    private final String chrom;
    private final RefGT fwdGT;
    private final GT revGT;
    private final WrappedIntArray basePos;
    private final WrappedIntArray reflectedBasePos;
    private final DoubleArray morganPos;
    private final DoubleArray reflectedMorganPos;

    private final IbsLengthProbs fwdIbsProbs;
    private final IbsLengthProbs revIbsProbs;

    /**
     * Constructs a new {@code ClustData} instance from the specified data.
     * The constructor will terminate the Java Virtual Machine with
     * an error message if an I/O error or file format error is detected.
     *
     * @param par the analysis parameters
     * @param refGT the phased genotype data
     * @param genMap the genetic map
     *
     * @throws IllegalArgumentException if a format error is detected in a
     * line of an input data file
     * @throws NumberFormatException if an input integer field is not
     * a parsable integer
     * @throws NullPointerException if
     * {@code (par == null) || (refGT == null) || (genMap == null)}
     */
    public ClustData(ClustPar par, RefGT refGT, PlinkGenMap genMap) {
        this.par = par;
        this.genMap = genMap;
        this.chrom = refGT.marker(0).chrom();
        this.fwdGT = refGT;
        this.revGT = new ReversedGT(fwdGT);

        this.basePos = basePos(fwdGT);
        this.reflectedBasePos = reflect(basePos);
        double[] cmPos = GeneticMap.genPos(genMap, MIN_CM_DIST, fwdGT.markers());
        this.morganPos = morganPos(cmPos);
        this.reflectedMorganPos = reflect(morganPos);

        GlobalIbsProbs gip = new GlobalIbsProbs(par, fwdGT, morganPos);
        IbsCounts fwdIbsCnts = new IbsCounts(par, fwdGT);
        IbsCounts revIbsCnts = fwdIbsCnts.reverse();
        this.fwdIbsProbs = new IbsLengthProbs(morganPos, fwdIbsCnts, gip);
        this.revIbsProbs = new IbsLengthProbs(reflectedMorganPos, revIbsCnts, gip);
    }

    private static WrappedIntArray basePos(GT gt) {
        return new WrappedIntArray(
                IntStream.range(0, gt.nMarkers())
                .parallel()
                .map(j -> gt.marker(j).pos())
        );
    }

    private static DoubleArray morganPos(double[] cmPos) {
        DoubleStream ds = Arrays.stream(cmPos)
                .parallel()
                .map(d -> 0.01*d);
        return new DoubleArray(ds);
    }

    /**
     * Returns the Morgan position of the specified base position.
     * @param morgan a Morgan position
     * @return the base position of the specified Morgan position
     */
    public int morganToBase(double morgan) {
        return ClustUtils.basePos(basePos, morganPos, morgan);
    }

    private static WrappedIntArray reflect(WrappedIntArray ia) {
        int sizeM1 = ia.size()-1;
        return new WrappedIntArray(
                IntStream.range(0, ia.size())
                .parallel()
                .map(j -> -ia.get(sizeM1 - j))
        );
    }

    private static DoubleArray reflect(DoubleArray da) {
        int sizeM1 = da.size()-1;
        return new DoubleArray(
                IntStream.range(0, da.size())
                .parallel()
                .mapToDouble(j -> -da.get(sizeM1 - j))
                .toArray()
        );
    }

    /**
     * Returns the Morgan length of the specified interval.
     * @param startPos the starting base coordinate
     * @param inclEndPos the ending base coordinate (inclusive).
     * @return the Morgan length of the specified interval
     */
    public double morganLength(int startPos, int inclEndPos) {
        return ClustUtils.morganPos(basePos, morganPos, inclEndPos)
                - ClustUtils.morganPos(basePos, morganPos, startPos);
    }

    /**
     * Returns the Morgan position of the specified base position.
     * @param base a base position
     * @return the Morgan position of the specified base position
     */
    public double baseToMorgan(int base) {
        return ClustUtils.morganPos(basePos, morganPos, base);
    }

    /**
     * Returns the analysis parameters.
     * @return the analysis parameters
     */
    public ClustPar par() {
        return par;
    }

    /**
     * Returns the genetic map.
     * @return the genMap
     */
    public GeneticMap genMap() {
        return genMap;
    }

    /**
     * Returns the chromosome identifier for the markers in {@code this.fwdGT()}
     * @return the chromosome identifier
     */
    public String chrom() {
        return chrom;
    }

    /**
     * Returns the input phased genotype data.
     * @return the input phased genotype data
     */
    public RefGT fwdGT() {
        return fwdGT;
    }

    /**
     * Returns the input phased genotype data with markers in reverse order.
     * @return the input phased genotype data with markers in reverse order
     */
    public GT revGT() {
        return revGT;
    }

    /**
     * Returns the list of marker base positions.
     * @return the list of marker base positions
     */
    public WrappedIntArray basePos() {
        return basePos;
    }

    /**
     * Returns the list of reflected marker base positions.
     *
     * @return the list of reflected marker base positions
     */
    public WrappedIntArray reflectedBasePos() {
        return reflectedBasePos;
    }

    /**
     * Returns the list of marker Morgan positions.
     * @return the list of marker Morgan positions
     */
    public DoubleArray morganPos() {
        return morganPos;
    }

    /**
     * Returns the list of reflected marker Morgan positions.
     * @return the list of reflected marker Morgan positions
     */
    public DoubleArray reflectedMorganPos() {
        return reflectedMorganPos;
    }

    /**
     * Returns the one-sided forward IBS length probabilities.
     * @return the one-sided forward IBS length probabilities
     */
    public IbsLengthProbs fwdIbsProbs() {
        return fwdIbsProbs;
    }

    /**
     * Returns the one-sided reverse IBS length probabilities.
     * @return the one-sided reverse IBS length probabilities
     */
    public IbsLengthProbs revIbsProbs() {
        return revIbsProbs;
    }
}
