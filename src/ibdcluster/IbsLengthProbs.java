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
import blbutil.FloatList;
import java.util.stream.IntStream;

/**
 * <p>Class {@code IbsLengthProbs} estimates the proportion of haplotype pairs
 * that are IBS on an interval.</p>
 *
 * <p>Instances of class {@code IbsLengthProbs} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbsLengthProbs {

    private final GlobalIbsProbs gip;
    private final float[][] probs;
    private final DoubleArray morgans;

    /**
     * Constructs a new {@code IbsLengthProbs} for the specified data.
     * @param morgans the genetic position for each marker in Morgans
     * @param ibsCnts IBS segment counts
     * @param gip the global IBS length probabilities
     * @throws IllegalArgumentException if
     * {@code  morganPos.size() != ibsCnts.nMarkers()}
     * @throws NullPointerException if any argument is {@code null}
     */
    public IbsLengthProbs(DoubleArray morgans, IbsCounts ibsCnts,
            GlobalIbsProbs gip) {
        if (gip==null) {
            throw new NullPointerException(GlobalIbsProbs.class.toString());
        }
        if (morgans.size() != ibsCnts.nMarkers()) {
            throw new IllegalArgumentException("inconsistent number of markers: "
                    + morgans.size() + " " + ibsCnts.nMarkers());
        }
        int n = ibsCnts.nHaps();
        double invPairsP1 = 1.0/(n*(n-1) + 1.0);
        this.gip = gip;
        this.morgans = morgans;
        this.probs = IntStream.range(0, ibsCnts.nMarkers())
                .parallel()
                .mapToObj(m -> ibsProbs(ibsCnts, m, invPairsP1))
                .toArray(float[][]::new);
    }

    private static float[] ibsProbs(IbsCounts ibsCnts, int start,
            double invPairsP1) {
        int n = ibsCnts.nHaps();
        FloatList probList = new FloatList(1<<8);
        int lastIbsPairs = n*(n-1);
        int end = ibsCnts.end(start);
        for (int m=start; m<end; ++m) {
            int ibsPairs = ibsCnts.counts(start, m);
            probList.add((float) ((lastIbsPairs - ibsPairs + 1)*invPairsP1));
            lastIbsPairs = ibsPairs;
        }
        if (end==ibsCnts.nMarkers()) {
            // store probability of IBS continuing to end of chromosome
            probList.add((float) ((lastIbsPairs + 1)*invPairsP1));
        }
        return probList.toArray();
    }

    /**
     * Returns the number of markers.
     * @return the number of markers
     */
    public int nMarkers() {
        return morgans.size();
    }

    /**
     * Returns the estimated proportion of haplotype pairs that have discordant
     * alleles at the specified end marker index, and that are IBS at all
     * markers between the specified start marker index (inclusive) and the
     * specified end marker index (exclusive). Returns 1.0 if
     * {@code (start == this.nMarkers() && end == this.nMarkers())} because all
     * haplotype pairs are assumed to have discordant alleles at a hypothetical
     * marker with index {@code this.nMarkers()}.
     * @param start the start marker
     * @param end the end marker
     * @return the estimated proportion of haplotype pairs that have discordant
     * alleles at the specified end marker index, and that are IBS at all
     * markers between the specified start marker index (inclusive) and the
     * specified end marker index (exclusive)
     * @throws IndexOutOfBoundsException if
     * {@code (start < 0 || start > this.nMarkers())}
     * @throws IndexOutOfBoundsException if
     * {@code (end > this.nMarkers() || start > end)}
     */
    public double ibsProb(int start, int end) {
        if (start==probs.length && end==probs.length) {
            return 1.0;
        }
        int index = end - start;
        if (index<probs[start].length) {
            assert probs[start][index]>0.0;
            return probs[start][index];
        }
        else if (end==morgans.size()) {
            double length = morgans.get(end-1) - morgans.get(start);
            return 1.0 - gip.cdf(length);
        }
        else {
            double x0 = morgans.get(start);
            double x1 = morgans.get(end-1);
            double x2 = morgans.get(end);
            double p1 = gip.cdf(x1 - x0);
            double p2 = gip.cdf(x2 - x0);
            return (p1==p2) ? 0.5/gip.nLengths() : p2-p1;
        }
    }
}
