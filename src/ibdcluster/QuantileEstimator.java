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

import blbutil.Const;
import blbutil.DoubleArray;
import ints.WrappedIntArray;
import java.util.Arrays;
import vcf.GT;

/**
 * <p>Class {@code QuantileEstimator} estimates a quantile of an IBD segment
 * end point distribution.</p>
 *
 * <p>Instances of class {@code QuantileEstimator} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class QuantileEstimator {

    private final ClustData data;
    private final int nMarkers;
    private final Data fwdData;
    private final Data revData;
    private final double ne;
    private final double err;
    private final double gc_err;
    private final int gc_bp;
    private final float minCdfRatio;

    // scratch space for storing CDF
    private final double[] cdf;
    private int cdfStart;
    private int cdfEnd;

    /**
     * Constructs a new {@code QuantileEstimator} instance from the
     * specified data.
     * @param data the analysis input data
     * @throws NullPointerException if {@code data == null}
     */
    public QuantileEstimator(ClustData data) {
        ClustPar par = data.par();
        this.data = data;
        this.nMarkers = data.fwdGT().nMarkers();
        this.fwdData = new Data(data, true);
        this.revData = new Data(data, false);
        this.ne = par.ne();
        this.err = par.discord();
        this.gc_bp = par.gc_bases();
        this.gc_err = par.gc_discord();
        this.minCdfRatio = par.min_cdf_ratio();
        this.cdf = new double[data.fwdGT().nMarkers()+1]; // includes hypothetical ending marker
    }

    /**
     * Returns the input data.
     * @return the input data
     */
    public ClustData data() {
        return data;
    }

    private double baseToMorgans(int basePos) {
        return ClustUtils.morganPos(data.basePos(), data.morganPos(), basePos);
    }

    /**
     * Estimates and returns the specified quantile of the ending point
     * of the specified IBD segment.
     * @param hap1 the first IBD haplotype
     * @param hap2 the second IBD haplotype
     * @param startMorgans an estimated genetic position in Morgans of the
     * start of the IBD segment
     * @param focusPos the base pair coordinate of the focus position in
     * the IBD segment
     * @param prob a probability
     * @return the end-point quantile corresponding to the specified probability
     * @throws IllegalArgumentException if
     * {@code Double.isNaN(ibdStartMorgans) == true}
     * @throws IllegalArgumentException if the focus position is less than or
     * equal to the estimated base pair position of the specified start of the
     * IBD segment
     * @throws IllegalArgumentException if
     * {@code (p <= 0.0 || p >= 1.0 || Double.isNaN(p))}
     * @throws IndexOutOfBoundsException if
     * {@code hap1 < 0 || hap1 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2 < 0 || hap2 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code (focusPos < data.fwdGT().marker(0).pos())}
     */
    int fwdQuantile(int hap1, int hap2, double startMorgans,
            int focusPos, double prob) {
        double focusMorgans = baseToMorgans(focusPos);
        setCDF(fwdData, hap1, hap2, startMorgans, focusPos, focusMorgans);
        return quantile(fwdData, startMorgans, focusPos, focusMorgans, prob);
    }

    /**
     * Estimates and returns the specified quantile of the ending point
     * of the specified IBD segment.
     * @param hap1 the first IBD haplotype
     * @param hap2 the second IBD haplotype
     * @param startMorgans an estimated genetic position in Morgans of the
     * start of the IBD segment
     * @param focusPos the base pair coordinate of the focus position in
     * the IBD segment
     * @param prob a probability
     * @return the end-point quantile corresponding to the specified probability
     * @throws IllegalArgumentException if
     * {@code Double.isNaN(ibdStartMorgans) == true}
     * @throws IllegalArgumentException if the focus position is less than or
     * equal to the estimated base pair position of the specified start of the
     * IBD segment
     * @throws IllegalArgumentException if
     * {@code (prob <= 0.0 || prob >= 1.0 || Double.isNaN(prob))}
     * @throws IllegalArgumentException if
     * {@code (Double.isFinite(trimMorgans)==false || trimMorgans < 0.0d)}
     * @throws IndexOutOfBoundsException if
     * {@code hap1 < 0 || hap1 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2 < 0 || hap2 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code (focusPos < data.fwdGT().marker(0).pos())}
     */
    double fwdMorganQuantile(int hap1, int hap2, double startMorgans,
            int focusPos, double prob) {
        double focusMorgans = baseToMorgans(focusPos);
        setCDF(fwdData, hap1, hap2, startMorgans, focusPos, focusMorgans);
        return morganQuantile(fwdData, startMorgans, focusMorgans, prob);
    }

    /**
     * Estimates and returns the specified quantile of the ending point
     * of the specified IBD segment.
     * @param hap1 the first IBD haplotype
     * @param hap2 the second IBD haplotype
     * @param inclEndMorgans an estimated genetic position in Morgans of the
     * end of the IBD segment
     * @param focusPos the base pair coordinate of the focus position in
     * the IBD segment
     * @param prob a probability
     * @param trimMorgans the trim length in Morgans
     * @return the end-point quantile corresponding to the specified probability
     * @throws IllegalArgumentException if
     * {@code Double.isNaN(ibdEndMorgans) == true}
     * @throws IllegalArgumentException if the focus position is greater than or
     * equal to the estimated base pair position of the specified end of the IBD
     * segment
     * @throws IllegalArgumentException if
     * {@code (prob <= 0.0 || prob >= 1.0 || Double.isNaN(prob))}
     * @throws IllegalArgumentException if
     * {@code (Double.isFinite(trimMorgans)==false || trimMorgans < 0.0d)}
     * @throws IndexOutOfBoundsException if
     * {@code hap1 < 0 || hap1 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2 < 0 || hap2 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code (focusPos > data.fwdGT().marker(data.fwdGT().nMarkers()-1).pos())}
     */
    int bwdQuantile(int hap1, int hap2, int focusPos,
            double inclEndMorgans, double prob) {
        double focusMorgans = baseToMorgans(focusPos);
        setCDF(revData, hap1, hap2, -inclEndMorgans, -focusPos, -focusMorgans);
        return -quantile(revData, -inclEndMorgans, -focusPos, -focusMorgans, prob);
    }

    /**
     * Estimates and returns the Morgan position of the specified quantile
     * of the end point distribution of the specified IBD segment.
     * @param hap1 the first IBD haplotype
     * @param hap2 the second IBD haplotype
     * @param inclEndMorgans an estimated genetic position in Morgans of the
     * end of the IBD segment
     * @param focusPos the base pair coordinate of the focus position in
     * the IBD segment
     * @param prob a probability
     * @return the end-point quantile corresponding to the specified probability
     * @throws IllegalArgumentException if
     * {@code Double.isNaN(ibdEndMorgans) == true}
     * @throws IllegalArgumentException if the focus position is greater than or
     * equal to the estimated base pair position of the specified end of the IBD
     * segment
     * @throws IllegalArgumentException if
     * {@code (prob <= 0.0 || prob >= 1.0 || Double.isNaN(prob))}
     * @throws IndexOutOfBoundsException if
     * {@code hap1 < 0 || hap1 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code hap2 < 0 || hap2 >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code (focusPos > data.fwdGT().marker(data.fwdGT().nMarkers()-1).pos())}
     */
    double bwdMorganQuantile(int hap1, int hap2, int focusPos,
            double inclEndMorgans, double prob) {
        double focusMorgans = baseToMorgans(focusPos);
        setCDF(revData, hap1, hap2, -inclEndMorgans, -focusPos, -focusMorgans);
        return -morganQuantile(revData, -inclEndMorgans, -focusMorgans, prob);
    }

    /*
     * Stores the prior probability distribution of the position of the end
     * of an IBD segment. For markers {@code m} in the interval bounded by
     * {@code this.cdfStart+1} (inclusive) and {@code this.cdfEnd} (exclusive),
     * the probability that the IBD segment end is between markers
     * {@code (m-1)} and {@code m} is stored in the {@code m}-th element
     * of the specified {@code cdf} array. The probability that an
     * IBD segment end is between the focus position and {@code this.cdfStart}
     * is stored in the {@code this.cdf[cdfStart]}.
     */
    private void setCDF(Data data, int h1, int h2, double startMorgans,
            int focusPos, double focusMorgans) {
        this.cdfStart = data.nextMarker(focusPos);
        cdf[cdfStart-1] = 0.0;
        double factor = 1.0;
        double F1 = ClustUtils.F((focusMorgans - startMorgans), ne);
        int start = cdfStart;
        int nextDiscord = data.nextDiscord(h1, h2, start);
        int minNextDiscordPos = data.pos(nextDiscord) + gc_bp;
        while (true) {
            cdfEnd = nextDiscord+1;
            for (int m=start; m<cdfEnd; ++m) {
                double F2 = ClustUtils.F((data.morgans(m) - startMorgans), ne);
                cdf[m] = cdf[m-1] + (F2-F1)*data.ibsProb(m, nextDiscord)*factor;
                F1 = F2;
            }
            if (finished(start)) {
                scale(cdf, cdfStart, cdfEnd, (1.0/cdf[cdfEnd-1]));
                return;
            }
            if (cdf[cdfEnd-1]>1e50) {
                double scaleFactor = (1.0/cdf[cdfEnd-1]);
                scale(cdf, cdfStart, cdfEnd, scaleFactor);
                factor *= scaleFactor;
            }
            start = cdfEnd;
            nextDiscord = data.nextDiscord(h1, h2, start);
            int discordPos = data.pos(nextDiscord);
            double num = gc_err;
            if (discordPos >= minNextDiscordPos) {
                num = err;
                minNextDiscordPos = discordPos + gc_bp;
            }
            factor *= (num/data.ibsProb(start, nextDiscord));
        }
    }

    private boolean finished(int lastEnd) {
        if (cdfEnd==cdf.length) {
            return true;
        }
        double lastValue = cdf[cdfEnd-1];
        return (lastValue - cdf[lastEnd-1])<(minCdfRatio*lastValue);
    }

    private static void scale(double[] da, int start, int end, double factor) {
        for (int j=start; j<end; ++j) {
            da[j] *= factor;
        }
    }

    private int quantile(Data data, double startMorgans, int focusPos,
            double focusMorgans, double p) {
        if (p<=0d || p>=1d || Double.isNaN(p)==true) {
            throw new IllegalArgumentException(String.valueOf(p));
        }
        int index = Arrays.binarySearch(cdf, cdfStart, cdfEnd, p);
        if (index<0) {
            index = -index - 1;
        }
        if (index==nMarkers) {
            return data.pos(nMarkers);
        }
        double p1 = cdf[index-1];
        double p2 = cdf[index];
        assert p1<=p && p<=p2;

        double x1 = (index==cdfStart) ? focusMorgans : data.morgans(index-1);
        double x2 = data.morgans(index);

        double F1 = ClustUtils.F((x1 - startMorgans), ne);
        double F2 = ClustUtils.F((x2 - startMorgans), ne);
        double pp = F1 + ((p-p1)/(p2-p1))*(F2-F1);
        double x = startMorgans + ClustUtils.invF(pp, ne);
        assert x1<=x2;
        assert (x<0  || ((0.999999999*x1)<=x && x<=(1.000000001*x2)));
        assert (x>=0 || ((1.000000001*x1)<=x && x<=(0.999999999*x2)));
        double delta = (x-x1)/(x2-x1);
        if (delta<0.0) {
            delta = Math.nextUp(0.0);
        }
        if (delta>1.0) {
            delta = Math.nextDown(1.0);
        }

        // minimum quantile needs to be focusPos+1 to avoid division by 0
        int y1 = index==cdfStart ? focusPos+1 : data.pos(index-1);
        int y2 = data.pos(index);
        int y = (int) Math.rint(y1 + delta*(y2-y1));
        assert y1<=y && y<=y2;
        return y;
    }

    private double morganQuantile(Data data, double startMorgans,
            double focusMorgans,  double p) {
        if (p<=0d || p>=1d || Double.isNaN(p)==true) {
            throw new IllegalArgumentException(String.valueOf(p));
        }
        int index = Arrays.binarySearch(cdf, cdfStart, cdfEnd, p);
        if (index<0) {
            index = -index - 1;
        }
        if (index==nMarkers) {
            return this.data.baseToMorgan(data.pos(nMarkers));
        }
        double p1 = cdf[index-1];
        double p2 = cdf[index];
        assert p1<=p && p<=p2;

        double x1 = (index==cdfStart) ? focusMorgans : data.morgans(index-1);
        double x2 = data.morgans(index);

        double F1 = ClustUtils.F((x1 - startMorgans), ne);
        double F2 = ClustUtils.F((x2 - startMorgans), ne);
        double pp = F1 + ((p-p1)/(p2-p1))*(F2-F1);
        return startMorgans + ClustUtils.invF(pp, ne);
    }

    private static class Data {

        private final GT gt;
        private final int nMarkers;
        private final int extBasePos;
        private final WrappedIntArray basePos;
        private final DoubleArray morganPos;
        private final IbsLengthProbs ibsProbs;

        // morgan position of hypothetical discordant marker with index morgans.size()
        private final double extMorganPos;

        private Data(ClustData data, boolean fwd) {
            double endMorgans = data.par().end_morgans();
            int lastIndex = data.fwdGT().nMarkers() - 1;
            if (fwd) {
                assert data.fwdGT().nMarkers()==data.basePos().size();
                assert data.fwdGT().nMarkers()==data.morganPos().size();
                this.gt = data.fwdGT();
                this.nMarkers = data.fwdGT().nMarkers();
                this.basePos = data.basePos();
                this.morganPos = data.morganPos();
                this.extBasePos = data.basePos().get(lastIndex);
                this.extMorganPos = data.morganPos().get(lastIndex) + endMorgans;
                this.ibsProbs = data.fwdIbsProbs();
            }
            else {
                assert data.revGT().nMarkers()==data.basePos().size();
                assert data.revGT().nMarkers()==data.morganPos().size();
                this.gt = data.revGT();
                this.nMarkers = data.revGT().nMarkers();
                this.basePos = data.reflectedBasePos();
                this.morganPos = data.reflectedMorganPos();
                this.extBasePos = data.reflectedBasePos().get(lastIndex);
                this.extMorganPos = data.reflectedMorganPos().get(lastIndex) + endMorgans;
                this.ibsProbs = data.revIbsProbs();
            }
        }

        private double morgans(int marker) {
            return marker==nMarkers ? extMorganPos : morganPos.get(marker);
        }

        private int pos(int marker) {
            return marker==nMarkers ? extBasePos: basePos.get(marker);
        }

        private int nextMarker(int position) {
            int insPt = basePos.binarySearch(position);
            return insPt<0 ? -insPt-1 : insPt+1;
        }

        private int nextDiscord(int hap1, int hap2, int start) {
            int end = gt.nMarkers();
            int m = start;
            while (m<end && gt.allele(m, hap1)==gt.allele(m, hap2)) {
                ++m;
            }
            return m;
        }

        private double ibsProb(int marker, int nextDiscord) {
            return ibsProbs.ibsProb(marker, nextDiscord);
        }

        private int morganToBase(double morgan) {
            return ClustUtils.basePos(basePos, morganPos, morgan);
        }
    }
}
