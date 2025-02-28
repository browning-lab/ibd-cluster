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
import java.util.Random;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * <p>Class {@code IbdEstimator} estimates IBD segment end-points.</p>
 *
 * <p>Instances of class {@code IbdEstimator} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class IbdEstimator {

    private final ClustData data;
    private final int chromStartPos;
    private final int chromInclEndPos;
    private final float prefocusQuantile;
    private final float quantile;
    private final QuantileEstimator quantEst;
    private final double minIbdMorgans;
    private final float trimMorgans;
    private final int maxIts;
    private final int maxItsM2;
    private final double maxRelChange;
    private final boolean fixFocus;

    private int h1;
    private int h2;
    private int startPos;
    private int inclEndPos;
    private int focusPos;
    private double startMorgans;
    private double inclEndMorgans;
    private double focusMorgans;

//    // xxx
//    private static final AtomicInteger index = new AtomicInteger(0);
//    private final Random random = new Random(index.getAndIncrement());
//    // xx

    /**
     * Constructs an {@code IbdEstimator} instance for the specified data.
     * @param data the input data
     * @param stats an object which stores analysis statistics
     * @throws NullPointerException if {@code (data == null) || (stats == null)}
     */
    public IbdEstimator(ClustData data, ClustStats stats) {
        if (stats==null) {
            throw new NullPointerException(ClustStats.class.toString());
        }
        ClustPar par = data.par();
        this.data = data;
        this.chromStartPos = data.fwdGT().marker(0).pos();
        this.chromInclEndPos = data.revGT().marker(0).pos();
        this.prefocusQuantile = par.prefocus_quantile();
        this.quantile = par.quantile();
        this.trimMorgans = 0.01f*par.trim();
        this.quantEst = new QuantileEstimator(data);
        this.maxIts = par.max_its()<<1; // doubled since there are two ends
        this.maxItsM2 = this.maxIts - 2;
        this.maxRelChange = par.max_rel_change();
        this.fixFocus = par.fix_focus();
        this.minIbdMorgans = 0.01*data.par().min_ibd_cm();
    }

    /**
     * Returns the input data.
     * @return the input data
     */
    public ClustData data() {
        return data;
    }

    /**
     * Returns an estimated IBD segment whose endpoints are estimated
     * from a focal point within the specified haplotype pair segment.
     * @param ibsSegment an haplotype pair segment
     * @return the estimated IBD segment
     *
     * @throws IllegalArgumentException if
     * {@code hps.startPos() < this.data().fwdGT().marker(0)}
     * @throws IllegalArgumentException if
     * {@code hps.inclEndPos() > this.data().revGT().marker(0)}
     * @throws IndexOutOfBoundsException if
     * {@code hps.hap1() < 0 || hps.hap1() >= this.data().fwdGT().nHaps()}
     * @throws IndexOutOfBoundsException if
     * {@code hps.hap2() < 0 || hps.hap2() >= this.data().fwdGT().nHaps()}
     * @throws NullPointerException if {@code hps == null}
     */
    public HapPairSegment ibdSegment(HapPairSegment ibsSegment) {
        checkSegment(ibsSegment);
        initializeFields(ibsSegment);
        int noUpdateCnt = 0;
        for (int j=0; noUpdateCnt<2 && j<maxItsM2; ++j) {
            if ((j & 1)==1) {
                int newStartPos = quantEst.bwdQuantile(h1, h2, focusPos,
                        inclEndMorgans, prefocusQuantile);
                boolean updated = updateStartPos(ibsSegment, newStartPos);
                noUpdateCnt = updated ? 0 : noUpdateCnt + 1;
            }
            else {
                int newInclEndPos = quantEst.fwdQuantile(h1, h2, startMorgans,
                        focusPos, prefocusQuantile);
                boolean updated = updateInclEndPos(ibsSegment, newInclEndPos);
                noUpdateCnt = updated ? 0 : noUpdateCnt + 1;
            }
        }
//        // xxx
//        if (random.nextDouble() < 0.0001) {
//            writeDebugInfo(ibsSegment);
//        }
//        // xx
        return trimmedIbdSegment(ibsSegment);
    }

    private HapPairSegment trimmedIbdSegment(HapPairSegment ibsSegment) {
        double ibdStartMorgans = quantEst.bwdMorganQuantile(h1, h2,
                focusPos, inclEndMorgans, quantile);
        double ibdEndMorgans = quantEst.fwdMorganQuantile(h1, h2,
                startMorgans, focusPos, quantile);
        double ibdLength = ibdEndMorgans - ibdStartMorgans;
        double trimmedStartMorgans = ibdStartMorgans + trimMorgans;
        double trimmedEndMorgans = ibdEndMorgans - trimMorgans;
        if ((ibdLength >= minIbdMorgans)
                && (trimmedStartMorgans <= trimmedEndMorgans)) {
            int ibdStartPos = data.morganToBase(trimmedStartMorgans);
            int ibdInclEndPos = data.morganToBase(trimmedEndMorgans);
            return new HapPairSegment(ibsSegment.hap1(), ibsSegment.hap2(),
                    ibdStartPos, ibdInclEndPos);
        }
        else {
            return HapPairSegment.ZERO_LENGTH_SEGMENT;
        }
    }

    // xxx
    private void writeDebugInfo(HapPairSegment ibsSegment) {
        double ibdStartMorgans = quantEst.bwdMorganQuantile(h1, h2,
                focusPos, inclEndMorgans, quantile);
        double ibdEndMorgans = quantEst.fwdMorganQuantile(h1, h2,
                startMorgans, focusPos, quantile);
        double ibdLength = ibdEndMorgans - ibdStartMorgans;

        int ibdStartPos = data.morganToBase(ibdStartMorgans);
        int ibdInclEndPos = data.morganToBase(ibdEndMorgans);
        HapPairSegment ibdSegment = new HapPairSegment(ibsSegment.hap1(),
                ibsSegment.hap2(), ibdStartPos, ibdInclEndPos);

        double trimmedStartMorgans = ibdStartMorgans + trimMorgans;
        double trimmedEndMorgans = ibdEndMorgans - trimMorgans;
        int trimmedIbdStartPos = data.morganToBase(trimmedStartMorgans);
        int trimmedIbdInclEndPos = data.morganToBase(trimmedEndMorgans);

        boolean failedLengthFilter = (ibdLength < minIbdMorgans)
                || (trimmedStartMorgans > trimmedEndMorgans);

        HapPairSegment trimmedIbd = HapPairSegment.ZERO_LENGTH_SEGMENT;
        if (trimmedIbdStartPos <= trimmedIbdInclEndPos) {
            trimmedIbd = new HapPairSegment(ibsSegment.hap1(), ibsSegment.hap2(),
                    trimmedIbdStartPos, trimmedIbdInclEndPos);
        }

        StringBuilder sb = new StringBuilder(128);
        sb.append(Const.nl);
        writeLine(sb, "IE 120: IBS segment:         ", ibsSegment, false);
        writeLine(sb, "IE 122: IBD segment:         ", ibdSegment, failedLengthFilter);
        writeLine(sb, "IE 124: trimmed IBD Segment: ", trimmedIbd, failedLengthFilter);
        System.out.println(sb);
    }

    private void writeLine(StringBuilder sb, String message, HapPairSegment hps,
            boolean failedLengthFilter) {
        sb.append(message);
        if (hps==HapPairSegment.ZERO_LENGTH_SEGMENT) {
            sb.append("zero-length-segment");
        }
        else {
            sb.append(hps.toString());
            if (failedLengthFilter) {
                sb.append(Const.tab);
                sb.append("FAIL-MIN_IBD_CM");
            }
        }
        sb.append(Const.nl);
    }
    // xx

    private void checkSegment(HapPairSegment hps) {
        if ((hps.startPos() < chromStartPos)
                || (hps.inclEndPos() > chromInclEndPos)) {
            StringBuilder sb = new StringBuilder();
            sb.append(Const.nl);
            sb.append("Error: haplotype segment extends beyond input markers");
            sb.append(Const.nl);
            sb.append("Marker interval:   ");
            sb.append(data.chrom());
            sb.append(Const.colon);
            sb.append(data.fwdGT().marker(0).pos());
            sb.append('-');
            sb.append(data.revGT().marker(0).pos());
            sb.append(Const.nl);
            sb.append("Haplotype segment: ");
            sb.append(data.chrom());
            sb.append(Const.colon);
            sb.append(hps.startPos());
            sb.append('-');
            sb.append(hps.inclEndPos());
            throw new IllegalArgumentException(sb.toString());
        }
    }

    private void initializeFields(HapPairSegment hps) {
        h1 = hps.hap1();
        h2 = hps.hap2();
        startPos = hps.startPos();
        inclEndPos = hps.inclEndPos();
        focusPos = (hps.startPos() + hps.inclEndPos()) >>> 1;
        startMorgans = data.baseToMorgan(hps.startPos());
        inclEndMorgans = data.baseToMorgan(hps.inclEndPos());
        focusMorgans = data.baseToMorgan(focusPos);
    }

    private boolean updateInclEndPos(HapPairSegment ibs, int newInclEndPos) {
        double newInclEndMorgans = data.baseToMorgan(newInclEndPos);
        boolean allowEndUpdate = allowEndUpdate(focusMorgans, inclEndMorgans, newInclEndMorgans);
        if (allowEndUpdate) {
            int newFocusPos = focusPos;
            double newFocusMorgans = focusMorgans;
            if (fixFocus==false) {
                newFocusPos = (startPos + newInclEndPos) >>> 1;
                if (newFocusPos <= ibs.startPos()) {
                    newFocusPos = ibs.startPos() + 1;
                }
                if (newFocusPos >= ibs.inclEndPos()) {
                    newFocusPos = ibs.inclEndPos() - 1;
                }
                newFocusMorgans = data.baseToMorgan(newFocusPos);
            }
            if ((newInclEndMorgans - newFocusMorgans)>0
                    && (newFocusMorgans - startMorgans)>0) {
                focusPos = newFocusPos;
                focusMorgans = newFocusMorgans;
                inclEndPos = newInclEndPos;
                inclEndMorgans = newInclEndMorgans;
                return true;
            }
            else {
                return false;
            }
        }
        return allowEndUpdate;
    }

    private boolean updateStartPos(HapPairSegment ibs, int newStartPos) {
        double newStartMorgans = data.baseToMorgan(newStartPos);
        boolean allowEndUpdate = allowEndUpdate(focusMorgans, startMorgans, newStartMorgans);
        if (allowEndUpdate) {
            int newFocusPos = focusPos;
            double newFocusMorgans = focusMorgans;
            if (fixFocus==false) {
                newFocusPos = (newStartPos + inclEndPos) >>> 1;
                if (newFocusPos <= ibs.startPos()) {
                    newFocusPos = ibs.startPos() + 1;
                }
                if (newFocusPos >= ibs.inclEndPos()) {
                    newFocusPos = ibs.inclEndPos() - 1;
                }
                newFocusMorgans = data.baseToMorgan(newFocusPos);
            }
            if ((newFocusMorgans - newStartMorgans)>0
                    && (inclEndMorgans - newFocusMorgans)>0) {
                startPos = newStartPos;
                startMorgans = newStartMorgans;
                focusPos = newFocusPos;
                focusMorgans = newFocusMorgans;
                return true;
            }
            else {
                return false;
            }
        }
        return allowEndUpdate;
    }

    private boolean allowEndUpdate(double focusMorgans,
            double oldEndpointMorgans, double newEndpointMorgans) {
        double oldDist = Math.abs(oldEndpointMorgans - focusMorgans);
        double newDist = Math.abs(newEndpointMorgans - focusMorgans);
        return oldDist==0.0 ? false
                : Math.abs((newDist - oldDist)/oldDist) > maxRelChange;
    }
}
