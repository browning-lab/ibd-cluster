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
import vcf.GT;

/**
 * <p>Class {@code IbdEndsUtils} contains static utility methods
 * used for estimating IBD segment endpoints.</p>
 *
 * <p>The static methods in class {@code IbdEndsUtils} are thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ClustUtils {

    private static final int BASE_POS_BACKOFF = 5_000_000;
    private static final float MORGAN_POS_BACKOFF = 0.05f;

    private ClustUtils() {
        // private constructor to prevent instantiation
    }

    /**
     * Returns the estimated base position of the specified Morgan
     * position.  The base position is estimated by linear interpolation.
     * The contract for this method is undefined if the {@code basePos}
     * or {@code morganPos} arrays are not sorted in increasing order.
     * @param basePos a list giving base positions corresponding to
     * the specified Morgan positions
     * @param morganPos a list giving Morgan positions corresponding to
     * the specified base positions
     * @param inputMorganPos a Morgan position
     * @return the estimated base position of the specified Morgan
     * position.
     * @throws IllegalArgumentException if {@code basePos.size() < 2}
     * @throws IllegalArgumentException if
     * {@code basePos.size() != morganPos.size()}
     * @throws NullPointerException if
     * {@code basePos == null || morganPos == null}
     */
    public static int basePos(WrappedIntArray basePos,
            DoubleArray morganPos, double inputMorganPos) {
        if (basePos.size()<2) {
            throw new IllegalArgumentException("insufficient data");
        }
        if (basePos.size()!=morganPos.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        int index = morganPos.binarySearch(inputMorganPos);
        if (index>=0) {
            return basePos.get(index);
        }
        else {
            int mapSizeM1 = morganPos.size() - 1;
            int insPt = -index - 1;
            int aIndex = insPt - 1;
            int bIndex = insPt;
            if (aIndex==mapSizeM1) {
                insPt = morganPos.binarySearch(morganPos.get(mapSizeM1) - MORGAN_POS_BACKOFF);
                if (insPt<0) {
                    insPt = -insPt - 2;
                }
                assert insPt<mapSizeM1;
                aIndex = Math.max(insPt, 0);
                bIndex = mapSizeM1;
            }
            else if (bIndex==0) {
                insPt = morganPos.binarySearch(morganPos.get(0) + MORGAN_POS_BACKOFF);
                if (insPt<0) {
                    insPt = -insPt - 1;
                }
                assert insPt>0;
                aIndex = 0;
                bIndex = Math.min(insPt, mapSizeM1);
            }
            double x = inputMorganPos;
            double a = morganPos.get(aIndex);
            double b = morganPos.get(bIndex);
            double fa = basePos.get(aIndex);
            double fb = basePos.get(bIndex);
            return (int) Math.round(fa + ( ((x-a)/(b-a)) * (fb-fa) ));
        }
    }

    /**
     * Returns the estimated Morgan position of the specified base
     * position.  The base position is estimated by linear interpolation.
     * The contract for this method is undefined if the {@code basePos}
     * or {@code morganPos} arrays are not sorted in increasing order.
     * @param basePos a list giving base positions corresponding to
     * the specified Morgan positions
     * @param morganPos a list giving Morgan positions corresponding to
     * the specified base positions
     * @param inputBasePos a base position
     * @return the estimated Morgan position of the specified base
     * position.
     * @throws IllegalArgumentException if {@code basePos.size() < 2}
     * @throws IllegalArgumentException if
     * {@code basePos.size() != morganPos.size()}
     * @throws NullPointerException if
     * {@code basePos == null || morganPos == null}
     */
    public static double morganPos(WrappedIntArray basePos,
            DoubleArray morganPos, int inputBasePos) {
        if (basePos.size()<2) {
            throw new IllegalArgumentException("insufficient data");
        }
        if (basePos.size()!=morganPos.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        int index = basePos.binarySearch(inputBasePos);
        if (index>=0) {
            return morganPos.get(index);
        }
        else {
            int mapSizeM1 = basePos.size() - 1;
            int insPt = -index - 1;
            int aIndex = insPt - 1;
            int bIndex = insPt;
            if (aIndex==mapSizeM1) {
                insPt = basePos.binarySearch(basePos.get(mapSizeM1) - BASE_POS_BACKOFF);
                if (insPt<0) {
                    insPt = -insPt - 2;
                }
                assert insPt<mapSizeM1;
                aIndex = Math.max(insPt, 0);
                bIndex = mapSizeM1;
            }
            else if (bIndex==0) {
                insPt = basePos.binarySearch(basePos.get(0) + BASE_POS_BACKOFF);
                if (insPt<0) {
                    insPt = -insPt - 1;
                }
                assert insPt>0;
                aIndex = 0;
                bIndex = Math.min(insPt, mapSizeM1);
            }
            int x = inputBasePos;
            int a = basePos.get(aIndex);
            int b = basePos.get(bIndex);
            double fa = morganPos.get(aIndex);
            double fb = morganPos.get(bIndex);
            return fa + ( ((double) (x-a)/(b-a)) * (fb-fa) );
        }
    }

    /**
     * Returns the probability that an IBD segment has its right endpoint
     * less than {@code y} Morgans from the left endpoint.
     * @param y a distance in Morgans
     * @param ne the constant effective population size
     * @return the probability that an IBD segment has its right endpoint
     * less than {@code y} Morgans from the left endpoint
     * @throws IllegalArgumentException if {@code y <= 0 || Double.isNaN(y)}
     * @throws IllegalArgumentException if
     * {@code ne <= 0.0 || Double.isFinite(ne) == false}
     */
    public static double F(double y, double ne) {
        if (y<=0 || Double.isNaN(y)) {
            throw new IllegalArgumentException(String.valueOf(y));
        }
        if (ne<=0.0 || Double.isFinite(ne)==false) {
            throw new IllegalArgumentException(String.valueOf(ne));
        }
        double den = 2*ne*Math.expm1(2*y) + 1d;
        return 1.0 - 1.0/den;
    }

    /**
     * Returns a value {@code y} such that {@code IbdEndsUtils.F(y, ne)} is
     * approximately equal to {@code p}.
     * @param p a probability satisfying {@code 0 < p && p < 1}
     * @param ne the effective population size
     * @return a value {@code y} such that {@code F(y, ne)} is approximately
     * equal to {@code p}
     * @throws IllegalArgumentException if
     * {@code p <= 0.0 || p >= 1.0 || Double.isNaN(p)}
     * @throws IllegalArgumentException if
     * {@code ne <= 0.0 || Double.isFinite(ne) == false}
     */
    public static double invF(double p, double ne) {
        if (p<=0d || p>=1d || Double.isNaN(p)) {
            throw new IllegalArgumentException(String.valueOf(p));
        }
        if (ne<=0.0 || Double.isFinite(ne)==false) {
            throw new IllegalArgumentException(String.valueOf(ne));
        }
        double d = 2*ne*(1 - p);
        return 0.5*Math.log((p + d)/d);
    }
}
