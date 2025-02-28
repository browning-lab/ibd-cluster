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
import java.util.Comparator;

/**
 * <p>Class {@code HapPairSegment} represents a shared chromosome segment for
 * a pair of haplotypes.</p>
 *
 * <p>Instances of class {@code HapPairSegment} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class HapPairSegment {

    /**
     * {@code ZERO_LENGTH_SEGMENT} represents a {@code HapPairSegment}
     * segment that does not have positive length.
     */
    public static final HapPairSegment ZERO_LENGTH_SEGMENT
            = new HapPairSegment(Integer.MAX_VALUE, Integer.MAX_VALUE,
                    Integer.MAX_VALUE, Integer.MAX_VALUE);

    private final int hap1;
    private final int hap2;
    private final int startPos;
    private final int inclEndPos;

    /**
     * Constructs a new {@code HapPairSegment} instance from the specified data
     * @param hap1 the first haplotype index
     * @param hap2 the second haplotype index
     * @param startPos the first base coordinate in the segment
     * @param inclEndPos the last base coordinate in the segment
     *
     * @throws IllegalArgumentException if
     * {@code  hap1 < 0 || hap2 < 0 || startPos > inclEndPos}
     */
    public HapPairSegment(int hap1, int hap2, int startPos, int inclEndPos) {
        if (hap1<0) {
            throw new IllegalArgumentException(String.valueOf(hap1));
        }
        if (hap2<0) {
            throw new IllegalArgumentException(String.valueOf(hap2));
        }
        else if (startPos>inclEndPos){
            String s = "startPos=" + startPos + " > inclEndPos=" + inclEndPos;
            throw new IllegalArgumentException(s);
        }
        this.hap1 = hap1;
        this.hap2 = hap2;
        this.startPos = startPos;
        this.inclEndPos = inclEndPos;
    }

    /**
     * Compares the specified object with this {@code HapPairSegment} for
     * equality.  Returns {@code  true} if the specified object is a
     * {@code HapPairSegment} instance whose {@code hap1()}, {@code hap2()},
     * {@code startPos()}, and {@code inclEndPos()} methods return the same
     * values as the corresponding methods for {@code  this}, and
     * returns {@code  false} otherwise.
     * @param obj the object to be compared for equality with {@code this}
     *
     * @return {@code  true} if the specified object is a
     * {@code HapPairSegment} instance whose {@code hap1()}, {@code hap2()},
     * {@code startPos()}, and {@code inclEndPos()} methods return the same
     * values as the corresponding methods for {@code  this}
     */
    @Override
    public boolean equals(Object obj) {
        if (this==obj) {
            return true;
        }
        if ((obj instanceof HapPairSegment)==false) {
            return false;
        }
        HapPairSegment other = (HapPairSegment) obj;
        if (this.hap1!=other.hap1) {
            return false;
        }
        if (this.hap2!=other.hap2) {
            return false;
        }
        if (this.startPos!=other.startPos) {
            return false;
        }
        return this.inclEndPos == other.inclEndPos;
    }

     /**
     * Returns a hash value for this object.
     * @return a hash value for this object
     */
    @Override
    public int hashCode() {
        int hash = 3;
        hash = 43 * hash + this.hap1;
        hash = 43 * hash + this.hap2;
        hash = 43 * hash + this.startPos;
        hash = 43 * hash + this.inclEndPos;
        return hash;
    }

    /**
     * Returns a {@code Comparator<HapPairSegment>} whose {@code compare()}
     * method orders {@code HapPairSegment} instances first by
     * {@code this.hap1()}, then by {@code this.hap2()},
     * then by {@code this.startPos()}, and finally by {@code this.inclEndPos()}.
     * @return a {@code Comparator} that orders {@code HapPairSegment}
     * instances by {@code this.hap1()}, then by {@code this.hap2()},
     * then by {@code this.startPos()}, and finally by {@code this.inclEndPos()}
     */
    public static Comparator<HapPairSegment> hapPairComparator() {
        return new Comparator<HapPairSegment>() {

            /**
             * Returns a negative integer, zero, or a positive integer
             * depending on whether {@code seg1} is less than, equal to, or
             * greater than {@code seg2}. The {@code HapPairSegment} instances
             * are ordered first by {@code this.hap1()}, then by
             * {@code this.hap2()}, then by {@code this.startPos()}, and
             * finally by {@code this.inclEndPos()}.
             * @param seg1 the first segment to be compared
             * @param seg2 the second segment to be compared
             * @return a negative integer, zero, or a positive integer depending
             * on whether {@code seg1} is less than, equal to, or greater than
             * {@code seg2}
             * @throws NullPointerException if {@code ((seg1 == null) || (seg2 == null))}
             */
            @Override
            public int compare(HapPairSegment seg1, HapPairSegment seg2) {
                if (seg1.hap1() != seg2.hap1()) {
                    return (seg1.hap1 < seg2.hap1) ? -1 : 1;
                }
                if (seg1.hap2() != seg2.hap2()) {
                    return (seg1.hap2 < seg2.hap2) ? -1 : 1;
                }
                if (seg1.startPos != seg2.startPos) {
                    return (seg1.startPos < seg2.startPos) ? -1 : 1;
                }
                if (seg1.inclEndPos != seg2.inclEndPos) {
                    return (seg1.inclEndPos < seg2.inclEndPos) ? -1 : 1;
                }
                return 0;
            }
        };
    }

    /**
     * Returns a {@code Comparator<HapPairSegment>} whose {@code compare()}
     * method orders {@code HapPairSegment} instances first by
     * {@code this.startPos()}, then by {@code this.inclEndPos()},
     * then by {@code this.hap1()}, and finally by {@code this.hap2()}.
     * @return a {@code Comparator} that orders {@code HapPairSegment}
     * instances by {@code this.startPos()}, then by {@code this.inclEndPos()},
     * then {@code this.hap1()}, and finally by {@code this.hap2()}
     */
    public static Comparator<HapPairSegment> intervalComparator() {
        return new Comparator<HapPairSegment>() {

            /**
             * Returns a negative integer, zero, or a positive integer
             * depending on whether {@code seg1} is less than, equal to, or
             * greater than {@code seg2}. The {@code HapPairSegment} instances
             * are ordered first by {@code this.startPos()}, then by
             * {@code this.inclEnd()}, then by {@code this.hap1()}, and
             * finally by {@code this.hap2()}.
             * @param seg1 the first segment to be compared
             * @param seg2 the second segment to be compared
             * @return a negative integer, zero, or a positive integer depending
             * on whether {@code seg1} is less than, equal to, or greater than
             * {@code seg2}
             * @throws NullPointerException if {@code ((seg1 == null) || (seg2 == null))}
             */
            @Override
            public int compare(HapPairSegment seg1, HapPairSegment seg2) {
                if (seg1.startPos != seg2.startPos) {
                    return (seg1.startPos < seg2.startPos) ? -1 : 1;
                }
                if (seg1.inclEndPos != seg2.inclEndPos) {
                    return (seg1.inclEndPos < seg2.inclEndPos) ? -1 : 1;
                }
                if (seg1.hap1() != seg2.hap1()) {
                    return (seg1.hap1 < seg2.hap1) ? -1 : 1;
                }
                if (seg1.hap2() != seg2.hap2()) {
                    return (seg1.hap2 < seg2.hap2) ? -1 : 1;
                }
                return 0;
            }
        };
    }

    /**
     * Returns the first haplotype index.
     * @return the first haplotype index
     */
    public int hap1() {
        return hap1;
    }

    /**
     * Returns the second haplotype index.
     * @return the second haplotype index
     */
    public int hap2() {
        return hap2;
    }

    /**
     * Returns the first base position in the segment.
     * @return the first base position in the segment
     */
    public int startPos() {
        return startPos;
    }

    /**
     * Returns the last base position in the segment.
     * @return the last base position in the segment
     */
    public int inclEndPos() {
        return inclEndPos;
    }

    /**
     * Returns a string representation of {@code this}.  The exact
     * details of the representation are unspecified and subject to change.
     * @return a string representation of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(hap1);
        sb.append(Const.tab);
        sb.append(hap2);
        sb.append(Const.tab);
        sb.append(startPos);
        sb.append(Const.tab);
        sb.append(inclEndPos);
        return sb.toString();
    }
}
