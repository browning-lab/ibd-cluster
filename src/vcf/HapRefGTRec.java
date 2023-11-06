/*
 * Copyright 2023 Brian L. Browning
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
package vcf;

import bref4.BrefRec;
import ints.IndexArray;
import ints.IntArray;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code HapRefGTRec}  represents phased, non-missing
 * genotypes for a list of reference samples at a single marker.
 * Genotype emission probabilities are determined by the sample
 * genotypes.
 * </p>
 * <p>Instances of class {@code HapRefGTRec} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class HapRefGTRec implements RefGTRec {

    private final Marker marker;
    private final Samples samples;
    private final IntArray hapToSeq;
    private final IntArray seqToAllele;

    @Override
    public long estBytes() {
        // Note: The same hapToSeq array can be counted in multiple records
        int overhead = 2*12;                    // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 4*8;         // assume 8 bytes per reference
        estBytes += (1 + hapToSeq.size())*4;    // 4 bytes per array to store length;
        estBytes += (1 + seqToAllele.size())*4; // 4 bytes per array to store length;
        return estBytes;
    }

    /**
     * Creates a new {@code HapRefGTRec} instance with phased,
     * non-missing genotypes from the specified marker, samples,
     * and haplotype alleles.  The contract for the constructed object
     * is undefined if any element of {@code hapToSeq} is negative or
     * greater than or equal to {@code hapToAllele.size()} or if any element
     * of {@code hapToAllele} is negative or greater than or equal to
     * {@code marker.nAlleles()}.
     *
     * @param marker the marker
     * @param samples the samples
     * @param hapToSeq an array whose {@code j}-th element is the index
     * of the distinct allele sequence carried by the {@code j}-th haplotype
     * @param seqToAllele an array whose {@code j}-th element is the marker
     * allele carried by the {@code j}-th distinct allele sequence
     *
     * @throws IllegalArgumentException if
     * {@code hapToSeq.size() != 2*samples.size()}
     * @throws NullPointerException if any parameter is {@code null}
     */
    public HapRefGTRec(Marker marker, Samples samples, IntArray hapToSeq,
        IntArray seqToAllele) {
        if (hapToSeq.size() != 2*samples.size()) {
            throw new IllegalArgumentException("inconsistent data");
        }
        this.marker = marker;
        this.samples = samples;
        this.hapToSeq = hapToSeq;
        this.seqToAllele = seqToAllele;
    }

    @Override
    public boolean isPhased(int sample) {
        if (sample < 0 || sample >= this.samples().size()) {
            throw new IndexOutOfBoundsException(String.valueOf(sample));
        }
        return true;
    }

    /**
     * Returns {@code true}.
     * @return {@code true}
     */
    @Override
    public boolean isPhased() {
        return true;
    }


    @Override
    public Samples samples() {
        return samples;
    }

    @Override
    public int size() {
        return hapToSeq.size();
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public int[][] alleleToHaps() {
        int[] alCnts = alleleCounts();
        int majAllele = 0;
        for (int al=1; al<alCnts.length; ++al) {
            if (alCnts[al] > alCnts[majAllele]) {
                majAllele = al;
            }
        }
        int[][] hapIndices = new int[alCnts.length][];
        for (int al=0; al<alCnts.length; ++al) {
            if (al!=majAllele) {
                hapIndices[al] = new int[alCnts[al]];
            }
        }
        Arrays.fill(alCnts, 0);
        for (int h=0, n=size(); h<n; ++h) {
            int al = get(h);
            if (al!=majAllele) {
                hapIndices[al][alCnts[al]++] = h;
            }
        }
        return hapIndices;
    }

    @Override
    public IndexArray hapToAllele() {
        int[] alleles = IntStream.range(0, size())
                .map(h -> seqToAllele.get(hapToSeq.get(h)))
                .toArray();
        return new IndexArray(alleles, marker.nAlleles());
    }

    @Override
    public int nAlleleCodedHaps() {
        return BrefRec.nonNullCnt(alleleToHaps());
    }

    @Override
    public boolean isAlleleCoded() {
        return false;
    }

    @Override
    public int majorAllele() {
        return majorAllele(alleleCounts());
    }

    private int majorAllele(int[] alCnts) {
        int majAl = 0;
        for (int al=1; al<alCnts.length; ++al) {
            if (alCnts[al]>alCnts[majAl]) {
                majAl = al;
            }
        }
        return majAl;
    }

    @Override
    public int[] alleleCounts() {
        int[] alCnts = new int[marker.nAlleles()];
        for (int h=0, n=size(); h<n; ++h) {
            ++alCnts[get(h)];
        }
        return alCnts;
    }

    @Override
    public int alleleCount(int allele) {
        int[] alCnts = alleleCounts();
        if (allele==majorAllele(alCnts)) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return alCnts[allele];
        }
    }

    @Override
    public int get(int hap) {
        return seqToAllele.get(hapToSeq.get(hap));
    }

    @Override
    public int hapIndex(int allele, int copy) {
        int[][] hapIndices = alleleToHaps();
        if (hapIndices[allele]==null) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return hapIndices[allele][copy];
        }
    }

    @Override
    public boolean isCarrier(int allele, int hap) {
        return get(hap)==allele;
    }

    /**
     * Returns the data represented by {@code this} as a VCF
     * record with a GT format field. The returned VCF record
     * will have missing QUAL and INFO fields, will have "PASS"
     * in the filter field, and will have a GT format field.
     * @return the data represented by {@code this} as a VCF
     * record with a GT format field
     */
    @Override
    public String toString() {
        return GTRec.toVcfRec(this);
    }

    @Override
    public int nMaps() {
        return 2;
    }

    @Override
    public IntArray[] maps() {
        return new IntArray[] {hapToSeq, seqToAllele};
    }

    @Override
    public IntArray map(int index) {
        switch (index) {
            case 0:
                return hapToSeq;
            case 1:
                return seqToAllele;
            default:
                throw new IndexOutOfBoundsException(String.valueOf(index));
        }
    }
}
