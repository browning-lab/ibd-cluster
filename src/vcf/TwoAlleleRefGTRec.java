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
package vcf;

import blbutil.Const;
import ints.IndexArray;
import ints.IntArray;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code TwoAlleleRefGTRec} represent represents phased,
 * non-missing genotypes for a list of reference samples at a single diallelic
 * marker.</p>
 *
 * <p>Class {@code TwoAlleleRefGTRec} stores haplotypes that carry
 * the minor allele.</p>
 *
 * <p>Instances of class {@code TwoAlleleRefGTRec} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class TwoAlleleRefGTRec implements RefGTRec {

    private final Marker marker;
    private final Samples samples;
    private final int nHaps;
    private final int majorAllele;
    private final int minorAllele;
    private final int[] minorAlleles;

    @Override
    public long estBytes() {
        int overhead = 12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + 3*8;  // assume 8 bytes per reference
        estBytes += 12; // 4 bytes per int
        estBytes += 4*(1 + minorAlleles.length); // 4 bytes to store array length
        return estBytes;
    }

    /**
     * Constructs a new {@code TwoAlleleRefGTRec} instance from the
     * specified data.
     *
     * @param rec the phased, non-missing genotype data
     * @throws IllegalArgumentException if {@code rec.marker().nAlleles() != 2}
     * @throws NullPointerException if {@code rec == null}
     */
    public TwoAlleleRefGTRec(RefGTRec rec) {
        if (rec.marker().nAlleles()!=2) {
            throw new IllegalArgumentException(
                    String.valueOf(rec.marker().nAlleles()!=2));
        }
        int[][] hapIndices = rec.alleleToHaps();
        int majAllele = 0;
        while (hapIndices[majAllele]!=null) {
            ++majAllele;
        }
        this.marker = rec.marker();
        this.samples = rec.samples();
        this.nHaps = rec.size();
        this.majorAllele = majAllele;
        this.minorAllele = 1 - majAllele;
        this.minorAlleles = hapIndices[minorAllele];
    }

    /**
     * Constructs a new {@code TwoAlleleRefGTRec} instance from the
     * specified data.
     *
     * @param gtp a VCF record parser that extracts sample genotypes
     * @throws IllegalArgumentException if the VCF record contains an
     * unphased genotype or missing allele
     * @throws IllegalArgumentException if {@code gtp.nAlleles() != 2}
     * @throws IllegalArgumentException if a format error is detected in the
     * VCF record
     * @throws NullPointerException if {@code gtp == null}
     */
    public TwoAlleleRefGTRec(VcfRecGTParser gtp) {
        if (gtp.nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(gtp.nAlleles()));
        }
        int[][] nonMajIndices = gtp.nonMajRefIndices();
        this.marker = gtp.marker();
        this.samples = gtp.samples();
        this.nHaps = 2*gtp.nSamples();
        this.majorAllele = nonMajIndices[0]==null ? 0 : 1;
        this.minorAllele = 1 - majorAllele;
        this.minorAlleles = nonMajIndices[minorAllele];
    }

    /**
     * Constructs a new {@code TwoAlleleRefGTRec} instance from the
     * specified data. The specified {@code hapIndices} array is required to
     * have length 2 and contain exactly one {@code null} element.
     * The {@code null} element should be the major allele because this is
     * most memory-efficient, but this requirement is not
     * enforced.

     * @param marker the marker
     * @param samples the samples
     * @param hapIndices whose {@code j}-th element is a list of haplotypes
     * sorted in increasing order that carry the {@code j}-th allele, or is
     * {@code null}
     *
     * @throws IllegalArgumentException if {@code marker.nAlleles() != 2}
     * @throws IllegalArgumentException if
     * {@code marker.nAlleles() != hapIndices.length}
     * @throws IllegalArgumentException if the {@code (hapIndices} does
     * not contain exactly one {@code null} element
     * @throws IllegalArgumentException if the non-null element of
     * {@code hapIndices} is not a sorted list of distinct haplotype indices
     * between 0 (inclusive) and {@code 2*samples.size()} (exclusive)
     * @throws NullPointerException if
     * {@code marker == null || samples == null || hapIndices == null}
     */
    public TwoAlleleRefGTRec(Marker marker, Samples samples,
            int[][] hapIndices) {
        if (marker.nAlleles()!=2) {
            throw new IllegalArgumentException(String.valueOf(marker.nAlleles()));
        }
        this.marker = marker;
        this.samples = samples;
        this.nHaps = 2*samples.size();
        this.majorAllele = AlleleRefGTRec.checkIndicesAndReturnNullIndex(hapIndices, nHaps);
        this.minorAllele = 1 - majorAllele;
        this.minorAlleles = hapIndices[minorAllele].clone();
    }

    @Override
    public int[][] alleleToHaps() {
        int[][] hapIndices = new int[2][];
        hapIndices[minorAllele] = minorAlleles.clone();
        return hapIndices;
    }

    @Override
    public IndexArray hapToAllele() {
        return new IndexArray(toIntArray(), 2);
    }

    private IntArray toIntArray() {
        int[] ia = IntStream.range(0, nHaps)
                .map(h -> majorAllele)
                .toArray();
        for (int h : minorAlleles) {
            ia[h] = minorAllele;
        }
        return IntArray.packedCreate(ia, 2);
    }

    @Override
    public int nAlleleCodedHaps() {
        return minorAlleles.length;
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
        return nHaps;
    }

    @Override
    public Marker marker() {
        return marker;
    }

    @Override
    public int get(int hap) {
        if (hap < 0 || hap >= nHaps) {
            throw new IndexOutOfBoundsException(String.valueOf(hap));
        }
        if (Arrays.binarySearch(minorAlleles, hap) >= 0) {
            return minorAllele;
        }
        else {
            return majorAllele;
        }
    }

    @Override
    public boolean isAlleleCoded() {
        return true;
    }

    @Override
    public int majorAllele() {
        return majorAllele;
    }

    @Override
    public int[] alleleCounts() {
        int[] alCnts = new int[2];
        alCnts[majorAllele] = nHaps - minorAlleles.length;
        alCnts[minorAllele] = minorAlleles.length;
        return alCnts;
    }

    @Override
    public int alleleCount(int allele) {
        if (allele==majorAllele) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return minorAlleles.length;
        }
    }

    @Override
    public int hapIndex(int allele, int copy) {
        if (allele==majorAllele) {
            throw new IllegalArgumentException("major allele");
        }
        else {
            return minorAlleles[copy];
        }
    }

    @Override
    public boolean isCarrier(int allele, int hap) {
        return get(hap)==allele;
    }

    /**
     * Returns the data represented by {@code this} as a string VCF record
     * with correct INFO/AN and INFO/AC fields and with FORMAT/GT as the
     * only FORMAT field.
     * @return the data represented by {@code this} as a string VCF record
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        marker.appendFirst8Fields(sb, nHaps, alleleCounts());
        sb.append(Const.tab);
        sb.append("GT");
        int index = 0;
        int nextMinorHap = (index<minorAlleles.length) ? minorAlleles[index++] : nHaps;
        for (int h=0; h<nHaps; ++h) {
            sb.append((h & 0b1)==0 ? Const.tab : Const.phasedSep);
            if (h==nextMinorHap) {
                sb.append(minorAllele);
                nextMinorHap = (index<minorAlleles.length) ? minorAlleles[index++] : nHaps;
            }
            else {
                sb.append(majorAllele);
            }
        }
        return sb.toString();
    }

    @Override
    public int nMaps() {
        return 1;
    }

    @Override
    public IntArray[] maps() {
        return new IntArray[] {toIntArray()};
    }

    @Override
    public IntArray map(int index) {
        if (index!=0) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return toIntArray();
    }
}
