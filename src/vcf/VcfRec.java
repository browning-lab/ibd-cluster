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
import blbutil.StringUtil;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.Map;

/**
 * <p>Class {@code VcfRec} represents a VCF record.  If one allele in a
 * diploid genotype is missing, then both alleles are set to missing.
 * </p>
 * <p>Instances of class {@code VcfRec} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfRec implements GTRec {

    private static final int SAMPLE_OFFSET = 9;

    private final VcfHeader vcfHeader;
    private final String vcfRecord;
    private final int[] delimiters;

    private final String[] formatFields;
    private final Map<String, Integer> formatMap;

    private final GTRec gtRec;

    @Override
    public long estBytes() {
        int overhead = (6 + formatFields.length + 2*formatMap.size())*12; // assume 12 bytes overhead per "owned" object
        long estBytes = overhead + (8 + formatFields.length + 2*formatMap.size())*8;  // assume 8 bytes per reference
        estBytes += (2 + vcfRecord.length())*2; // 4 bytees per string to store length
        estBytes += (1 + delimiters.length)*4; // 4 bytes per array to store length;
        for (String s : formatFields) {
            estBytes += (2 + s.length())*2;
        }
        for (String s : formatMap.keySet()) {
            estBytes += (2 + s.length() + 2)*4; // 4 bytes per string to store length
        }
        if (gtRec!=null) {
            estBytes += gtRec.estBytes(); // 4 bytes per array to store length
        }
        return estBytes;
    }

    /**
     * Returns the VCF genotype index for the specified pair of alleles.
     * @param a1 the first allele
     * @param a2 the second allele
     * @return the VCF genotype index for the specified pair of alleles
     * @throws IllegalArgumentException if {@code a1 < 0 || a2 < 0}
     */
    public static int gtIndex(int a1, int a2) {
        if (a1 < 0) {
            throw new IllegalArgumentException("a1<0: " + a1);
        }
        if (a2 < 0) {
            throw new IllegalArgumentException("a2<0: " + a2);
        } else if (a1 < a2) {
            return (a2 * (a2 + 1)) / 2 + a1;
        } else {
            return (a1 * (a1 + 1)) / 2 + a2;
        }
    }

    /**
     * Constructs a new {@code VcfRec} instance from a VCF record.
     *
     * @param vcfHeader meta-information lines and header line for the
     * specified VCF record.
     * @param vcfRecord a VCF record with a GL format field corresponding to
     * the specified {@code vcfHeader} object
     * @param stripOptionalFields {@code true} if the VCF record's ID,
     * QUAL, FILTER, and INFO subfields should be discarded
     *
     * @throws IllegalArgumentException if the VCF record does not have a
     * GT format field
     * @throws IllegalArgumentException if a VCF record format error is
     * detected
     * @throws IllegalArgumentException if there are not
     * {@code vcfHeader.nHeaderFields()} tab-delimited fields in the
     * specified VCF record
     * @throws NullPointerException if
     * {@code (vcfHeader == null) || (vcfRecord == null)}
     */
    public VcfRec(VcfHeader vcfHeader, String vcfRecord,
            boolean stripOptionalFields) {
        this.vcfHeader = vcfHeader;
        this.vcfRecord = vcfRecord;
        this.delimiters = delimiters(vcfHeader, vcfRecord);
        this.formatFields = formats(format());
        this.formatMap = formatToIndexMap(vcfHeader, vcfRecord, formatFields);
        if (formatMap.containsKey("GT")==false) {
            String s = "Missing FORMAT/GT field: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        this.gtRec = new BasicGTRec(
                new VcfRecGTParser(vcfHeader, vcfRecord, stripOptionalFields));
    }

    private static int[] delimiters(VcfHeader vcfHeader, String vcfRecord) {
        int nFields = vcfHeader.nHeaderFields();
        int[] delimiters = new int[nFields + 1];
        delimiters[0] = -1;
        for (int j=1; j<nFields; ++j) {
            delimiters[j] = vcfRecord.indexOf(Const.tab, delimiters[j-1] + 1);
            if (delimiters[j] == -1) {
                fieldCountError(vcfHeader, vcfRecord);
            }
        }
        if (vcfRecord.indexOf(Const.tab, delimiters[nFields-1] + 1) != -1) {
            fieldCountError(vcfHeader, vcfRecord);
        }
        delimiters[nFields] = vcfRecord.length();
        return delimiters;
    }

    private static void fieldCountError(VcfHeader vcfHeader, String vcfRecord) {
        String src = "File source: " + vcfHeader.src();
        String[] fields = StringUtil.getFields(vcfRecord, Const.tab);
        String s = "VCF header line has " + vcfHeader.nHeaderFields()
                + " fields, but data line has " + fields.length + " fields"
                + Const.nl + "File source:" + src
                + Const.nl + Arrays.toString(fields);
        throw new IllegalArgumentException(s);
    }

    private String[] formats(String formats) {
        if (formats.equals(Const.MISSING_DATA_STRING) || formats.isEmpty()) {
            String s = "missing format field: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        String[] fields =  StringUtil.getFields(formats, Const.colon);
        for (String f : fields) {
            if (f.isEmpty()) {
                String s = "missing format in format subfield list: " + vcfRecord;
                throw new IllegalArgumentException(s);
            }
        }
        return fields;
    }

    private static Map<String, Integer> formatToIndexMap(VcfHeader vcfHeader,
            String vcfRecord, String[] formatFields) {
        if (vcfHeader.nSamples()==0) {
            return Collections.emptyMap();
        }
        Map<String, Integer> map = new HashMap<>(formatFields.length);
        for (int j=0; j<formatFields.length; ++j) {
            map.put(formatFields[j], j);
        }
        if (map.containsKey("GT") && map.get("GT")!=0) {
            String s = "GT format is not first format: " + vcfRecord;
            throw new IllegalArgumentException(s);
        }
        return map;
    }

    /* returns exclusive end */
    private int formatSubfieldEnd(int start) {
        while (start < vcfRecord.length()) {
            char c = vcfRecord.charAt(start);
            if (c == Const.colon || c == Const.tab) {
                return start;
            }
            ++start;
        }
        return start;
    }

    /**
     * Returns the QUAL field.
     * @return the QUAL field
     */
    public String qual() {
        return vcfRecord.substring(delimiters[5] + 1, delimiters[6]);
    }

    /**
     * Returns the FILTER field.
     * @return the FILTER field
     */
    public String filter() {
        return vcfRecord.substring(delimiters[6] + 1, delimiters[7]);
    }

    /**
     * Returns the INFO field.
     * @return the INFO field
     */
    public String info() {
        return vcfRecord.substring(delimiters[7] + 1, delimiters[8]);
    }

    /**
     * Returns the FORMAT field.  Returns the empty string ("") if the FORMAT
     * field is missing.
     * @return the FORMAT field
     */
    public String format() {
        if (delimiters.length > 9) {
            return vcfRecord.substring(delimiters[8] + 1, delimiters[9]);
        }
        else {
            return "";
        }
    }

    /**
     * Returns the number of FORMAT subfields.
     * @return the number of FORMAT subfields
     */
    public int nFormatSubfields() {
        return formatFields.length;
    }

    /**
     * Returns the specified FORMAT subfield.
     * @param subfieldIndex a FORMAT subfield index
     * @return the specified FORMAT subfield
     *
     * @throws IndexOutOfBoundsException if
     * {@code subfieldIndex < 0 || subfieldIndex >= this.nFormatSubfields()}
     */
    public String formatSubfield(int subfieldIndex) {
        if (formatFields==null) {
            throw new IllegalArgumentException("No format exists");
        }
        return formatFields[subfieldIndex];
    }

    /**
     * Returns {@code true} if the specified FORMAT subfield is
     * present, and returns {@code false} otherwise.
     * @param formatCode a FORMAT subfield code
     * @return {@code true} if the specified FORMAT subfield is
     * present
     */
    public boolean hasFormat(String formatCode) {
        return formatMap.get(formatCode)!=null;
    }

    /**
     * Returns the index of the specified FORMAT subfield if the
     * specified subfield is defined for this VCF record, and returns -1
     * otherwise.
     * @param formatCode the format subfield code
     * @return the index of the specified FORMAT subfield if the
     * specified subfield is defined for this VCF record, and {@code -1}
     * otherwise
     */
    public int formatIndex(String formatCode) {
        Integer index = formatMap.get(formatCode);
        return (index==null) ? -1 : index;
    }

    /**
     * Returns the data for the specified sample.
     * @param sample a sample index
     * @return the data for the specified sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.size()}
     */
    public String sampleData(int sample) {
        int index = vcfHeader.unfilteredSampleIndex(sample);
        return vcfRecord.substring(delimiters[index + SAMPLE_OFFSET] + 1,
                delimiters[index + SAMPLE_OFFSET + 1]);
    }

    /**
     * Returns the specified data for the specified sample.
     * @param sample a sample index
     * @param formatCode a FORMAT subfield code
     * @return the specified data for the specified sample
     *
     * @throws IllegalArgumentException if
     * {@code this.hasFormat(formatCode)==false}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.size()}
     */
    public String sampleData(int sample, String formatCode) {
        Integer formatIndex = formatMap.get(formatCode);
        if (formatIndex==null) {
            String s = "missing format data: " + formatCode;
            throw new IllegalArgumentException(s);
        }
        return VcfRec.this.sampleData(sample, formatIndex);
    }

    /**
     * Returns the specified data for the specified sample.
     * @param sample a sample index
     * @param subfieldIndex a FORMAT subfield index
     * @return the specified data for the specified sample
     *
     * @throws IndexOutOfBoundsException if
     * {@code field < 0 || field >= this.nFormatSubfields()}
     * @throws IndexOutOfBoundsException if
     * {@code sample < 0 || sample >= this.size()}
     */
    public String sampleData(int sample, int subfieldIndex) {
        if (subfieldIndex < 0 || subfieldIndex >= formatFields.length) {
            throw new IndexOutOfBoundsException(String.valueOf(subfieldIndex));
        }
        int index = SAMPLE_OFFSET + vcfHeader.unfilteredSampleIndex(sample);
        int start = delimiters[index] + 1;
        for (int j = 0; j < subfieldIndex; ++j) {
            int end = formatSubfieldEnd(start);
            if (end==vcfRecord.length() || vcfRecord.charAt(end)==Const.tab) {
                return ".";
            }
            else {
                start = end + 1;
            }
        }
        int end = formatSubfieldEnd(start);
        if (end==start) {
            return ".";
        }
        else {
            return vcfRecord.substring(start, end);
        }
    }

    /**
     * Returns an array of length {@code this.size()}
     * containing the specified FORMAT subfield data for each sample.  The
     * {@code k}-th element of the array is the specified FORMAT subfield data
     * for the {@code k}-th sample.
     * @param formatCode a format subfield code
     * @return an array of length {@code this.size()}
     * containing the specified FORMAT subfield data for each sample
     *
     * @throws IllegalArgumentException if
     * {@code this.hasFormat(formatCode) == false}
     */
    public String[] formatData(String formatCode) {
        Integer formatIndex = formatMap.get(formatCode);
        if (formatIndex==null) {
            String s = "missing format data: " + formatCode;
            throw new IllegalArgumentException(s);
        }
        String[] sa = new String[vcfHeader.nSamples()];
        for (int j=0; j<sa.length; ++j) {
            sa[j] = sampleData(j, formatIndex);
        }
        return sa;
    }

    @Override
    public Samples samples() {
        return vcfHeader.samples();
    }



    /**
     * Returns the VCF meta-information lines and the VCF header line.
     * @return the VCF meta-information lines and the VCF header line
     */
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    @Override
    public Marker marker() {
        return gtRec.marker();
    }

    @Override
    public int get(int hap) {
        return gtRec == null ? -1 : gtRec.get(hap);
    }

    @Override
    public boolean isPhased(int sample) {
        return gtRec == null ? false : gtRec.isPhased(sample);
    }

    @Override
    public boolean isPhased() {
        return  gtRec == null ? false : gtRec.isPhased();
    }

    @Override
    public int size() {
        return 2*vcfHeader.nSamples();
    }

    /**
     * Returns the VCF record.
     * @return the VCF record
     */
    @Override
    public String toString() {
        return vcfRecord;
    }
}
