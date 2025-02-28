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

import blbutil.FileIt;
import blbutil.VcfFileIt;
import java.io.File;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Deque;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.BiFunction;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * <p>Class {@code VcfIt} represents  an iterator whose {@code next()}
 * method returns an object storing data from a VCF record.
 * </p>
 * <p>Instances of class {@code VcfIt} are not thread-safe.
 * </p>
 * <p>Methods of this class will terminate the Java Virtual Machine with
 * an error message if an I/O error or file format error is detected.
 * </p>
 * @param <E> the type parameter
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class VcfIt<E extends GTRec> implements VcfFileIt<E> {

    private final VcfHeader vcfHeader;
    private final FileIt<String> it;
    private final Function<String, E> mapper;
    private String next;
    private final Predicate<Marker> markerFilter;

    private final int bufferSize;
    private final String[] stringBuffer;
    private final Deque<E> recBuffer;

    /**
     * The default number of VCF records stored in a buffer.
     */
    public static final int DEFAULT_BUFFER_SIZE = 1<<10;

   /**
     * A function mapping a string VCF record with a FORMAT/GT field
     * to a memory-efficient {@code GTRec} object. Phase status is stored
     * per record.  All genotypes are considered to be unphased if any
     * genotype is unphased or if any allele is missing.
     */
    public static final BiFunction<VcfHeader, String, GTRec> TO_LOWMEM_GT_REC
            = (VcfHeader h, String vcfRec) -> {
                return lowMemGTRec(new VcfRecGTParser(h, vcfRec, false));
            };

   /**
     * A function that strips the VCF ID, QUAL, FILTER, and INFO, and FORMAT
     * fields from a VCF record with a FORMAT/GT field and maps the
     * stripped VCF record to a memory-efficient {@code GTRec} object.
     * Phase status is stored per VCF record.  All genotypes are considered
     * to be unphased if any genotype is unphased or if any allele is missing.
     */
    public static final BiFunction<VcfHeader, String, GTRec> TO_STRIPPED_LOWMEM_GT_REC
        = (VcfHeader h, String vcfRec) -> {
            return lowMemGTRec(new VcfRecGTParser(h, vcfRec, true));
        };

        private static GTRec lowMemGTRec(VcfRecGTParser recParser) {
            VcfRecGTParser.HapListRep hlr = recParser.hapListRep();
            int nonMajorAlleleThreshold = (hlr.samples().size()>>7);
            if (hlr.nonmajorAlleleCnt()<=nonMajorAlleleThreshold) {
                if (hlr.marker().nAlleles()==2) {
                    return new LowMafDiallelicGTRec(hlr);
                }
                else {
                    return new LowMafGTRec(hlr);
                }
            }
            else {
                return new BitArrayGTRec(hlr);
            }
        }

   /**
     * A function mapping a string VCF record with a FORMAT/GT field
     * to a {@code BasicGTRec} object. Phase status is stored per-genotype.
     */
    public static final BiFunction<VcfHeader, String, GTRec> TO_BASIC_GT_REC
        = (VcfHeader h, String vcfRec) -> new BasicGTRec(new VcfRecGTParser(h, vcfRec, false));

   /**
     * A function that strips the VCF ID, QUAL, FILTER, and INFO, and FORMAT
     * fields from a VCF record with a FORMAT/GT field and maps the
     * stripped VCF record to a {@code BasicGTRec} object. Phase status is
     * stored per-genotype.
     */
    public static final BiFunction<VcfHeader, String, GTRec> TO_STRIPPED_BASIC_GT_REC
        = (VcfHeader h, String vcfRec) -> new BasicGTRec(new VcfRecGTParser(h, vcfRec, true));

    /**
     * A function mapping a string VCF record with a FORMAT/GT field to
     * a {@code VcfRec} object.
     */
    public static final BiFunction<VcfHeader, String, VcfRec> TO_VCF_REC
            = (VcfHeader h, String vcfRec) -> new VcfRec(h, vcfRec, false);

   /**
     * A function that strips the VCF ID, QUAL, FILTER, and INFO, and FORMAT
     * fields from a VCF record with a FORMAT/GT field and maps the stripped
     * VCF record to a {@code VcfRec} object. Phase status is stored
     * per-genotype.
     */
    public static final BiFunction<VcfHeader, String, VcfRec> TO_STRIPPED_VCF_REC
            = (VcfHeader h, String vcfRec) -> new VcfRec(h, vcfRec, true);

    /**
     * Returns an array containing VCF meta-information lines, the
     * VCF header line, and the first VCF data line.  The returned array
     * will contain all initial lines that begin with the '#' character
     * and the next line.
     * @param src a string describing the source of the VCF file
     * @param it an iterator that returns the lines of a VCF file
     * @return an array containing VCF meta-information lines, the
     * VCF header line, and the first VCF data line
     * @throws NullPointerException if {@code it == null}
     * @throws IllegalArgumentException if all lines returned by the iterator
     * begin with the '#' character
     */
    public static String[] head(String src, FileIt<String> it) {
        String hash = "#";
        List<String> lines = new ArrayList<>(32);
        String line =  it.hasNext() ? it.next() : null;
        while (line!=null && line.startsWith(hash)) {
            lines.add(line);
            line = it.hasNext() ? it.next() : null;
        }
        if (line==null) {
            throw new IllegalArgumentException("ERROR: missing VCF data lines ("
                    + src + ")");
        }
        lines.add(line);
        return lines.toArray(new String[0]);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param recMapper a function mapping string VCF records to
     * {@code GTRec} objects
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if
     * {@code strIt == null || mapFactory == null}
     */
    public static <R extends GTRec> VcfIt<R> create(
            FileIt<String> strIt, BiFunction<VcfHeader, String, R> recMapper) {
        return VcfIt.create(strIt, FilterUtil.acceptAllPredicate(),
                FilterUtil.acceptAllPredicate(), recMapper);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param recMapper a function mapping string VCF records to
     * {@code GTRec} objects
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws NullPointerException if
     * {@code strIt == null || mapFactory == null}
     */
    public static <R extends GTRec> VcfIt<R> create(
            FileIt<String> strIt, Predicate<String> sampleFilter,
            Predicate<Marker> markerFilter,
            BiFunction<VcfHeader, String, R> recMapper) {
        return VcfIt.create(strIt, sampleFilter, markerFilter, recMapper,
                DEFAULT_BUFFER_SIZE);
    }

    /**
     * Create and returns a new {@code VcfIt} instance from the specified
     * objects.
     * @param <R> the type returned by the returned {@code VcfIt}
     * @param strIt an iterator that returns lines of a VCF file
     * @param sampleFilter a sample filter or {@code null}
     * @param markerFilter a marker filter or {@code null}
     * @param recMapper a function mapping string VCF records to
     * {@code GTRec} objects
     * @param bufferSize the requested buffer size
     * @return a new {@code VcfIt} instance
     * @throws IllegalArgumentException if a format error is detected in a
     * line of a VCF file returned by {@code strIt}
     * @throws IllegalArgumentException if {@code bufferSize < 1}
     * @throws NullPointerException if
     * {@code strIt == null || mapFactory == null}
     */
    public static <R extends GTRec> VcfIt<R> create(
            FileIt<String> strIt, Predicate<String> sampleFilter,
            Predicate<Marker> markerFilter,
            BiFunction<VcfHeader, String, R> recMapper, int bufferSize) {
        return new VcfIt<>(strIt, sampleFilter, markerFilter,
                recMapper, bufferSize);
    }

    private VcfIt(FileIt<String> it, Predicate<String> sampleFilter,
            Predicate<Marker> markerFilter,
            BiFunction<VcfHeader, String, E> recMapper,
            int bufferSize) {
        if (bufferSize < 1) {
            throw new IllegalArgumentException(String.valueOf(bufferSize));
        }
        if (markerFilter==null) {
            markerFilter = FilterUtil.acceptAllPredicate();
        }
        String src = it.file()==null ? "stdin" : it.file().getName();
        String[] head = head(src, it);
        String[] nonDataLines = Arrays.copyOf(head, head.length-1);
        String firstDataLine = head[head.length-1];
        boolean[] isDiploid = VcfHeader.isDiploid(firstDataLine);
        this.it = it;
        this.vcfHeader = new VcfHeader(src, nonDataLines, isDiploid, sampleFilter);
        this.mapper = (String s) -> recMapper.apply(vcfHeader, s);
        this.next = firstDataLine;
        this.markerFilter = markerFilter;
        this.bufferSize = bufferSize;
        this.stringBuffer = new String[stringBufferSize(next, bufferSize)];
        this.recBuffer = new ArrayDeque<>(bufferSize);
        fillEmissionBuffer();
    }

    private static int stringBufferSize(String line, int maxBufferSize) {
        long maxMem = Runtime.getRuntime().maxMemory();
        long nBytesPerLine = 2*(line==null ? 0 : line.length());
        long bufferSize = 1 + (maxMem >> 4)/(nBytesPerLine);
        return (int) Math.min(bufferSize, maxBufferSize);
    }

    private void fillEmissionBuffer() {
        assert recBuffer.isEmpty();
        int size = -1;
        while (size!=0 && recBuffer.size()<bufferSize) {
            size = fillStringBuffer(stringBuffer.length);
            if (size>0) {
                List<E> list = IntStream.range(0, size)
                        .parallel()
                        .mapToObj(j -> stringBuffer[j])
                        .map(mapper)
                        .filter(e -> markerFilter.test(e.marker()))
                        .collect(Collectors.toList());
                recBuffer.addAll(list);
            }
        }
    }

    private int fillStringBuffer(int maxSize) {
        int size = 0;
        if (next != null) {
            while (next!=null && size<maxSize) {
                stringBuffer[size++] = next;
                next = readLine(it);
            }
        }
        return size;
    }

    private static String readLine(FileIt<String> it) {
        if (it.hasNext()==false) {
            return null;
        }
        String line = it.next();
        while (line.trim().isEmpty() && it.hasNext()) {
            line = it.next();
        }
        return line;
    }

    @Override
    public void close() {
        next = null;
        it.close();
        recBuffer.clear();
        Arrays.fill(stringBuffer, null);

    }

    /**
     * Returns {@code true} if the iteration has more elements, and returns
     * {@code false} otherwise.
     * @return {@code true} if the iteration has more elements
     */
    @Override
    public boolean hasNext() {
        return !recBuffer.isEmpty();
    }

    /**
     * Returns the next element in the iteration.
     * @return the next element in the iteration
     * @throws NoSuchElementException if the iteration has no more elements.
     */
    @Override
    public E next() {
        if (!hasNext()) {
            throw new NoSuchElementException();
        }
        E first = recBuffer.removeFirst();
        if (recBuffer.isEmpty()) {
            fillEmissionBuffer();
        }
        return first;
    }

    /**
     * The {@code remove} method is not supported by this iterator.
     * @throws UnsupportedOperationException if this method is invoked
     */
    @Override
    public void remove() {
        throw new UnsupportedOperationException(this.getClass().toString());
    }

    @Override
    public File file() {
        return it.file();
    }

    @Override
    public Samples samples() {
        return vcfHeader.samples();
    }

    @Override
    public VcfHeader vcfHeader() {
        return vcfHeader;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(80);
        sb.append(this.getClass().toString());
        sb.append(" : ");
        sb.append(it.file()==null ? "stdin" : it.file().toString());
        return sb.toString();
    }
}
