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
package blbutil;

import vcf.GTRec;
import vcf.VcfHeader;

/**
 * <p>An iterator for records in a VCF file.  Each record contains
 * data for the same set of samples.
 *</p>
 * Instances of class {@code VcfFileIt} are not thread-safe.
 *
 * @param <E> the type of the elements returned by this iterator's
 * {@code next()} method.
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface VcfFileIt<E extends GTRec> extends SampleFileIt<E> {

    /**
     * Returns the VCF meta-information lines and header line.
     * @return the VCF meta-information lines and header line
     */
    VcfHeader vcfHeader();
}
