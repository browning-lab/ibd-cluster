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
package ibdcluster;

import ints.IndexArray;
import ints.IntArray;

/**
 * <p>Class {@code AggregateMarker} represents an aggregate marker
 * whose alleles are haplotype sequences.</p>
 *
 * <p>Instances of class {@code AggregateMarker} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class AggMarker {

    private final Position position;
    private final IndexArray alleles;

    /**
     * Constructs a new {@code AggMarker} instance from the specified data.
     * @param position the genomic position
     * @param alleles the map from haplotype to aggregate marker allele index
     * @throws NullPointerException if
     * {@code (position == null) || (alleles == null)}
     */
    public AggMarker(Position position, IndexArray alleles) {
        if (alleles == null) {
            throw new NullPointerException(IntArray.class.toString());
        }
        if (position==null) {
            throw new NullPointerException(Position.class.toString());
        }
        this.position = position;
        this.alleles = alleles;
    }

    /**
     * Returns the genomic position.
     * @return the genomic position
     */
    public Position position() {
        return position;
    }

    /**
     * Returns the map from haplotype to aggregate marker allele index.
     * @return the map from haplotype to aggregate marker allele index
     */
    public IndexArray alleles() {
        return alleles;
    }
}
