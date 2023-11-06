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

/**
 * <p>Class {@code Position} stores a genomic position.</p>
 *
 * <p>Instances of class {@code Position} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Position {

    private final int chrom;
    private final int pos;
    private final double genPos;

    /**
     * Constructs a new {@code Position} from the specified data.
     *
     * @param chromIndex a chromosome index
     * @param pos a position
     * @param genPos a genetic position
     */
    public Position(int chromIndex, int pos, double genPos) {
        this.chrom = chromIndex;
        this.pos = pos;
        this.genPos = genPos;
    }

    /**
     * Return the chromosome index.
     * @return the chromosome index
     */
    public int chrom() {
        return chrom;
    }

    /**
     * Return the position.
     * @return the position
     */
    public int pos() {
        return pos;
    }

    /**
     * Return the genetic position
     * @return the genetic position
     */
    public double genPos() {
        return genPos;
    }
}
