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

import beagleutil.ChromIds;
import blbutil.Const;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code Partition} implements a disjoint union data
 * structure that stores the partition of IBD clusters.</p>
 *
 * <p>Instances of class {@code Partition} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Partition {

    private final Position position;
    private final int[] parent;
    private final int[] rank;
    private int nSets;

    /**
     * Constructs a new {@code Partition} instance for the specified data.
     * The haplotypes will be initially partitioned into {@code nHaps}
     * singleton sets.
     *
     * @param position a genomic position
     * @param nHaps the number of elements that are partitioned
     * @throws NegativeArraySizeException if {@code nHaps < 0}
     * @throws NullPointerException if {@code position == null}
     */
    public Partition(Position position, int nHaps) {
        if (position==null) {
            throw new NullPointerException(Position.class.toString());
        }
        this.position = position;
        this.parent = IntStream.range(0, nHaps).toArray();
        this.rank = new int[nHaps];
        this.nSets = nHaps;
    }

    /**
     * Returns the number of haplotypes.
     * @return the number of haplotypes
     */
    public int nHaps() {
        return parent.length;
    }

    /**
     * Returns the genomic position.
     * @return the genomic position
     */
    public Position position() {
        return position;
    }

    /**
     * Returns the representative member of the set containing the specified
     * haplotype.  This method performs path compression.
     * @param hap a haplotype index
     * @return the representative member of the set to which the specified
     * haplotype belongs
     * @throws IndexOutOfBoundsException if {@code x < 0 || x >= this.nHaps()}
     */
    public int find(int hap) {
        if (parent[hap]!=hap) {
            parent[hap] = find(parent[hap]);
        }
        return parent[hap];
    }

    /**
     * Merges the two sets in the partition that contain the specified
     * haplotypes if the two haplotypes are in different sets. No merging
     * is performed if the two haplotypes are in the same set when the
     * method is invoked.
     * @param x a haplotype index
     * @param y a haplotype index
     * @throws IndexOutOfBoundsException if {@code x < 0 || x >= this.nHaps()}
     * @throws IndexOutOfBoundsException if {@code y < 0 || y >= this.nHaps()}
     */
    public void union(int x, int y) {
        int xRoot = find(x);
        int yRoot = find(y);
        if (xRoot!=yRoot) {
            --nSets;
            if (rank[xRoot] <= rank[yRoot]) {
                if (rank[xRoot]==rank[yRoot]) {
                    ++rank[yRoot];
                }
                parent[xRoot] = yRoot;
            }
            else {
                parent[yRoot] = xRoot;
            }
        }
    }

    /**
     * Returns the number of sets in the partition
     * @return the number of sets in the partition
     */
    public int nSets() {
        return nSets;
    }

    /**
     * Returns a string representation of {@code this}. The string
     * representation does not end with a line separator string.
     * The exact details of the representation are unspecified and
     * subject to change.
     */
    @Override
    public String toString() {
        DecimalFormat df = new DecimalFormat("0.0000");
        int clustIndex = 0;
        Arrays.fill(rank, -1);   // temporarily use rank field to store cluster indices
        for (int j=0; j<parent.length; ++j) {
            parent[j] = find(j);
            if (rank[parent[j]]==-1) {
                rank[parent[j]] = clustIndex++;
            }
            rank[j] = rank[parent[j]];
        }
        StringBuilder sb = new StringBuilder(4*parent.length + 40);
        sb.append(ChromIds.instance().id(position.chrom()));
        sb.append(Const.tab);
        sb.append(position.pos());
        sb.append(Const.tab);
        sb.append(df.format(position.genPos()));
        for (int h=0; h<parent.length; h+=2) {
            sb.append(Const.tab);
            sb.append(rank[h]);
            sb.append(Const.phasedSep);
            sb.append(rank[h+1]);
        }
        Arrays.fill(rank, 1);   // set rank field after path compression
        return sb.toString();
    }

    /**
     * Writes {@code this.toString()} followed by the system line separator
     * string to the specified {@code PrintWriter}.
     * @param out the {@code PrintWriter} to which the string representation
     * will be written
     * @throws NullPointerException if {@code out == null}
     */
    public void write(PrintWriter out) {
        DecimalFormat df = new DecimalFormat("0.0000");
        int clustIndex = 0;
        Arrays.fill(rank, -1);   // temporarily use rank field to store cluster indices
        for (int j=0; j<parent.length; ++j) {
            parent[j] = find(j);
            if (rank[parent[j]]==-1) {
                rank[parent[j]] = clustIndex++;
            }
            rank[j] = rank[parent[j]];
        }
        out.print(ChromIds.instance().id(position.chrom()));
        out.print(Const.tab);
        out.print(position.pos());
        out.print(Const.tab);
        out.print(df.format(position.genPos()));
        for (int h=0; h<parent.length; h+=2) {
            out.print(Const.tab);
            out.print(rank[h]);
            out.print(Const.phasedSep);
            out.print(rank[h+1]);
        }
        out.println();
        Arrays.fill(rank, 1);   // set rank field after path compression
    }
}
