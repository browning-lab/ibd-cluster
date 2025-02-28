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
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * <p>Class {@code Partition} implements a disjoint union data
 * structure that stores IBD clusters.</p>
 *
 * <p>Instances of class {@code Partition} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class Partition {

    private final String chrom;
    private final int pos;
    private final double cmPos;
    private final int[] parent;
    private final int[] rank;
    private int nSets;

    /**
     * Constructs a new {@code Partition} instance for the specified data.
     * Each haplotype is contained in a singleton set within the partition.
     *
     * @param chrom the chromosome
     * @param pos a base position
     * @param cmPos the genetic map cM position
     * @param nHaps the number of elements that are partitioned
     * @throws NegativeArraySizeException if {@code nHaps < 0}
     * @throws NullPointerException if {@code chrom == null}
     */
    public Partition(String chrom, int pos, double cmPos, int nHaps) {
        if (chrom==null) {
            throw new NullPointerException("chrom==null");
        }
        this.chrom = chrom;
        this.pos = pos;
        this.cmPos = cmPos;
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
     * Returns the chromosome.
     * @return the chromosome
     */
    public String chrom() {
        return chrom;
    }

    /**
     * Returns the base position.
     * @return the base position
     */
    public int pos() {
        return pos;
    }

    /**
     * Returns the genetic map cM position.
     * @return the genetic map cM position
     */
    public double cmPos() {
        return cmPos;
    }

    /**
     * Returns the representative haplotype of the set in the partition
     * that contains the specified haplotype. This method performs path
     * compression.
     * @param hap a haplotype index
     * @return the representative haplotype of the set in the partition
     * that contains the specified haplotype
     * @throws IndexOutOfBoundsException if {@code x < 0 || x >= this.nHaps()}
     */
    public int find(int hap) {
        if (parent[hap]!=hap) {
            parent[hap] = find(parent[hap]);
        }
        return parent[hap];
    }

    /**
     * Merges the sets in the partition that contain the specified
     * haplotypes if the two haplotypes are in different sets.
     * @param x a haplotype index
     * @param y a haplotype index
     * @throws IndexOutOfBoundsException if {@code (x < 0 || x >= this.nHaps())}
     * @throws IndexOutOfBoundsException if {@code (y < 0 || y >= this.nHaps())}
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
     * Appends a string representation of {@code this} to the
     * specified {@code StringBuilder}. The details of the string
     * representation are unspecified and subject to change.
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
        sb.append(chrom);
        sb.append(Const.tab);
        sb.append(pos);
        sb.append(Const.tab);
        sb.append(df.format(cmPos));
        for (int h=0; h<parent.length; h+=2) {
            sb.append(Const.tab);
            sb.append(rank[h]);
            sb.append(Const.phasedSep);
            sb.append(rank[h+1]);
        }
        Arrays.fill(rank, 1);   // set rank field after path compression
        sb.append(Const.nl);
        return sb.toString();
    }

    /**
     * Writes a string representation of {@code this} to the specified
     * {@code PrintWriter}. The details of the string representation 
     * are unspecified and subject to change.
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
        out.print(chrom);
        out.print(Const.tab);
        out.print(pos);
        out.print(Const.tab);
        out.print(df.format(cmPos));
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
