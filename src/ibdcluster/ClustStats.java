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

import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;
import java.util.concurrent.atomic.LongAdder;

/**
 * <p>Class {@code ClustStats} stores statistics from an ibd-cluster
 * analysis.</p>
 *
 * <p>Instances of class {@code ClustStats} are thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class ClustStats {

    private final AtomicInteger nSamples;
    private final AtomicLong nMarkers;
    private final AtomicLong nFilteredMarkers;
    private final AtomicLong nIbdSets;
    private final AtomicInteger nOutputPositions;
    private final LongAdder discordCnt;
    private final LongAdder totalCnt;

    /**
     * Constructs a new {@code ClustStats} instance.
     */
    public ClustStats(){
        this.nSamples = new AtomicInteger(0);
        this.nMarkers = new AtomicLong(0);
        this.nFilteredMarkers = new AtomicLong(0);
        this.nIbdSets = new AtomicLong(0);
        this.nOutputPositions = new AtomicInteger(0);
        this.discordCnt = new LongAdder();
        this.totalCnt = new LongAdder();
    }

    /**
     * Sets the number of samples to the specified value
     * @param nSamples the number of samples
     */
    public void setNSamples(int nSamples) {
        this.nSamples.set(nSamples);
    }

    /**
     * Returns the number of samples, or {@code 0} if {@code this.setNSamples()}
     * has not yet been invoked.
     * @return the number of samples
     */
    public int nSamples() {
        return nSamples.get();
    }

    /**
     * Adds the specified number of markers to the cumulative number
     * of markers.
     * @param cnt the number of markers to add
     */
    public void addMarkers(long cnt) {
        nMarkers.addAndGet(cnt);
    }

    /**
     * Adds the specified number of filtered markers to the cumulative
     * number of filtered markers.
     * @param cnt the number of filtered markers to add
     */
    public void addFilteredMarkers(long cnt) {
        nFilteredMarkers.addAndGet(cnt);
    }

    /**
     * Adds the specified number of IBD sets to the cumulative
     * number of IBD sets.
     * @param cnt the number of IBD sets to add
     */
    public void addIbdSets(long cnt) {
        nIbdSets.addAndGet(cnt);
    }

    /**
     * Adds the specified number of output positions to the cumulative
     * number of output positions.
     * @param cnt the number of output positions to add
     */
    public void addOutputPositions(int cnt) {
        nOutputPositions.addAndGet(cnt);
    }

    /**
     * Returns the cumulative number of markers.
     * @return the cumulative number of markers
     */
    public long nMarkers() {
        return nMarkers.get();
    }

    /**
     * Returns the cumulative number of filtered markers.
     * @return the cumulative number of filtered markers
     */
    public long nFilteredMarkers() {
        return nFilteredMarkers.get();
    }

    /**
     * Returns the number of IBD sets.
     * @return the number of IBD sets
     */
    public long nIbdSets() {
        return nIbdSets.get();
    }

    /**
     * Returns the number of output positions.
     * @return the number of output positions
     */
    public int nOutputPositions() {
        return nOutputPositions.get();
    }

    /**
     * Adds the specified number of discordant alleles and number
     * of examined alleles to the cumulative number of discordant
     * alleles and number of examined alleles.
     * @param discordant the number of discordant alleles
     * @param total the total number of examined alleles
     * @throws IllegalArgumentException if {@code discordant > total}
     */
    public void updateDiscordRate(int discordant, int total) {
        if (discordant > total) {
            throw new IllegalArgumentException(discordant + ">" + total);
        }
        discordCnt.add(discordant);
        totalCnt.add(total);
    }

    /**
     * Returns the IBD segment allele discordance rate. Invocation in the
     * absence of concurrent updates returns an accurate result, but
     * concurrent updates that occur while the count is being calculated might
     * not be incorporated.
     * @return the IBD segment allele discordance rate
     */
    public float discordRate() {
        long num = discordCnt.sum();
        long den = totalCnt.sum();
        return (float) ((double) num / den);
    }
}
