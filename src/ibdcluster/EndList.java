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

import java.util.Arrays;
import java.util.NoSuchElementException;

/**
 * <p>Class {@code EndList} represents the elements at the end of a list.</p>
 *
 * <p>Instances of class {@code EndList} are not thread-safe.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class EndList {

    private int capacity;
    private int mask;
    private Partition[] list;

    private int start;
    private int end;

    /**
     * Constructs a new {@code EndList} object with the specified minimum
     * capacity.  The capacity of the {@code EndList} is guaranteed to
     * be no larger than the smallest positive integral power of 2 that is
     * greater than or equal to {@code minCapacity}.
     * @param minCapacity the minimum capacity
     * @throws NegativeArraySizeException if {@code minCapacity > (1<<30)}
     */
    public EndList(int minCapacity) {
        if (minCapacity < 1) {
            minCapacity = 1;
        }
        this.capacity = Integer.highestOneBit(minCapacity);
        if (capacity<minCapacity) {
            capacity<<=1;
        }
        this.mask = capacity - 1;
        this.list = new Partition[capacity];
        this.start = 0;
        this.end = 0;
    }

    /**
     * Adds the specified object to the end of this list.
     * @param p the partition to be added
     * @throws NegativeArraySizeException if
     * {@code (this.end() - this.start() + 1) == (1<<30)}
     */
    public void add(Partition p) {
        if ((end - start) == mask) {
            int oldMask = mask;
            Partition[] oldList = list;
            capacity <<=1;
            mask = capacity - 1;
            list = new Partition[capacity];
            for (int j=start; j<end; ++j) {
                list[j & mask] = oldList[j & oldMask];
            }
        }
        list[end++ & mask] = p;
    }

    /**
     * Returns the specified object.
     * @param index an index
     * @return the specified object
     * @throws IndexOutOfBoundsException if
     * {@code index < this.start() || index >= this.end()}
     */
    public Partition get(int index) {
        if (index < start || index >= end) {
            throw new IndexOutOfBoundsException(String.valueOf(index));
        }
        return list[index & mask];
    }

    /**
     * Removes and returns the partition with index {@code this.start()}
     * and increments the start index.
     * @return the partition with index {@code this.start()}
     * @throws NoSuchElementException if {@code this.start() == this.end()}
     */
    public Partition removeStart() {
        if (start >= end) {
            throw new NoSuchElementException();
        }
        Partition p = list[start & mask];
        list[start++ & mask] = null;
        return p;
    }

    /**
     * Removes all elements from the list
     */
    public void clear() {
        Arrays.fill(list, null);
        this.start = 0;
        this.end = 0;
    }

    /**
     * Returns the index of the first remaining element modulo
     * {@code 0x7fffffff}. Returns {@code this.end()} if there are no
     * remaining elements.
     * @return the index of the first remaining element modulo
     * {@code 0x7fffffff}
     */
    public int start() {
        return start & 0x7f_ff_ff_ff;
    }

    /**
     * Returns the cumulative number of elements that have been stored in
     * this list modulo {@code 0x7fffffff}.
     * @return the cumulative number of elements that have been stored in
     * this list modulo {@code 0x7fffffff}
     */
    public int end() {
        return end & 0x7f_ff_ff_ff;
    }
}
