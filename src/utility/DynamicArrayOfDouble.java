/*
 * Copyright (c) 2015, Asad
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package utility;

/**
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class DynamicArrayOfDouble {

    private double[] data;
    int position = 0;

    /**
     * Dynamic Double Array
     */
    public DynamicArrayOfDouble() {
        data = new double[1];
    }

    /**
     *
     * @param position
     * @return
     */
    public double get(int position) {
        if (position >= data.length) {
            return 0;
        } else {
            return data[position];
        }
    }

    /**
     *
     * @param position
     * @param value
     */
    public void put(int position, double value) {
        this.position = position;
        if (position >= data.length) {
            int newSize = 2 * data.length;
            if (position >= newSize) {
                newSize = 2 * position;
            }
            double[] newData = new double[newSize];
            System.arraycopy(data, 0, newData, 0, data.length);
            data = newData;
//            System.out.println("Size of dynamic array increased to " + newSize);
        }
        data[position] = value;
    }

    /**
     *
     * @return length of the array
     */
    public int getSize() {
        return (position + 1);
    }

    /**
     *
     * @return length of the array
     */
    public double[] getArray() {
        double[] newData = new double[position + 1];
        System.arraycopy(data, 0, newData, 0, getSize());
        return newData;
    }
}
