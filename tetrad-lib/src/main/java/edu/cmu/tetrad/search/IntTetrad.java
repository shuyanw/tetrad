///////////////////////////////////////////////////////////////////////////////
// For information as to what this class does, see the Javadoc, below.       //
// Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006,       //
// 2007, 2008, 2009, 2010, 2014, 2015 by Peter Spirtes, Richard Scheines, Joseph   //
// Ramsey, and Clark Glymour.                                                //
//                                                                           //
// This program is free software; you can redistribute it and/or modify      //
// it under the terms of the GNU General Public License as published by      //
// the Free Software Foundation; either version 2 of the License, or         //
// (at your option) any later version.                                       //
//                                                                           //
// This program is distributed in the hope that it will be useful,           //
// but WITHOUT ANY WARRANTY; without even the implied warranty of            //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             //
// GNU General Public License for more details.                              //
//                                                                           //
// You should have received a copy of the GNU General Public License         //
// along with this program; if not, write to the Free Software               //
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA //
///////////////////////////////////////////////////////////////////////////////

package edu.cmu.tetrad.search;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * Represents a tetrad of variables.
 */
public class IntTetrad {
    private int i;
    private int j;
    private int k;
    private int l;
    private double pValue;

    public IntTetrad(int i, int j, int k, int l) {
        testDistinctness(i, j, k, l);
        this.i = i;
        this.j = j;
        this.k = k;
        this.l = l;
        this.pValue = Double.NaN;
    }

    public IntTetrad(int i, int j, int k, int l, double pValue) {
        testDistinctness(i, j, k, l);
        this.i = i;
        this.j = j;
        this.k = k;
        this.l = l;
        this.pValue = pValue;
    }

    private void testDistinctness(int i, int j, int k, int l) {
        if (i == j || i == k || i == l) {
            throw new IllegalArgumentException("ints not distinct.");
        }

        if (j == k || j == l) {
            throw new IllegalArgumentException("ints not distinct.");
        }

        if (k == l) {
            throw new IllegalArgumentException("ints not distinct.");
        }
    }

    public int getI() {
        return i;
    }

    public int getJ() {
        return j;
    }

    public int getK() {
        return k;
    }

    public int getL() {
        return l;
    }

    public int hashCode() {
        return i * j + k * l;
    }

    public boolean equals(Object o) {
        if (o == null) return false;
        IntTetrad tetrad = (IntTetrad) o;
        return (i == tetrad.i && j == tetrad.j && k == tetrad.k && l == tetrad.l)
                || (i == tetrad.j && j == tetrad.i && k == tetrad.k && l == tetrad.l)
                || (i == tetrad.i && j == tetrad.j && k == tetrad.l && l == tetrad.k)
                || (i == tetrad.j && j == tetrad.i && k == tetrad.l && l == tetrad.k)
                || (i == tetrad.k && j == tetrad.l && k == tetrad.i && l == tetrad.j)
                || (i == tetrad.k && j == tetrad.l && k == tetrad.j && l == tetrad.i)
                || (i == tetrad.l && j == tetrad.k && k == tetrad.i && l == tetrad.j)
                || (i == tetrad.l && j == tetrad.k && k == tetrad.j && l == tetrad.i);
    }

    public String toString() {
        NumberFormat nf = new DecimalFormat("0.0000");

        if (Double.isNaN(pValue)) {
            return "s(" + i + "," + j + ")*s(" + k + "," + l + ")-s(" + i + "," + k + ")*s(" + j + "," + l + ")";
        } else {
            return "<" + i + ", " + j + ", " + k + ", " + l + ", " + nf.format(pValue) + ">";
        }
    }

    public double getPValue() {
        return pValue;
    }

    public List<Integer> getInts() {
        List<Integer> ints = new ArrayList<>();
        ints.add(i);
        ints.add(j);
        ints.add(k);
        ints.add(l);
        return ints;
    }
}


