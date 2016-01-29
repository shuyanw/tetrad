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

import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.ProbUtils;
import org.apache.commons.math3.special.Gamma;

import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Calculates the BDeu score.
 */
public class BDeuScore2 implements LocalDiscreteScore, GesScore {
    private final int[] reuse;
    private List<Node> variables;
    private List<DiscreteVariable> discreteVariables;
    private int[][] data;
    //    private final LocalScoreCache localScoreCache = new LocalScoreCache();
    private int sampleSize;

    private double samplePrior = 1;
    private double structurePrior = 1;

    private int[] numCategories;

    private double lastBumpThreshold = 0.0;

    private VariableBox[][] boxes;

    private int[] save1;
    private int[] save2;


    public BDeuScore2(DataSet dataSet) {
        if (dataSet == null) {
            throw new NullPointerException();
        }

        if (dataSet instanceof BoxDataSet) {
            DataBox dataBox = ((BoxDataSet) dataSet).getDataBox();

            this.variables = dataSet.getVariables();
            this.discreteVariables = discreteVariables(dataSet.getVariables());

            if (!(((BoxDataSet) dataSet).getDataBox() instanceof VerticalIntDataBox)) {
                throw new IllegalArgumentException();
            }

            VerticalIntDataBox box = (VerticalIntDataBox) dataBox;

            data = box.getVariableVectors();
            this.sampleSize = dataSet.getNumRows();
        } else {
            data = new int[dataSet.getNumColumns()][];
            this.variables = dataSet.getVariables();
            this.discreteVariables = discreteVariables(dataSet.getVariables());

            for (int j = 0; j < dataSet.getNumColumns(); j++) {
                data[j] = new int[dataSet.getNumRows()];

                for (int i = 0; i < dataSet.getNumRows(); i++) {
                    data[j][i] = dataSet.getInt(i, j);
                }
            }

            this.sampleSize = dataSet.getNumRows();
        }

        save1 = new int[data[0].length + 1];
        save2 = new int[data[0].length + 1];


        final List<Node> variables = dataSet.getVariables();
        numCategories = new int[variables.size()];
        for (int i = 0; i < variables.size(); i++) {
            numCategories[i] = (int) (getVariable(i)).getNumCategories();
        }

        reuse = new int[data[0].length];

        countData();
    }

    private List<DiscreteVariable> discreteVariables(List<Node> variables) {
        List<DiscreteVariable> discreteVariables = new ArrayList<>();

        for (Node node : variables) {
            discreteVariables.add((DiscreteVariable) node);
        }

        return discreteVariables;
    }

    private DiscreteVariable getVariable(int i) {
        return (DiscreteVariable) variables.get(i);
    }

    public double localScore(int node, int parents[]) {

        // Number of categories for node.
        int r = numCategories[node];

        // Numbers of categories of parents.
        int[] dims = new int[parents.length];

        for (int p = 0; p < parents.length; p++) {
            dims[p] = numCategories[parents[p]];
        }

        // Number of parent states.
        int q = 1;

        for (int p = 0; p < parents.length; p++) {
            q *= dims[p];
        }

        // Conditional cell coefs of data for node given parents(node).
        int n_jk[][] = new int[q][r];
        int n_j[] = new int[q];

        for (int _q = 0; _q < q; _q++) {
            int[] parentValues = getParentValues(_q, dims);

            for (int _r = 0; _r < r; _r++) {
                int[] intersection = null;

                if (parentValues.length == 0) {
                    int count = (int) this.boxes[node][node].getList(_r, _r).length;
                    n_jk[_q][_r] = count;
                    n_j[_q] += count;
                } else {
                    for (int i = 0; i < parentValues.length; i++) {
                        // Need the intersection of Pi = _pi x C = _c for each i.
                        int Pi = parents[i];
                        int pi = parentValues[i];

                        this.boxes[Pi][node].getList(pi, _r);

                        if (intersection == null) {
                            intersection = this.boxes[Pi][node].getList(pi, _r);
                        } else {
                            intersection = findIntersection(intersection, this.boxes[Pi][node].getList(pi, _r));
                        }
                    }

                    if (intersection == null) {
                        throw new NullPointerException();
                    }

                    n_jk[_q][_r] = (int) intersection.length;
                    n_j[_q] += intersection.length;
                }
            }
        }

        //Finally, compute the score
        double score = 0.0;

        score += (r - 1) * q * Math.log(getStructurePrior());

        final double cellPrior = getSamplePrior() / (r * q);
        final double rowPrior = getSamplePrior() / q;

        for (int j = 0; j < q; j++) {
            score -= Gamma.logGamma(rowPrior + n_j[j]);

            for (int k = 0; k < r; k++) {
                score += Gamma.logGamma(cellPrior + n_jk[j][k]);
            }
        }

        score += q * Gamma.logGamma(rowPrior);
        score -= r * q * Gamma.logGamma(cellPrior);

        lastBumpThreshold = ((r - 1) * q * Math.log(getStructurePrior()));

        return score;
    }

    @Override
    public double localScoreDiff(int i, int[] parents, int extra) {
        return localScore(i, append(parents, extra)) - localScore(i, parents);
    }

    int[] append(int[] parents, int extra) {
        int[] all = new int[parents.length + 1];
        System.arraycopy(parents, 0, all, 0, parents.length);
        all[parents.length] = extra;
        return all;
    }

    @Override
    public double localScore(int node, int parent) {
        return localScore(node, new int[]{parent});
    }

    @Override
    public double localScore(int node) {
        return localScore(node, new int[0]);
    }

    @Override
    public List<Node> getVariables() {
        return this.variables;
    }

    @Override
    public int getSampleSize() {
        return sampleSize;
    }

    /**
     * Must be called directly after the corresponding scoring call.
     */
    public boolean isEffectEdge(double bump) {
        return bump > lastBumpThreshold;
    }

    @Override
    public boolean isDiscrete() {
        return true;
    }

    @Override
    public DataSet getDataSet() {
        throw new UnsupportedOperationException();
    }

    private static int getRowIndex(int[] dim, int[] values) {
        int rowIndex = 0;
        for (int i = 0; i < dim.length; i++) {
            rowIndex *= dim[i];
            rowIndex += values[i];
        }
        return rowIndex;
    }

    public int[] getParentValues(int rowIndex, int[] parentDims) {
        int[] values = new int[parentDims.length];

        for (int i = parentDims.length - 1; i >= 0; i--) {
            values[i] = rowIndex % parentDims[i];
            rowIndex /= parentDims[i];
        }

        return values;
    }

    /**
     * Finds the intersection of l1 and l2, placing the result in intersect. It is assumed that l1 and l2 are
     * sorted. Intersect is reused, is why it's a parameter.
     */
    private int[] findIntersection(int[] l1, int[] l2) {
        int i = 0;
        int j = 0;
        List<Integer> intersection = new ArrayList<>();
        boolean column = true;

        while (i < l1.length && j < l2.length) {
            if (l1[i] == l2[j]) {
                intersection.add(l1[i]);
                i++;
                j++;
            } else if (column ? l1[i] > l2[j] : l2[j] > l1[i]) {
                column = !column;
            } else {
                if (column) i++;
                else j++;
            }
        }

        int[] ret = new int[intersection.size()];
        for (int k = 0; k < intersection.size(); k++) ret[k] = intersection.get(k);
        return ret;
    }

    /**
     * Finds the intersection of A and B, placing the result in intersect. It is assumed that A and B are
     * sorted. Intersect is reused, is why it's a parameter.
     */
    private int[] findIntersection2(int[] A, int[] B) {
        if (B.length > A.length) {
            int[] temp = A;
            A = B;
            B = temp;
        }

        int n = A.length;
        int m = B.length;
        int y = 0;

        int[] intersect = new int[Math.min(m, n)];

        int low = 1;

        for (int i = 0; i < m; i++) {
            int diff = 1;

            while (low + diff <= n && A[low + diff - 1] < B[i]) {
                diff *= 2;
            }

            int high = Math.min(n, low + diff);

            // Returns -1 if not found.
            int k = binarySearch(A, low, high, B[i]);

            if (k != -1) {
                intersect[y++] = A[k];
            }

            low = k + 1;
        }

        int[] ret = new int[y];
        System.arraycopy(intersect, 0, ret, 0, y);
        return ret;
    }

    /**
     * Searches a range of the sorted array of for the specified value using the
     * binary search algorithm. Bounds are not checked.
     *
     * @param a         the array to be searched
     * @param fromIndex the index of the first element (inclusive) to be
     *                  searched
     * @param toIndex   the index of the last element (exclusive) to be searched
     * @param key       the value to be searched for
     * @return index of the search key, if it is contained in the array
     * within the specified range;
     * otherwise, -1 if the key is not found.
     * @since 1.6
     */
    private static int binarySearch(int[] a, int fromIndex, int toIndex,
                                    int key) {
        int low = fromIndex;
        int high = toIndex - 1;

        while (low <= high) {
            int mid = (low + high) >>> 1;
            int midVal = a[mid];

            if (midVal < key)
                low = mid + 1;
            else if (midVal > key)
                high = mid - 1;
            else
                return mid; // key found
        }
//        return (low + 1);  // key not found.
        return -1;
    }

    private int findIntersection2(int[] l1, int[] l2, int[] intersection) {
        int i = 0;
        int j = 0;
        int k = 0;
        boolean column = true;

        while (l1[i] != -99 && l2[j] != -99) {
            if (l1[i] == l2[j]) {
                intersection[k] = l1[i];
                i++;
                j++;
                k++;
            } else if (column ? l1[i] > l2[j] : l2[j] > l1[i]) {
                column = !column;
            } else {
                if (column) i++;
                else j++;
            }
        }

        intersection[k] = -99;

        return k;
    }

    public double getStructurePrior() {
        return structurePrior;
    }

    public double getSamplePrior() {
        return samplePrior;
    }

    public void setStructurePrior(double structurePrior) {
        this.structurePrior = structurePrior;
    }

    public void setSamplePrior(double samplePrior) {
        this.samplePrior = samplePrior;
    }

    private void countData() {
        boxes = new VariableBox[discreteVariables.size()][discreteVariables.size()];

        for (int i = 0; i < discreteVariables.size(); i++) {
            System.out.println("counting i = " + i);

            for (int j = i; j < discreteVariables.size(); j++) {
                DiscreteVariable v1 = discreteVariables.get(i);
                DiscreteVariable v2 = discreteVariables.get(j);

                VariableBox variableBox = new VariableBox(v1, v2, data[i], data[j]);

                if (i == j) {
                    boxes[i][j] = variableBox;
                } else {
                    boxes[i][j] = variableBox;
                    boxes[j][i] = boxes[i][j].reverse();
                }
            }
        }
    }

    private class VariableBox {
        private int numCategories1;
        private int numCategories2;

        private int[][][] lists;

        private VariableBox() {

        }

        public VariableBox(DiscreteVariable var1, DiscreteVariable var2, int[] data1, int[] data2) {
            this.numCategories1 = var1.getNumCategories();
            this.numCategories2 = var2.getNumCategories();

            lists = new int[numCategories1][numCategories2][data1.length];

            for (int i = 0; i < numCategories1; i++) {
                for (int j = 0; j < numCategories2; j++) {
                    lists[i][j] = new int[data1.length];
                }
            }

            int[][] bounds = new int[numCategories1][numCategories2];

            for (int i = 0; i < data1.length; i++) {
                int d1 = data1[i];
                int d2 = data2[i];
                lists[d1][d2][bounds[d1][d2]] = i;
                bounds[d1][d2]++;
            }

            for (int d1 = 0; d1 < numCategories1; d1++) {
                for (int d2 = 0; d2 < numCategories2; d2++) {
                    int[] list = lists[d1][d2];
                    int[] interList = new int[bounds[d1][d2]];
                    System.arraycopy(list, 0, interList, 0, bounds[d1][d2]);
                    lists[d1][d2] = interList;
                }
            }
        }

        public int[] getList(int d1, int d2) {
            return lists[d1][d2];
        }

        public VariableBox reverse() {
            VariableBox box2 = new VariableBox();

            box2.numCategories1 = numCategories2;
            box2.numCategories2 = numCategories1;
            box2.lists = transpose(lists);

            return box2;
        }

        private int[][][] transpose(int[][][] lists) {
            int[][][] lists2 = new int[lists[0].length][lists.length][];

            for (int i = 0; i < lists[0].length; i++) {
                for (int j = 0; j
                        < lists.length; j++) {
                    lists2[i][j] = lists[j][i];
                }
            }

            return lists2;
        }
    }
}



