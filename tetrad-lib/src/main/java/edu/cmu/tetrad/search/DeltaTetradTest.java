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
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;

import java.util.*;

/**
 * Implements a test for simultaneously zero tetrads in Bollen, K. (1990). "Outlier screening and distribution-free test
 * for vanishing tetrads." Sociological Methods and Research 19, 80-92 and Bollen and Ting, Confirmatory Tetrad
 * Analysis.
 *
 * @author Joseph Ramsey
 */
public class DeltaTetradTest {
    private int N;
    private ICovarianceMatrix cov;
    private final CorrelationMatrix corr;
    private int df;
    private double chisq;
    private double[][] centeredData;
    private double[][] standardizedData;
    private List<Node> variables;
    private Map<Node, Integer> variablesHash;


    // As input we require a data set and a list of non-redundant Tetrads.

    // Need a method to remove Tetrads from the input list until what's left is
    // non-redundant. Or at least to check for non-redundancy. Maybe simply
    // checking to see whether a matrix exception is thrown is sufficient.
    // Peter suggests looking at Modern Factor Analysis by Harmon, at triplets.

    /**
     * Constructs a test using a given standardizedData set. If a standardizedData set is provided (that is, a tabular standardizedData set), fourth moment
     * statistics can be calculated (p. 160); otherwise, it must be assumed that the standardizedData are multivariate Gaussian.
     */
    public DeltaTetradTest(DataSet dataSet) {
        System.out.println(dataSet);

        if (dataSet == null) {
            throw new NullPointerException();
        }

        if (!dataSet.isContinuous()) {
            throw new IllegalArgumentException();
        }

        this.cov = new CovarianceMatrix(dataSet);
        this.corr = new CorrelationMatrix(dataSet);

        this.centeredData = DataUtils.centerData(dataSet.getDoubleData()).toArray();
        this.standardizedData = DataUtils.standardizeData(dataSet.getDoubleData()).transpose().toArray();
        this.N = dataSet.getNumRows();
        this.variables = dataSet.getVariables();

        this.variablesHash = new HashMap<>();

        for (int i = 0; i < variables.size(); i++) {
            variablesHash.put(variables.get(i), i);
        }
    }

    /**
     * Constructs a test using the given covariance matrix. Fourth moment statistics are not caculated; it is assumed
     * that the standardizedData are distributed as multivariate Gaussian.
     */
    public  DeltaTetradTest(ICovarianceMatrix cov) {
        if (cov == null) {
            throw new NullPointerException();
        }

        this.cov = cov;
        this.corr = new CorrelationMatrix(cov.copy());
        this.N = cov.getSampleSize();
        this.variables = cov.getVariables();

        this.variablesHash = new HashMap<>();

        for (int i = 0; i < variables.size(); i++) {
            variablesHash.put(variables.get(i), i);
        }
    }

//    private void initializeForthMomentMatrix(List<Node> variables) {
//        int n = variables.size();
//        fourthMoment = new double[n][n][n][n];
//    }

    /**
     * Takes a list of tetrads for the given data set and returns the chi square value for the test. We assume that the
     * tetrads are non-redundant; if not, a matrix exception will be thrown.
     * <p>
     * Calculates the T statistic (Bollen and Ting, p. 161). This is significant if tests as significant using the Chi
     * Square distribution with degrees of freedom equal to the number of nonredundant tetrads tested.
     */
    public double calcChiSquare(Tetrad... tetrads) {
        this.df = tetrads.length;

        // Need a list of symbolic covariances--i.e. covariances that appear in tetrads.
        Set<Sigma> boldSigmaSet = new LinkedHashSet<>();
        List<Sigma> boldSigma = new ArrayList<>();

        for (Tetrad tetrad : tetrads) {
            boldSigmaSet.add(new Sigma(tetrad.getI(), tetrad.getK()));
            boldSigmaSet.add(new Sigma(tetrad.getI(), tetrad.getL()));
            boldSigmaSet.add(new Sigma(tetrad.getJ(), tetrad.getK()));
            boldSigmaSet.add(new Sigma(tetrad.getJ(), tetrad.getL()));
        }

        for (Sigma sigma : boldSigmaSet) {
            boldSigma.add(sigma);
        }

        // Need a matrix of variances and covariances of sample covariances.
        TetradMatrix sigma_ss = new TetradMatrix(boldSigma.size(), boldSigma.size());

        for (int i = 0; i < boldSigma.size(); i++) {
            for (int j = 0; j < boldSigma.size(); j++) {
                Sigma sigmaef = boldSigma.get(i);
                Sigma sigmagh = boldSigma.get(j);

                Node e = sigmaef.getA();
                Node f = sigmaef.getB();
                Node g = sigmagh.getA();
                Node h = sigmagh.getB();

                double _ss;

                if (centeredData != null) {
                    _ss = s(e, f, g, h) - s(e, f) * s(g, h);

                    // General.
//                   _ss = r(e, f, g, h) + 0.25 * r(e, f) * r(g, h) *
//                            (r(e, e, g, g) * r(f, f, g, g) + r(e, e, h, h) + r(f, f, h, h))
//                            - 0.5 * r(e, f) * (r(e, e, g, h) + r(f, f, g, h))
//                            - 0.5 * r(g, h) * (r(e, f, g, g) + r(e, f, h, h));

                } else if (cov != null) {
                    if (cov instanceof CorrelationMatrix) {
//                        Assumes multinormality. Using formula 23. (Not implementing formula 22 because that case
//                        does not come up.)
                        _ss = 0.5 * (r(e, f) * r(g, h))
                                * (r(e, g) * r(e, g) + r(e, h) * r(e, h) + r(f, g) * r(f, g) + r(f, h) * r(f, h))
                                + r(e, g) * r(f, h) + r(e, h) * r(f, g)
                                - r(e, f) * (r(f, g) * r(f, h) + r(e, g) * r(e, h))
                                - r(g, h) * (r(f, g) * r(e, g) + r(f, h) * r(e, h));

                    } else {
                        // Assumes multinormality--see p. 160.
                        _ss = s(e, g) * s(f, h) + s(e, h) * s(f, g);
                    }
                } else {
                    throw new IllegalStateException();
                }

                sigma_ss.set(i, j, _ss);
            }
        }

        // Need a matrix of of population estimates of partial derivatives of tetrads
        // with respect to covariances in boldSigma.w
        TetradMatrix del = new TetradMatrix(boldSigma.size(), tetrads.length);

        for (int i = 0; i < boldSigma.size(); i++) {
            for (int j = 0; j < tetrads.length; j++) {
                Sigma sigma = boldSigma.get(i);
                Tetrad tetrad = tetrads[j];

                Node e = tetrad.getI();
                Node f = tetrad.getJ();
                Node g = tetrad.getK();
                Node h = tetrad.getL();

                double derivative = getDerivative(e, f, g, h, sigma.getA(), sigma.getB());
                del.set(i, j, derivative);
            }
        }

        // Need a vector of population estimates of the tetrads.
        TetradMatrix t = new TetradMatrix(tetrads.length, 1);

        for (int i = 0; i < tetrads.length; i++) {
            Tetrad tetrad = tetrads[i];

            Node e = tetrad.getI();
            Node f = tetrad.getJ();
            Node g = tetrad.getK();
            Node h = tetrad.getL();

            double d1 = s(e, f);
            double d2 = s(g, h);
            double d3 = s(e, g);
            double d4 = s(f, h);

            double value = d1 * d2 - d3 * d4;
            t.set(i, 0, value);
        }

        // Now multiply to get Sigma_tt
        TetradMatrix w1 = del.transpose().times(sigma_ss);
        TetradMatrix sigma_tt = w1.times(del);

        // And now invert and multiply to get T.
        TetradMatrix v0 = sigma_tt.inverse();
        TetradMatrix v1 = t.transpose().times(v0);
        TetradMatrix v2 = v1.times(t);
        double chisq = N * v2.get(0, 0);

        this.chisq = chisq;

        return chisq;
    }

    /**
     * @return the p value for the most recent test.
     */
    public double getPValue() {
        double cdf = new ChiSquaredDistribution(this.df).cumulativeProbability(this.chisq);
        return 1.0 - cdf;
    }

    public double getPValue(Tetrad... tetrads) {
        calcChiSquare(tetrads);
        return getPValue();
    }

    private double r(Node a, Node b) {
        int i = variablesHash.get(a);
        int j = variablesHash.get(b);

        if (cov != null) {
            return corr.getValue(i, j);
        } else {
            throw new IllegalStateException();
        }
    }

    private double s(Node a, Node b) {
        int i = variablesHash.get(a);
        int j = variablesHash.get(b);

        if (cov != null) {
            return cov.getValue(i, j);
        }

        throw new IllegalStateException();
    }

    // Assumes standardizedData are mean-centered.
    private double r(Node x, Node y, Node z, Node w) {
        int i = variablesHash.get(x);
        int j = variablesHash.get(y);
        int k = variablesHash.get(z);
        int l = variablesHash.get(w);

        double sxyzw = 0.0;

        double[] _x = standardizedData[i];
        double[] _y = standardizedData[j];
        double[] _z = standardizedData[k];
        double[] _w = standardizedData[l];

        int N = _x.length;

        for (int g = 0; g < N; g++) {
            sxyzw += _x[g] * _y[g] * _z[g] * _w[g];
        }

        return (1.0 / N) * sxyzw;
    }

    private double s(Node x, Node y, Node z, Node w) {
        int i = variablesHash.get(x);
        int j = variablesHash.get(y);
        int k = variablesHash.get(z);
        int l = variablesHash.get(w);

        double sxyzw = 0.0;

        double[] _x = centeredData[i];
        double[] _y = centeredData[j];
        double[] _z = centeredData[k];
        double[] _w = centeredData[l];

        int N = _x.length;

        for (int g = 0; g < N; g++) {
            sxyzw += _x[g] * _y[g] * _z[g] * _w[g];
        }

        return (1.0 / N) * sxyzw;
    }

    private double getDerivative(Node node1, Node node2, Node node3, Node node4, Node a, Node b) {
        if (node1 == a && node2 == b) {
            return s(node3, node4);
        }

        if (node1 == b && node2 == a) {
            return s(node3, node4);
        }

        if (node3 == a && node4 == b) {
            return s(node1, node2);
        }

        if (node3 == b && node4 == a) {
            return s(node1, node2);
        }

        if (node1 == a && node3 == b) {
            return -s(node2, node4);
        }

        if (node1 == b && node3 == a) {
            return -s(node2, node4);
        }

        if (node2 == a && node4 == b) {
            return -s(node1, node3);
        }

        if (node2 == b && node4 == a) {
            return -s(node1, node3);
        }

        return 0.0;
    }

    private static class Sigma {
        private Node a;
        private Node b;

        public Sigma(Node a, Node b) {
            this.a = a;
            this.b = b;
        }

        public Node getA() {
            return a;
        }

        public Node getB() {
            return b;
        }

        public boolean equals(Object o) {
            if (!(o instanceof Sigma)) {
                throw new IllegalArgumentException();
            }

            Sigma _o = (Sigma) o;
            return (_o.getA().equals(getA()) && _o.getB().equals(getB())) || (_o.getB().equals(getA()) && _o.getA().equals(getB()));
        }

        public int hashCode() {
            return a.hashCode() + b.hashCode();
        }

        public String toString() {
            return "Sigma(" + getA() + ", " + getB() + ")";
        }
    }
}



