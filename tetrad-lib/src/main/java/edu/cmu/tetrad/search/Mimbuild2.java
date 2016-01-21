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

import edu.cmu.tetrad.data.CovarianceMatrix;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;

import java.util.ArrayList;
import java.util.List;

/**
 * An implemetation of Mimbuild based on the Fgsl score. The search will attempt a GES search first and if that
 * throws and exception then a CPC search. The penalty discount parameter is for the GES search; the alpha
 * value is for the CPC search. Or you can just grab the latent covariance matrix and run whatever search you
 * want to. (I don't know why GES sometimes fails, it is a mystery.)
 * </p>
 * Uses a different (better) algorithm from Mimbuild. Preferable.
 *
 * @author Joseph Ramsey
 */
public class Mimbuild2 {

    /**
     *
     */
    private Graph structureGraph;

    /**
     * The alpha level used for CPC
     */
    private double alpha = 0.001;

    /**
     * Background knowledge for CPC.
     */
    private IKnowledge knowledge = new Knowledge2();

    /**
     * The estimated covariance matrix over the latents.
     */
    private ICovarianceMatrix latentsCov;

    /**
     * The p value of the optimization.
     */
    private double pValue;
    private List<Node> latents;
    private int minClusterSize = 4;
    private boolean searchForMax = false;
    private List<List<Node>> clustering;
    double threshold = 1e-7;

    public Mimbuild2() {
    }

    //=================================== PUBLIC METHODS =========================================//

    public Graph search(List<List<Node>> clustering, List<String> latentNames, ICovarianceMatrix measuresCov) {
        int N = measuresCov.getSampleSize();
        List<String> _latentNames = new ArrayList<>(latentNames);
        clustering = new ArrayList<>(clustering);

        List<Node> latents = defineLatents(_latentNames);
        this.latents = latents;

        // This removes the small clusters and their names.
        removeSmallClusters(latents, clustering, getMinClusterSize());

        if (searchForMax) {
            double bestMaxP = -1;
            List<List<Node>> bestClustering = null;
            Graph bestSructureGraph = null;
            List<Node> bestLatents = null;
            ICovarianceMatrix bestLatentsCov = null;

            int min = 1;
            int max = 5;

            for (int i = min; i <= Math.min(clustering.size(), max); i++) {
                double maxP = -1;
                List<List<Node>> maxClustering = null;
                Graph maxStructureGraph = null;
                List<Node> maxLatents = null;
                ICovarianceMatrix maxLatentsCov = null;

                ChoiceGenerator gen = new ChoiceGenerator(clustering.size(), i);
                int[] choice;

                while ((choice = gen.next()) != null) {
                    List<List<Node>> subClusters = new ArrayList<>();
                    List<Node> subLatents = new ArrayList<>();
                    for (int aChoice : choice) {
                        subClusters.add(clustering.get(aChoice));
                        subLatents.add(latents.get(aChoice));
                    }
                    ICovarianceMatrix measuresCovSub = measuresCov.getSubmatrix(getAllVarNames(subClusters));
                    Pair pair = getCov(subClusters, subLatents, measuresCovSub);
                    CovarianceMatrix latentsCov = getLatentsCov(latents, pair.getMatrix(), measuresCovSub.getSampleSize());
                    Graph graph = getGraph(latentsCov);
                    double p = getP(graph, measuresCovSub, subLatents, N, pair.getGls());

                    System.out.println("P = " + p + "  graph = " + graph);

                    if (p > maxP) {
                        maxP = p;
                        maxClustering = new ArrayList<>(subClusters);
                        maxStructureGraph = new EdgeListGraph(graph);
                        maxLatents = new ArrayList<>(subLatents);
                        maxLatentsCov = latentsCov;
                    }
                }

                if (maxClustering != null) {
                    bestMaxP = maxP;
                    bestClustering = maxClustering;
                    bestSructureGraph = maxStructureGraph;
                    bestLatents = maxLatents;
                    bestLatentsCov = maxLatentsCov;
                }
            }

            this.structureGraph = bestSructureGraph;
            this.pValue = bestMaxP;
            this.clustering = bestClustering;
            this.latents = bestLatents;
            this.latentsCov = bestLatentsCov;

            System.out.println("Final P = " + pValue + " graph = " + bestSructureGraph);
            return bestSructureGraph;
        } else {
            List<String> allVarNames = getAllVarNames(clustering);
            ICovarianceMatrix measuresCovSub = measuresCov.getSubmatrix(allVarNames);
            Pair pair = getCov(clustering, latents, measuresCovSub);
            CovarianceMatrix _cov = getLatentsCov(latents, pair.getMatrix(), measuresCovSub.getSampleSize());
            Graph graph = getGraph(_cov);
            this.pValue = getP(graph, measuresCovSub, latents, N, pair.getGls());
            this.clustering = clustering;
            this.latentsCov = _cov;
            this.structureGraph = graph;
            return graph;
        }
    }

    private List<String> getAllVarNames(List<List<Node>> _clustering) {
        List<String> allVarNames = new ArrayList<>();

        for (List<Node> cluster : _clustering) {
            for (Node node : cluster) allVarNames.add(node.getName());
        }

        return allVarNames;
    }

    private double getP(Graph graph, ICovarianceMatrix measuresCov, List<Node> latents, int n, double gls) {
        int m = measuresCov.getDimension();
        double dof = (m * (m + 1)) / 2 - 2 * m; // - graph.getNumEdges();
        double p = 1.0 - new ChiSquaredDistribution(dof).cumulativeProbability((n - 1) * gls);

        int l = latents.size();
        double N = (n - 1) * gls;

        System.out.println("gls = " + gls + " N = " + N + " dof = " + dof + " p = " + p + " m = " + m + " l = " + l);

        return p;
    }

    private Graph getGraph(ICovarianceMatrix _cov) {
        Graph graph;

//        Cpc search = new Cpc(new IndTestFisherZ(latentscov, alpha));
//        search.setKnowledge(knowledge);
//        graph = search.search();
//
//        PcMax search = new PcMax(new IndTestFisherZ(latentscov, .001));
//        search.setKnowledge(knowledge);
//        graph = search.search();

        Fgs search = new Fgs(_cov);
        search.setPenaltyDiscount(2);
        search.setKnowledge(knowledge);
        graph = search.search();
//
//
//        Ges search = new Ges(latentscov);
//        search.setPenaltyDiscount(2);
//        search.setKnowledge(knowledge);
//        graph = search.search();

        graph = new EdgeListGraph(graph);
        GraphUtils.fruchtermanReingoldLayout(graph);
        return graph;
    }

    private CovarianceMatrix getLatentsCov(List<Node> variables, TetradMatrix latentsCov, int sampleSize) {
        return new CovarianceMatrix(variables, latentsCov, sampleSize);
    }

    private static class Pair {
        private TetradMatrix matrix;
        private double gls;

        public Pair(TetradMatrix matrix, double gls) {
            this.matrix = matrix;
            this.gls = gls;
        }

        public TetradMatrix getMatrix() {
            return matrix;
        }

        public double getGls() {
            return gls;
        }
    }

    private Pair getCov(List<List<Node>> clustering, List<Node> latents, ICovarianceMatrix measuresCov) {
        List<List<Node>> _clustering = new ArrayList<>();

        for (List<Node> cluster : clustering) {
            List<Node> _cluster = new ArrayList<Node>();

            for (Node node : cluster) {
                _cluster.add(measuresCov.getVariable(node.getName()));
            }

            _clustering.add(_cluster);
        }

        TetradMatrix measurescov = measuresCov.getMatrix();
        TetradMatrix latentscov = new TetradMatrix(latents.size(), latents.size());

        for (int i = 0; i < latentscov.rows(); i++) {
            for (int j = 0; j < latentscov.columns(); j++) {
                latentscov.set(i, j, (i == j) ? 1.0 : 0.0);
            }
        }

        double[][] loadings = new double[_clustering.size()][];

        for (int i = 0; i < _clustering.size(); i++) {
            loadings[i] = new double[_clustering.get(i).size()];
        }

        for (int i = 0; i < _clustering.size(); i++) {
            loadings[i] = new double[_clustering.get(i).size()];

            for (int j = 0; j < clustering.get(i).size(); j++) {
                loadings[i][j] = 0.0;
            }
        }

        int[][] indicatorIndices = new int[_clustering.size()][];
        List<Node> measures = measuresCov.getVariables();

        for (int i = 0; i < _clustering.size(); i++) {
            indicatorIndices[i] = new int[_clustering.get(i).size()];

            for (int j = 0; j < _clustering.get(i).size(); j++) {
                indicatorIndices[i][j] = measures.indexOf(_clustering.get(i).get(j));
            }
        }

        // Variances of the measures.
        double[] delta = new double[measurescov.columns()];

        for (int i = 0; i < delta.length; i++) {
            delta[i] = 1;
        }

        optimizeQuick(measurescov, latentscov, loadings, indicatorIndices);
        double gls = optimizeAllParameters(measurescov, indicatorIndices, latentscov, loadings, delta);
//        optimizeDeltaParameters(measurescov, indicatorIndices, latentscov, loadings, delta);

        return new Pair(latentscov, gls);
    }

    public List<List<Node>> getClustering() {
        return clustering;
    }

    public double getAlpha() {
        return alpha;
    }

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }

    public IKnowledge getKnowledge() {
        return knowledge;
    }

    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public ICovarianceMatrix getLatentsCov() {
        return this.latentsCov;
    }

    public double getPValue() {
        return pValue;
    }

    /**
     * @return the full discovered graph, with latents and indicators.
     */
    public Graph getFullGraph() {
        Graph graph = new EdgeListGraph(structureGraph);

        for (int i = 0; i < this.latents.size(); i++) {
            Node latent = this.latents.get(i);
            List<Node> measuredGuys = getClustering().get(i);

            for (Node measured : measuredGuys) {
                if (!graph.containsNode(measured)) {
                    graph.addNode(measured);
                }

                graph.addDirectedEdge(latent, measured);
            }
        }

        return graph;
    }

    //=================================== PRIVATE METHODS =========================================//

    private List<Node> defineLatents(List<String> names) {
        List<Node> latents = new ArrayList<>();

        for (String name : names) {
            Node node = new GraphNode(name);
            node.setNodeType(NodeType.LATENT);
            latents.add(node);
        }

        return latents;
    }

    private void removeSmallClusters(List<Node> latents, List<List<Node>> clustering, int minimumSize) {
        List<List<Node>> _clustering = new ArrayList<>(clustering);
        List<Node> _latents = new ArrayList<>(latents);

        for (int i = 0; i >= _clustering.size(); i++) {
            if (_clustering.get(i).size() < minimumSize) {
                clustering.remove(_clustering.get(i));
                latents.remove(_latents.get(i));
            }
        }
    }

    private void optimizeQuick(TetradMatrix measurescov, TetradMatrix latentscov,
                               double[][] loadings, int[][] indicatorIndices) {
        Function1 function = new Function1(indicatorIndices, measurescov, loadings, latentscov);
        MultivariateOptimizer search = new PowellOptimizer(threshold, threshold);

        search.optimize(
                new InitialGuess(getAllValues(latentscov, loadings)),
                new ObjectiveFunction(function),
                GoalType.MINIMIZE,
                new MaxEval(100000));
    }

    /**
     * @return the minimum GLS value.
     */
    private double optimizeDeltaParameters(TetradMatrix measurescov,
                                           int[][] indicatorIndices, TetradMatrix latentscov, double[][] loadings,
                                           double[] delta) {
        Function4 function = new Function4(indicatorIndices, measurescov, loadings, latentscov);
        MultivariateOptimizer search = new PowellOptimizer(threshold, threshold);

        PointValuePair pair = search.optimize(
                new InitialGuess(getAllValues(delta)),
                new ObjectiveFunction(function),
                GoalType.MINIMIZE,
                new MaxEval(1000000));

        return pair.getValue();
    }

    /**
     * @return the minimum GLS value.
     */
    private double optimizeAllParameters(TetradMatrix measurescov,
                                         int[][] indicatorIndices, TetradMatrix latentscov, double[][] loadings,
                                         double[] delta) {
        Function5 function = new Function5(indicatorIndices, measurescov, loadings, latentscov, delta);
        MultivariateOptimizer search = new PowellOptimizer(threshold, threshold);

        PointValuePair pair = search.optimize(
                new InitialGuess(getAllValues(latentscov, loadings, getAllValues(latentscov, loadings, delta))),
                new ObjectiveFunction(function),
                GoalType.MINIMIZE,
                new MaxEval(1000000));

        return pair.getValue();
    }


    private double[] getAllValues(TetradMatrix latentscov, double[][] loadings, double[] delta) {
        List<Double> _values = new ArrayList<>();

        for (int i = 0; i < loadings.length; i++) {
            _values.add(latentscov.get(i, i));
        }

        for (int i = 0; i < loadings.length; i++) {
            for (int j = i + 1; j < loadings.length; j++) {
                _values.add(latentscov.get(i, j));
            }
        }

        for (double[] loading : loadings) {
            for (double aLoading : loading) {
                _values.add(aLoading);
            }
        }

        for (double aDelta : delta) {
            _values.add(aDelta);
        }

        double[] values = new double[_values.size()];

        for (int i = 0; i < _values.size(); i++) {
            values[i] = _values.get(i);
        }

        return values;
    }

    private double[] getAllValues(TetradMatrix latentscov, double[][] loadings) {
        List<Double> _values = new ArrayList<>();

        for (int i = 0; i < loadings.length; i++) {
            _values.add(latentscov.get(i, i));
        }

        for (int i = 0; i < loadings.length; i++) {
            for (int j = i + 1; j < loadings.length; j++) {
                _values.add(latentscov.get(i, j));
            }
        }

        for (double[] loading : loadings) {
            for (double aLoading : loading) {
                _values.add(aLoading);
            }
        }
        double[] values = new double[_values.size()];

        for (int i = 0; i < _values.size(); i++) {
            values[i] = _values.get(i);
        }

        return values;
    }

    private double[] getAllValues(double[] delta) {
        List<Double> _values = new ArrayList<>();

        for (double aDelta : delta) {
            _values.add(aDelta);
        }

        double[] values = new double[_values.size()];

        for (int i = 0; i < _values.size(); i++) {
            values[i] = _values.get(i);
        }

        return values;
    }

    /**
     * Clusters smaller than this size will be tossed out.
     */
    public int getMinClusterSize() {
        return minClusterSize;
    }

    public void setMinClusterSize(int minClusterSize) {
        if (minClusterSize < 3)
            throw new IllegalArgumentException("Minimum cluster size must be >= 3: " + minClusterSize);
        this.minClusterSize = minClusterSize;
    }

    public void setSearchForMax(boolean searchForMax) {
        this.searchForMax = searchForMax;
    }

    private class Function1 implements org.apache.commons.math3.analysis.MultivariateFunction {
        private final int[][] indicatorIndices;
        private final TetradMatrix measurescov;
        private final double[][] loadings;
        private final TetradMatrix latentscov;

        public Function1(int[][] indicatorIndices, TetradMatrix measurescov, double[][] loadings,
                         TetradMatrix latentscov) {
            this.indicatorIndices = indicatorIndices;
            this.measurescov = measurescov;
            this.loadings = loadings;
            this.latentscov = latentscov;
        }

        @Override
        public double value(double[] values) {
            int count = 0;

            for (int i = 0; i < loadings.length; i++) {
                latentscov.set(i, i, values[count]);
                count++;
            }

            for (int i = 0; i < loadings.length; i++) {
                for (int j = i + 1; j < loadings.length; j++) {
                    latentscov.set(i, j, values[count]);
                    latentscov.set(j, i, values[count]);
                    count++;
                }
            }

            for (int i = 0; i < loadings.length; i++) {
                for (int j = 0; j < loadings[i].length; j++) {
                    loadings[i][j] = values[count];
                    count++;
                }
            }

            for (int i = 0; i < loadings.length; i++) {
                if (values[i] <= 0) return Double.POSITIVE_INFINITY;
            }

            return sumOfDifferences(indicatorIndices, measurescov, loadings, latentscov);
        }
    }

    private class Function4 implements org.apache.commons.math3.analysis.MultivariateFunction {
        private final int[][] indicatorIndices;
        private final TetradMatrix measurescov;
        private TetradMatrix measuresCovInverse;
        private final double[][] loadings;
        private final TetradMatrix latentscov;

        public Function4(int[][] indicatorIndices, TetradMatrix measurescov, double[][] loadings, TetradMatrix latentscov) {
            this.indicatorIndices = indicatorIndices;
            this.measurescov = measurescov;
            this.loadings = loadings;
            this.latentscov = latentscov;
            this.measuresCovInverse = measurescov.inverse();
        }

        @Override
        public double value(double[] delta) {

            // The delta here are just the measured variances; other delta are being held fixed.
            TetradMatrix implied = impliedCovariance(indicatorIndices, loadings, measurescov, latentscov, delta);

            TetradMatrix I = TetradMatrix.identity(implied.rows());
            TetradMatrix diff = I.minus((implied.times(measuresCovInverse)));  // time hog. times().

            return 0.5 * (diff.times(diff)).trace();
        }
    }

    private class Function5 implements org.apache.commons.math3.analysis.MultivariateFunction {
        private final int[][] indicatorIndices;
        private final TetradMatrix measurescov;
        private TetradMatrix measuresCovInverse;
        private final double[][] loadings;
        private final TetradMatrix latentscov;
        private final double[] delta;

        public Function5(int[][] indicatorIndices, TetradMatrix measurescov, double[][] loadings, TetradMatrix latentscov,
                         double[] delta) {
            this.indicatorIndices = indicatorIndices;
            this.measurescov = measurescov;
            this.loadings = loadings;
            this.latentscov = latentscov;
            this.delta = delta;
            this.measuresCovInverse = measurescov.inverse();
        }

        @Override
        public double value(double[] values) {
            int count = 0;

            for (int i = 0; i < loadings.length; i++) {
                latentscov.set(i, i, values[count]);
                count++;
            }

            for (int i = 0; i < loadings.length; i++) {
                for (int j = i + 1; j < loadings.length; j++) {
                    latentscov.set(i, j, values[count]);
                    latentscov.set(j, i, values[count]);
                    count++;
                }
            }

            for (int i = 0; i < loadings.length; i++) {
                for (int j = 0; j < loadings[i].length; j++) {
                    loadings[i][j] = values[count];
                    count++;
                }
            }

            for (int i = 0; i < delta.length; i++) {
                delta[i] = values[count];
                count++;
            }

            TetradMatrix implied = impliedCovariance(indicatorIndices, loadings, measurescov, latentscov, delta);

            TetradMatrix I = TetradMatrix.identity(implied.rows());
            TetradMatrix diff = I.minus((implied.times(measuresCovInverse)));  // time hog. times().

            return 0.5 * (diff.times(diff)).trace();
        }
    }

    private TetradMatrix impliedCovariance(int[][] indicatorIndices, double[][] loadings, TetradMatrix cov, TetradMatrix loadingscov,
                                           double[] delta) {
        TetradMatrix implied = new TetradMatrix(cov.rows(), cov.columns());

        for (int i = 0; i < loadings.length; i++) {
            for (int k = 0; k < loadings[i].length; k++) {
                for (int l = k + 1; l < loadings[i].length; l++) {
                    double prod = loadings[i][k] * loadings[i][l] * loadingscov.get(i, i);
                    implied.set(indicatorIndices[i][k], indicatorIndices[i][l], prod);
                    implied.set(indicatorIndices[i][l], indicatorIndices[i][k], prod);
                }
            }
        }

        for (int i = 0; i < loadings.length; i++) {
            for (int j = i + 1; j < loadings.length; j++) {
                for (int k = 0; k < loadings[i].length; k++) {
                    for (int l = 0; l < loadings[j].length; l++) {
                        double prod = loadings[i][k] * loadings[j][l] * loadingscov.get(i, j);
                        implied.set(indicatorIndices[i][k], indicatorIndices[j][l], prod);
                        implied.set(indicatorIndices[j][l], indicatorIndices[i][k], prod);
                    }
                }
            }
        }

        for (int i = 0; i < implied.rows(); i++) {
            implied.set(i, i, delta[i]);
        }

        return implied;
    }

    private double sumOfDifferences(int[][] indicatorIndices, TetradMatrix cov, double[][] loadings, TetradMatrix loadingscov) {
        double sum = 0;

        for (int i = 0; i < loadings.length; i++) {
            for (int k = 0; k < loadings[i].length; k++) {
                for (int l = k + 1; l < loadings[i].length; l++) {
                    double _cov = cov.get(indicatorIndices[i][k], indicatorIndices[i][l]);
                    double prod = loadings[i][k] * loadings[i][l] * loadingscov.get(i, i);
                    double diff = _cov - prod;
                    sum += diff * diff;
                }
            }
        }

        for (int i = 0; i < loadings.length; i++) {
            for (int j = i + 1; j < loadings.length; j++) {
                for (int k = 0; k < loadings[i].length; k++) {
                    for (int l = 0; l < loadings[j].length; l++) {
                        double _cov = cov.get(indicatorIndices[i][k], indicatorIndices[j][l]);
                        double prod = loadings[i][k] * loadings[j][l] * loadingscov.get(i, j);
                        double diff = _cov - prod;
                        sum += diff * diff;
                    }
                }
            }
        }

        return sum;
    }

}



