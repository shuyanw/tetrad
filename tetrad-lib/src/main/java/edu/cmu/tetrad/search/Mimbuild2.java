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
import edu.cmu.tetrad.util.MatrixUtils;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.distribution.ChiSquaredDistribution;
import org.apache.commons.math3.optim.InitialGuess;
import org.apache.commons.math3.optim.MaxEval;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.GoalType;
import org.apache.commons.math3.optim.nonlinear.scalar.MultivariateOptimizer;
import org.apache.commons.math3.optim.nonlinear.scalar.ObjectiveFunction;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.PowellOptimizer;

import java.util.*;

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
    private List<List<Node>> clustering;
    double threshold = 1e-7;

    public Method getMethod() {
        return method;
    }

    public void setMethod(Method method) {
        this.method = method;
    }

    public enum Method {MAX, MAX_LARGE, STRAIGHT, CLUSTERS_FIRST}

    private Method method = Method.STRAIGHT;

    public Mimbuild2() {
    }

    //=================================== PUBLIC METHODS =========================================//

    public Graph search(List<List<Node>> _clustering, List<String> latentNames, ICovarianceMatrix measuresCov) {
        List<List<Node>> clustering = convertClusters(_clustering, measuresCov);
        measuresCov = measuresCov.getSubmatrix(getAllVarNames(clustering));
        int N = measuresCov.getSampleSize();
        List<String> _latentNames = new ArrayList<>(latentNames);
        clustering = new ArrayList<>(clustering);

        List<Node> latents = defineLatents(_latentNames);
        this.latents = latents;

        // This removes the small clusters and their names.
        removeSmallClusters(latents, clustering, getMinClusterSize());

        if (method == Method.MAX_LARGE) {
            double bestMaxP = -1;
            List<List<Node>> bestClustering = null;
            Graph bestSructureGraph = null;
            List<Node> bestLatents = null;
            ICovarianceMatrix bestLatentsCov = null;

            int min = 2;
            int max = clustering.size();

            for (int i = Math.min(clustering.size(), max); i >= min; i--) {
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
                    CovRet pair = getCov(subClusters, subLatents, measuresCov);
                    CovarianceMatrix latentsCov = getLatentsCov(subLatents, pair.getLatentscov(),
                            measuresCovSub.getSampleSize());
                    Graph graph = getGraph(latentsCov);
                    double p = getP(measuresCovSub, subLatents, N, pair.getGls());


                    if (p > alpha) {
                        System.out.println("P = " + p + "  graph = " + graph);
                    }

                    if (p > maxP) {
                        maxP = p;
                        maxClustering = new ArrayList<>(subClusters);
                        maxStructureGraph = new EdgeListGraph(graph);
                        maxLatents = new ArrayList<>(subLatents);
                        maxLatentsCov = latentsCov;
                    }
                }

                if (maxP > alpha && maxClustering != null) {
                    bestMaxP = maxP;
                    bestClustering = maxClustering;
                    bestSructureGraph = maxStructureGraph;
                    bestLatents = maxLatents;
                    bestLatentsCov = maxLatentsCov;
                    break;
                }
            }

            this.structureGraph = bestSructureGraph;
            this.pValue = bestMaxP;
            this.clustering = bestClustering;
            this.latents = bestLatents;
            this.latentsCov = bestLatentsCov;

            if (structureGraph != null) {
                System.out.println("Final P = " + pValue + " graph = " + bestSructureGraph);
                printDagP(SearchGraphUtils.dagFromPattern(this.structureGraph), this.clustering, this.latents,
                        measuresCov.getSubmatrix(getAllVarNames(this.clustering)));
            }
            return bestSructureGraph;
        } else if (method == Method.MAX) {
            int min = 1;
            int max = clustering.size();
            List<List<Node>> maxClustering = null;
            Graph maxStructureGraph = null;
            List<Node> maxLatents = null;
            ICovarianceMatrix maxLatentsCov = null;
            double maxP = -1;

            for (int i = min; i <= Math.min(clustering.size(), max); i++) {

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
                    CovRet pair = getCov(subClusters, subLatents, measuresCovSub);
                    CovarianceMatrix latentsCov = getLatentsCov(subLatents, pair.getLatentscov(), measuresCovSub.getSampleSize());
                    Graph graph = getGraph(latentsCov);
                    double p = getP(measuresCovSub, subLatents, N, pair.getGls());

                    System.out.println("P = " + p + "  graph = " + graph);

                    if (p > maxP) {
                        maxP = p;
                        maxClustering = new ArrayList<>(subClusters);
                        maxStructureGraph = new EdgeListGraph(graph);
                        maxLatents = new ArrayList<>(subLatents);
                        maxLatentsCov = latentsCov;
                        System.out.println("NEW MAX P = " + p + "  graph = " + graph);
                    }
                }
            }

            this.structureGraph = maxStructureGraph;
            this.pValue = maxP;
            this.clustering = maxClustering;
            this.latents = maxLatents;
            this.latentsCov = maxLatentsCov;

            System.out.println("Final P = " + pValue + " graph = " + maxStructureGraph);
            printDagP(SearchGraphUtils.dagFromPattern(this.structureGraph), this.clustering, this.latents,
                    measuresCov.getSubmatrix(getAllVarNames(this.clustering)));
            return maxStructureGraph;
        } else if (method == Method.STRAIGHT) {
            List<String> allVarNames = getAllVarNames(clustering);
            ICovarianceMatrix measuresCovSub = measuresCov.getSubmatrix(allVarNames);
            CovRet pair = getCov(clustering, latents, measuresCovSub);
            CovarianceMatrix _cov = getLatentsCov(latents, pair.getLatentscov(), measuresCovSub.getSampleSize());
            Graph graph = getGraph(_cov);
            this.pValue = getP(measuresCovSub, latents, N, pair.getGls());
            this.clustering = clustering;
            this.latentsCov = _cov;
            this.structureGraph = graph;
            printDagP(SearchGraphUtils.dagFromPattern(this.structureGraph), this.clustering, this.latents,
                    measuresCov.getSubmatrix(getAllVarNames(this.clustering)));
            return graph;
        } else if (method == Method.CLUSTERS_FIRST) {
            List<String> allVarNames = getAllVarNames(clustering);
            ICovarianceMatrix measuresCovSub = measuresCov.getSubmatrix(allVarNames);
            CovRet pair = getCov2(clustering, latents, measuresCov);
            CovarianceMatrix _cov = getLatentsCov(latents, pair.getLatentscov(), measuresCov.getSampleSize());
            Graph graph = getGraph(_cov);
            this.pValue = getP(measuresCovSub, latents, N, pair.getGls());
            this.clustering = clustering;
            this.latentsCov = _cov;
            this.structureGraph = graph;

            System.out.println("CLUSTER FIRST P = " + this.pValue + "  graph = " + graph);

            printDagP(SearchGraphUtils.dagFromPattern(this.structureGraph), this.clustering, this.latents,
                    measuresCov.getSubmatrix(getAllVarNames(this.clustering)));
            return graph;
        } else {
            throw new IllegalStateException();
        }
    }

    private void printDagP(Graph dag, List<List<Node>> clustering, List<Node> latents,
                           ICovarianceMatrix measuresCov) {
        CovRet ret = getDagCov2(clustering, latents, dag, measuresCov);
        int n = measuresCov.getSampleSize();
        double p = getP(measuresCov, latents, n, ret.getGls());
        System.out.println("P for DAG = " + p + " DAG = " + dag);
    }

    private List<String> getAllVarNames(List<List<Node>> _clustering) {
        List<String> allVarNames = new ArrayList<>();

        for (List<Node> cluster : _clustering) {
            for (Node node : cluster) allVarNames.add(node.getName());
        }

        return allVarNames;
    }

    private List<Node> getAllVars(List<List<Node>> _clustering) {
        List<Node> allVars = new ArrayList<>();

        for (List<Node> cluster : _clustering) {
            allVars.addAll(cluster);
        }

        return allVars;
    }


    private double getP(ICovarianceMatrix measuresCov, List<Node> latents, int n, double gls) {
        int m = measuresCov.getDimension();
        int l = latents.size();
        double dof = (m * (m + 1)) / 2 - 2 * m - (l * (l + 1) / 2);
        double p = 1.0 - new ChiSquaredDistribution(dof).cumulativeProbability((n - 1) * gls);
        double N = (n - 1) * gls;
        System.out.println("gls = " + gls + " N = " + N + " dof = " + dof + " p = " + p + " m = " + m + " l = " + l);
        return p;
    }

    private Graph getGraph(ICovarianceMatrix _cov) {
        Graph graph;

//        Cpc search = new Cpc(new IndTestFisherZ(_cov, alpha));
//        search.setKnowledge(knowledge);
//        graph = search.search();

//        PcMax search = new PcMax(new IndTestFisherZ(_cov, .001));
//        search.setKnowledge(knowledge);
//        graph = search.search();

        Fgs search = new Fgs(_cov);
        search.setPenaltyDiscount(2);
        search.setKnowledge(knowledge);
        graph = search.search();
//
//
//        Ges search = new Ges(_cov);
//        search.setPenaltyDiscount(2);
//        search.setKnowledge(knowledge);
//        graph = search.search();

        graph = new EdgeListGraph(graph);
        GraphUtils.fruchtermanReingoldLayout(graph);
        return graph;
    }

    private CovarianceMatrix getLatentsCov(List<Node> latentVars, TetradMatrix latentsCov, int sampleSize) {
        return new CovarianceMatrix(latentVars, latentsCov, sampleSize);
    }

    private CovRet getCov(List<List<Node>> _clustering, List<Node> latents, ICovarianceMatrix _measuresCov) {
        ICovarianceMatrix measuresCovSub = _measuresCov.getSubmatrix(getAllVarNames(_clustering));
        TetradMatrix measurescov = measuresCovSub.getMatrix();
        TetradMatrix latentscov = getLatentsCovInit(latents);

        double[][] loadings = getLoadingInit(_clustering);
        int[][] indicatorIndices = getIndicatorIndicesInit(_clustering);

        // Variances of the measures.
        double[] delta = initDelta(getAllVarNames(_clustering).size());

        optimizeQuick(measurescov, latentscov, loadings, indicatorIndices);
        double gls = optimizeAllParameters(measurescov, latentscov, indicatorIndices, loadings, delta);
//        optimizeDeltaParameters(measurescov, indicatorIndices, latentscov, loadings, deltas);

        return new CovRet(latentscov, loadings, indicatorIndices, delta, gls, getAllVars(_clustering));
    }

    private CovRet getDagCov(List<List<Node>> _clustering, List<Node> latents, Graph dag,
                             ICovarianceMatrix _measuresCov) {
        ICovarianceMatrix measuresCovSub = _measuresCov.getSubmatrix(getAllVarNames(_clustering));
        TetradMatrix measurescov = measuresCovSub.getMatrix();
        boolean[][] mask = getLatentEdgeCoefMask(latents, dag);
        TetradMatrix latentCoef = getLatentEdgeCoefInit(latents);
        boolean[][] maskErr = getLatentErrCovarMask(latents);
        TetradMatrix latentErrCovar = getLatentErrCovarInit(latents);
        double[][] loadings = getLoadingInit(_clustering);
        int[][] indicatorIndices = getIndicatorIndicesInit(_clustering);

        // Variances of the measures.
        double[] delta = initDelta(getAllVarNames(_clustering).size());

        optimizeQuickDag(measurescov, mask, latentCoef, maskErr, latentErrCovar, loadings, indicatorIndices);
        double gls = optimizeAllParametersDag(measurescov, mask, latentCoef,
                maskErr, latentErrCovar, indicatorIndices, loadings, delta);
//        optimizeDeltaParameters(measurescov, indicatorIndices, latentscov, loadings, deltas);

        TetradMatrix latentscov = MatrixUtils.impliedCovar(latentCoef, latentErrCovar);

        return new CovRet(latentscov, loadings, indicatorIndices, delta, gls, getAllVars(_clustering));
    }

    private CovRet getCov2(List<List<Node>> _clustering, List<Node> latents, ICovarianceMatrix _measuresCov) {
        ICovarianceMatrix measuresCov = _measuresCov.getSubmatrix(getAllVarNames(_clustering));
        double[][] loadings = getLoadingInit(_clustering);
        int[][] indicatorIndices = getIndicatorIndicesInit(_clustering);
        double[][] deltas = getDeltasInit(_clustering);
        List<Node> allVars = getAllVars(_clustering);

        for (int i = 0; i < latents.size(); i++) {
            CovRet ret = getSingleFactorCov(i, _clustering, latents, measuresCov);
            List<Node> theseVars = ret.getVariables();
            loadings[i] = ret.getLoadings()[0];
            indicatorIndices[i] = translate(ret.getIndicatorIndices()[0], theseVars, allVars);
            deltas[i] = ret.getDeltas();
        }

        // Variances of the measures.
        double[] allDelta = combineDeltas(deltas);
        TetradMatrix measurescov = measuresCov.getMatrix();
        TetradMatrix latentscov = getLatentsCovInit(latents);

        optimizeQuick(measurescov, latentscov, loadings, indicatorIndices);
        double gls = optimizeAllParameters(measurescov, latentscov, indicatorIndices, loadings, allDelta);
//        optimizeQuickLatent(measurescov, latentscov, loadings, indicatorIndices);
//        double gls = optimizeLatentDeltaParams(measurescov, indicatorIndices, latentscov, loadings, allDelta);
//        optimizeDeltaParameters(measurescov, indicatorIndices, latentscov, loadings, deltas);

        return new CovRet(latentscov, loadings, indicatorIndices, allDelta, gls, getAllVars(_clustering));
    }

    private CovRet getDagCov2(List<List<Node>> _clustering, List<Node> latents, Graph dag, ICovarianceMatrix measuresCov) {
        ICovarianceMatrix _measurescov = measuresCov.getSubmatrix(getAllVarNames(_clustering));
        TetradMatrix measurescov = _measurescov.getMatrix();
        boolean[][] mask = getLatentEdgeCoefMask(latents, dag);
        TetradMatrix latentCoef = getLatentEdgeCoefInit(latents);
        boolean[][] maskErr = getLatentErrCovarMask(latents);
        TetradMatrix latentErrCovar = getLatentErrCovarInit(latents);
        double[][] loadings = getLoadingInit(_clustering);
        int[][] indicatorIndices = getIndicatorIndicesInit(_clustering);
        double[][] deltas = getDeltasInit(_clustering);
        List<Node> allVars = getAllVars(_clustering);

        for (int i = 0; i < latents.size(); i++) {
            CovRet ret = getSingleFactorCov(i, _clustering, latents, measuresCov);
            List<Node> theseVars = ret.getVariables();
            loadings[i] = ret.getLoadings()[0];
            indicatorIndices[i] = translate(ret.getIndicatorIndices()[0], theseVars, allVars);
            deltas[i] = ret.getDeltas();
        }

        // Variances of the measures.
        double[] allDelta = combineDeltas(deltas);

        optimizeQuickDag(measurescov, mask, latentCoef, maskErr, latentErrCovar, loadings, indicatorIndices);
        double gls = optimizeAllParametersDag(measurescov, mask, latentCoef, maskErr,
                latentErrCovar, indicatorIndices, loadings, allDelta);
//        optimizeDeltaParameters(measurescov, indicatorIndices, latentscov, loa22dings, deltas);

        TetradMatrix latentscov = MatrixUtils.impliedCovar(latentCoef, latentErrCovar);

        return new CovRet(latentscov, loadings, indicatorIndices, allDelta, gls, getAllVars(_clustering));
    }

    private int[][] getIndicatorIndicesInit(List<List<Node>> _clustering) {
        int[][] indicatorIndices = new int[_clustering.size()][];
        List<Node> measures = getAllVars(_clustering);

        for (int i = 0; i < _clustering.size(); i++) {
            indicatorIndices[i] = new int[_clustering.get(i).size()];

            for (int j = 0; j < _clustering.get(i).size(); j++) {
                indicatorIndices[i][j] = measures.indexOf(_clustering.get(i).get(j));
            }
        }

        return indicatorIndices;
    }

    private int[] translate(int[] indices, List<Node> theseVars, List<Node> allVars) {
        int[] translated = new int[indices.length];

        for (int i = 0; i < indices.length; i++) {
            int ind = allVars.indexOf(theseVars.get(i));
            translated[i] = ind;
        }

        return translated;
    }

    private double[] combineDeltas(double[][] deltas) {
        List<Double> allDeltas = new ArrayList<>();

        for (int i = 0; i < deltas.length; i++) {
            for (int j = 0; j < deltas[i].length; j++) {
                allDeltas.add(deltas[i][j]);
            }
        }

        double[] _allDeltas = new double[allDeltas.size()];
        for (int i = 0; i < allDeltas.size(); i++) _allDeltas[i] = allDeltas.get(i);
        return _allDeltas;
    }

    private Map<List<Node>, CovRet> singleFactorResults = new HashMap<>();

    private CovRet getSingleFactorCov(int i, List<List<Node>> clustering, List<Node> latents, ICovarianceMatrix measuresCov) {
        List<Node> key = clustering.get(i);

        if (singleFactorResults.get(key) == null) {
            List<List<Node>> clusteringSub = Collections.singletonList(key);
            List<Node> latentsSub = Collections.singletonList(latents.get(i));
            CovRet covRet = getCov(clusteringSub, latentsSub, measuresCov.getSubmatrix(getAllVarNames(clustering)));
            singleFactorResults.put(key, covRet);
        }

        return singleFactorResults.get(key);
    }

    private List<List<Node>> convertClusters(List<List<Node>> clustering, ICovarianceMatrix measuresCov) {
        List<List<Node>> _clustering = new ArrayList<>();

        for (List<Node> cluster : clustering) {
            List<Node> _cluster = new ArrayList<Node>();

            for (Node node : cluster) {
                _cluster.add(measuresCov.getVariable(node.getName()));
            }

            _clustering.add(_cluster);
        }
        return _clustering;
    }

    private double[][] getDeltasInit(List<List<Node>> _clustering) {
        double[][] deltas = new double[_clustering.size()][];

        for (int i = 0; i < _clustering.size(); i++) {
            deltas[i] = new double[_clustering.get(i).size()];

            for (int j = 0; j < _clustering.get(i).size(); j++) {
                deltas[i][j] = 1;
            }
        }
        return deltas;
    }

    private double[][] getLoadingInit(List<List<Node>> clustering) {
        double[][] loadings = new double[clustering.size()][];

        for (int i = 0; i < clustering.size(); i++) {
            loadings[i] = new double[clustering.get(i).size()];
        }

        for (int i = 0; i < clustering.size(); i++) {
            loadings[i] = new double[clustering.get(i).size()];

            for (int j = 0; j < clustering.get(i).size(); j++) {
                loadings[i][j] = 0.0;
            }
        }

        return loadings;
    }

    private double[] initDelta(int numMeasures) {
        double[] delta = new double[numMeasures];

        for (int i = 0; i < delta.length; i++) {
            delta[i] = 1;
        }

        return delta;
    }

    private TetradMatrix getLatentsCovInit(List<Node> latents) {
        TetradMatrix latentscov = new TetradMatrix(latents.size(), latents.size());

        for (int i = 0; i < latentscov.rows(); i++) {
            latentscov.set(i, i, 1.0);
        }

        return latentscov;
    }

    private boolean[][] getLatentEdgeCoefMask(List<Node> latents, Graph dag) {
        boolean[][] mask = new boolean[latents.size()][latents.size()];

        for (int i = 0; i < latents.size(); i++) {
            for (int j = 0; j < latents.size(); j++) {
                Edge edge = dag.getEdge(latents.get(i), latents.get(j));
                if (edge != null && edge.pointsTowards(latents.get(j))) {
                    mask[j][i] = true;
                } else {
                    mask[j][i] = false;
                }
            }
        }

        return mask;
    }

    private boolean[][] getLatentErrCovarMask(List<Node> latents) {
        boolean[][] mask = new boolean[latents.size()][latents.size()];

        for (int i = 0; i < latents.size(); i++) {
            mask[i][i] = true;
        }

        return mask;
    }

    private TetradMatrix getLatentEdgeCoefInit(List<Node> latents) {
        return new TetradMatrix(latents.size(), latents.size());
    }

    private TetradMatrix getLatentErrCovarInit(List<Node> latents) {
        TetradMatrix errCovar = new TetradMatrix(latents.size(), latents.size());

        for (int i = 0; i < latents.size(); i++) {
            errCovar.set(i, i, 1.0);
        }

        return errCovar;
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
        List<List<Node>> _clustering = new ArrayList<>();
        List<Node> _latents = new ArrayList<>();

        for (int i = 0; i < clustering.size(); i++) {
            if (clustering.get(i).size() >= minimumSize) {
                _clustering.add(clustering.get(i));
                _latents.add(latents.get(i));
            }
        }

        clustering.clear();
        latents.clear();
        clustering.addAll(_clustering);
        latents.addAll(_latents);

        System.out.println();
    }

    private void optimizeQuick(TetradMatrix measurescov, TetradMatrix latentscov,
                               double[][] loadings, int[][] indicatorIndices) {
        Function function = new Function(measurescov, indicatorIndices, loadings, latentscov, null, true);
        MultivariateOptimizer search = new PowellOptimizer(threshold, threshold);

        search.optimize(
                new InitialGuess(getValues(latentscov, loadings, null)),
                new ObjectiveFunction(function),
                GoalType.MINIMIZE,
                new MaxEval(100000));
    }

    private void optimizeQuickDag(TetradMatrix measurescov,
                                  boolean[][] mask, TetradMatrix latentCoef,
                                  boolean[][] maskErr, TetradMatrix latentErrCovar,
                                  double[][] loadings, int[][] indicatorIndices) {
        FunctionDag function = new FunctionDag(measurescov, indicatorIndices, loadings,
                mask, latentCoef, maskErr, latentErrCovar,
                null, true);
        MultivariateOptimizer search = new PowellOptimizer(threshold, threshold);

        search.optimize(
                new InitialGuess(getValuesDag(mask, latentCoef, maskErr, latentErrCovar, loadings, null)),
                new ObjectiveFunction(function),
                GoalType.MINIMIZE,
                new MaxEval(100000));
    }


    /**
     * @return the minimum GLS value.
     */
    private double optimizeAllParameters(TetradMatrix measurescov,
                                         TetradMatrix latentscov, int[][] indicatorIndices, double[][] loadings,
                                         double[] delta) {
        Function function = new Function(measurescov, indicatorIndices, loadings, latentscov, delta, false);
        MultivariateOptimizer search = new PowellOptimizer(threshold, threshold);

        PointValuePair pair = search.optimize(
                new InitialGuess(getValues(latentscov, loadings, delta)),
                new ObjectiveFunction(function),
                GoalType.MINIMIZE,
                new MaxEval(1000000));

        return pair.getValue();
    }

    /**
     * @return the minimum GLS value.
     */
    private double optimizeAllParametersDag(TetradMatrix measurescov,
                                            boolean[][] mask, TetradMatrix latentCoef,
                                            boolean[][] maskErr, TetradMatrix latentErrCovar,
                                            int[][] indicatorIndices, double[][] loadings,
                                            double[] delta) {
        FunctionDag function = new FunctionDag(measurescov, indicatorIndices, loadings,
                mask, latentCoef, maskErr, latentErrCovar,
                delta, false);
        MultivariateOptimizer search = new PowellOptimizer(threshold, threshold);

        PointValuePair pair = search.optimize(
                new InitialGuess(getValuesDag(mask, latentCoef, maskErr, latentErrCovar, loadings, delta)),
                new ObjectiveFunction(function),
                GoalType.MINIMIZE,
                new MaxEval(1000000));

        return pair.getValue();
    }

    private double[] getValues(TetradMatrix latentscov, double[][] loadings, double[] delta) {
        List<Double> _values = new ArrayList<>();

        if (latentscov != null) {
            for (int i = 0; i < latentscov.columns(); i++) {
                _values.add(latentscov.get(i, i));
            }

            for (int i = 0; i < latentscov.columns(); i++) {
                for (int j = i + 1; j < latentscov.columns(); j++) {
                    _values.add(latentscov.get(i, j));
                }
            }
        }

        if (loadings != null) {
            for (double[] loading : loadings) {
                for (double aLoading : loading) {
                    _values.add(aLoading);
                }
            }
        }

        if (delta != null) {
            for (double aDelta : delta) {
                _values.add(aDelta);
            }
        }

        double[] values = new double[_values.size()];

        for (int i = 0; i < _values.size(); i++) {
            values[i] = _values.get(i);
        }

        return values;
    }

    private double[] getValuesDag(boolean[][] mask, TetradMatrix latentCoef,
                                  boolean[][] maskErr, TetradMatrix latentErrCovar,
                                  double[][] loadings, double[] delta) {
        List<Double> _values = new ArrayList<>();


        if (maskErr != null) {
            for (int i = 0; i < maskErr.length; i++) {
                if (maskErr[i][i]) {
                    _values.add(latentErrCovar.get(i, i));
                }
            }

            for (int i = 0; i < mask.length; i++) {
                for (int j = 0; j < mask.length; j++) {
                    if (mask[i][j]) {
                        _values.add(latentCoef.get(i, j));
                    }
                }
            }
        }

        if (loadings != null) {
            for (double[] loading : loadings) {
                for (double aLoading : loading) {
                    _values.add(aLoading);
                }
            }
        }

        if (delta != null) {
            for (double aDelta : delta) {
                _values.add(aDelta);
            }
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

    private class Function implements org.apache.commons.math3.analysis.MultivariateFunction {
        private final int[][] indicatorIndices;
        private final TetradMatrix measurescov;
        private TetradMatrix measuresCovInverse;
        private final double[][] loadings;
        private final TetradMatrix latentscov;
        private final double[] delta;
        private final boolean quick;

        public Function(TetradMatrix measurescov, int[][] indicatorIndices, double[][] loadings, TetradMatrix latentscov,
                        double[] delta, boolean quick) {
            this.indicatorIndices = indicatorIndices;
            this.measurescov = measurescov;
            this.loadings = loadings;
            this.latentscov = latentscov;
            this.delta = delta;
            this.measuresCovInverse = measurescov.inverse();
            this.quick = quick;
        }

        @Override
        public double value(double[] values) {
            int count = 0;

            if (latentscov != null) {
                for (int i = 0; i < latentscov.columns(); i++) {
                    latentscov.set(i, i, values[count]);
                    count++;
                }

                for (int i = 0; i < latentscov.columns(); i++) {
                    for (int j = i + 1; j < loadings.length; j++) {
                        latentscov.set(i, j, values[count]);
                        latentscov.set(j, i, values[count]);
                        count++;
                    }
                }
            }

            if (!MatrixUtils.isPositiveDefinite(latentscov)) {
                return 1e30;
            }

            if (loadings != null) {
                for (int i = 0; i < loadings.length; i++) {
                    for (int j = 0; j < loadings[i].length; j++) {
                        loadings[i][j] = values[count];
                        count++;
                    }
                }
            }

            if (delta != null) {
                for (int i = 0; i < delta.length; i++) {
                    delta[i] = values[count];
                    count++;
                }
            }

            if (quick) {
                return sumOfDifferences(indicatorIndices, measurescov, loadings, latentscov);
            } else {
                TetradMatrix implied = impliedCovariance(indicatorIndices, loadings, measurescov, latentscov, delta);

                TetradMatrix I = TetradMatrix.identity(implied.rows());
                TetradMatrix diff = I.minus((implied.times(measuresCovInverse)));

                return 0.5 * (diff.times(diff)).trace();
            }
        }
    }

    private class FunctionDag implements org.apache.commons.math3.analysis.MultivariateFunction {
        private final int[][] indicatorIndices;
        private final TetradMatrix measurescov;
        private final boolean[][] mask;
        private final boolean[][] maskErr;
        private TetradMatrix latentscov;
        private final TetradMatrix latentCoef;
        private final TetradMatrix latentErrCovar;
        private TetradMatrix measuresCovInverse;
        private final double[][] loadings;
        private final double[] delta;
        private final boolean quick;

        public FunctionDag(TetradMatrix measurescov, int[][] indicatorIndices, double[][] loadings,
                           boolean[][] mask, TetradMatrix latentCoef,
                           boolean[][] maskErr, TetradMatrix latentErrCovar,
                           double[] delta, boolean quick) {
            this.indicatorIndices = indicatorIndices;
            this.measurescov = measurescov;
            this.loadings = loadings;
            this.mask = mask;
            this.latentCoef = latentCoef;
            this.maskErr = maskErr;
            this.latentErrCovar = latentErrCovar;
            this.delta = delta;
            this.measuresCovInverse = measurescov.inverse();
            this.quick = quick;
            latentscov = MatrixUtils.impliedCovar(latentCoef, latentErrCovar);
        }

        @Override
        public double value(double[] values) {
            int count = 0;

            if (maskErr != null) {
                for (int i = 0; i < maskErr.length; i++) {
                    if (values[i] <= 0) return 1e30;
                }

                for (int i = 0; i < maskErr.length; i++) {
                    if (maskErr[i][i]) {
                        latentErrCovar.set(i, i, values[count]);
                        count++;
                    }
                }

                for (int i = 0; i < mask.length; i++) {
                    for (int j = 0; j < mask.length; j++) {
                        if (mask[i][j]) {
                            latentCoef.set(i, j, values[count]);
                            count++;
                        }
                    }
                }
            }

            latentscov = MatrixUtils.impliedCovar(latentCoef, latentErrCovar);

            if (loadings != null) {
                for (int i = 0; i < loadings.length; i++) {
                    for (int j = 0; j < loadings[i].length; j++) {
                        loadings[i][j] = values[count];
                        count++;
                    }
                }
            }

            if (delta != null)

            {
                for (int i = 0; i < delta.length; i++) {
                    delta[i] = values[count];
                    count++;
                }
            }

            if (quick)

            {
                return sumOfDifferences(indicatorIndices, measurescov, loadings, latentscov);
            } else

            {
                TetradMatrix implied = impliedCovariance(indicatorIndices, loadings, measurescov, latentscov, delta);
                TetradMatrix I = TetradMatrix.identity(implied.rows());
                TetradMatrix diff = I.minus((implied.times(measuresCovInverse)));
                return 0.5 * (diff.times(diff)).trace();
            }
        }
    }

    private TetradMatrix impliedCovariance(int[][] indicatorIndices, double[][] loadings,
                                           TetradMatrix measuredCov, TetradMatrix loadingscov,
                                           double[] delta) {
        TetradMatrix implied = new TetradMatrix(measuredCov.rows(), measuredCov.columns());

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

    private class CovRet {
        private TetradMatrix latentscov;
        private double[][] loadings;
        private int[][] indicatorIndices;
        private double[] deltas;
        private double gls;
        private List<Node> variables;

        public CovRet(TetradMatrix latentscov, double[][] loadings, int[][] indicatorIndices,
                      double[] deltas, double gls, List<Node> variables) {
            this.latentscov = latentscov;
            this.loadings = loadings;
            this.indicatorIndices = indicatorIndices;
            this.deltas = deltas;
            this.gls = gls;
            this.variables = variables;
        }

        public TetradMatrix getLatentscov() {
            return latentscov;
        }

        public double[][] getLoadings() {
            return loadings;
        }

        public double[] getDeltas() {
            return deltas;
        }

        public double getGls() {
            return gls;
        }

        public List<Node> getVariables() {
            return variables;
        }

        public int[][] getIndicatorIndices() {
            return indicatorIndices;
        }
    }
}



