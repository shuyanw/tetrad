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

import edu.cmu.tetrad.data.CorrelationMatrix;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.*;

import static java.lang.Math.abs;
import static java.lang.Math.sqrt;


/**
 * Implements FindTwoFactorClusters by Erich Kummerfeld.
 *
 * @author Joseph Ramsey
 */
public class FindTwoFactorClusters2 {

    public boolean isZeroCorrChecked() {
        return zeroCorrChecked;
    }

    public void setZeroCorrChecked(boolean zeroCorrChecked) {
        this.zeroCorrChecked = zeroCorrChecked;
    }

    public boolean isScoreChecked() {
        return scoreChecked;
    }

    public void setScoreChecked(boolean scoreChecked) {
        this.scoreChecked = scoreChecked;
    }

    public enum Algorithm {SAG, GAP}

    public void setAlgorithm(FindTwoFactorClusters.Algorithm algorithm) {
        this.algorithm = algorithm;
    }

    // Correlation matrix, for scoring.
    private CorrelationMatrix corr;

    // The list of all variables.
    private List<Node> variables;

    // Indices for all of the variables.
    private List<Integer> variableIndices;

    // True if pairwise zero correlations should be checked.
    private boolean zeroCorrChecked = false;

    // True if scores should be checked.
    private boolean scoreChecked = false;

    // Lower bound on the average pairwise correlation score.
    double scoreCutoff = Double.NEGATIVE_INFINITY;

    // The significance level.
    private double alpha;

    // The Delta test. Testing two sextads simultaneously.
    private DeltaSextadTest2 test;

    private List<List<Node>> clusters;

    // True if verbose output should be printed.
    private boolean verbose = true;

    // The algorithm to run, GAP or SAP.
    private FindTwoFactorClusters.Algorithm algorithm = FindTwoFactorClusters.Algorithm.GAP;

    //========================================PUBLIC METHODS====================================//

    public FindTwoFactorClusters2(ICovarianceMatrix cov, FindTwoFactorClusters.Algorithm algorithm, double alpha) {
        this.variables = cov.getVariables();
        this.alpha = alpha;
        this.test = new DeltaSextadTest2(cov);
        this.algorithm = algorithm;
        this.corr = new CorrelationMatrix(cov);
        this.scoreCutoff = Math.log(alpha);
    }

    public FindTwoFactorClusters2(DataSet dataSet, FindTwoFactorClusters.Algorithm algorithm, double alpha) {
        this.variables = dataSet.getVariables();
        this.alpha = alpha;
        this.test = new DeltaSextadTest2(dataSet);
        this.algorithm = algorithm;
        this.corr = new CorrelationMatrix(dataSet);
        this.scoreCutoff = Math.log(alpha);
    }

    public Graph search() {
        Set<List<Integer>> allClusters;

        if (algorithm == FindTwoFactorClusters.Algorithm.GAP) {
            allClusters = estimateClustersGAP();
        } else if (algorithm == FindTwoFactorClusters.Algorithm.SAG) {
            allClusters = estimateClustersSAG();
        } else {
            throw new IllegalStateException("Expected GAP or SAG: " + algorithm);
        }

        this.clusters = variablesForIndices(allClusters);
        return convertToGraph(allClusters);
    }

    public FindTwoFactorClusters.Algorithm getAlgorithm() {
        return algorithm;
    }

    //========================================PRIVATE METHODS====================================//

    private Set<List<Integer>> estimateClustersSAG() {
        List<Integer> _variables = new ArrayList<>(allVariables());
        Set<List<Integer>> clusters = new HashSet<>();

        int clusterSize = 5;

        for (int k = 6; k >= clusterSize; k--) {

            VARIABLES:
            while (!_variables.isEmpty()) {
                if (verbose) {
                    log(_variables.toString());
                }

                if (_variables.size() < clusterSize) break;

                ChoiceGenerator gen = new ChoiceGenerator(_variables.size(), clusterSize);
                int[] choice;

                while ((choice = gen.next()) != null) {
                    int n1 = _variables.get(choice[0]);
                    int n2 = _variables.get(choice[1]);
                    int n3 = _variables.get(choice[2]);
                    int n4 = _variables.get(choice[3]);
                    int n5 = _variables.get(choice[4]);

                    List<Integer> cluster = cluster(n1, n2, n3, n4, n5);

                    if (scoreFails(cluster)) continue;
                    if (zeroCorr(cluster, 4)) continue;

                    // Note that purity needs to be assessed with respect to all of the variables in order to
                    // remove all latent-measure impurities between pairs of latents.
                    if (purePentad(cluster)) {
                        if (verbose) {
                            log("Found a pure: " + variablesForIndices(cluster) + " score = " + score(cluster));
                        }

                        List<Integer> nodes = new ArrayList<>(cluster);

                        addOtherVariablesStrict(_variables, nodes);

                        if (nodes.size() < k) continue;

                        if (verbose) {
                            log("Cluster found: " + variablesForIndices(cluster));
                        }

                        clusters.add(nodes);
                        _variables.removeAll(nodes);

                        continue VARIABLES;
                    }
                }

                break;
            }
        }

        return clusters;

    }

    private void addOtherVariablesStrict(List<Integer> _variables, List<Integer> cluster) {

        O:
        for (int o : _variables) {
            if (cluster.contains(o)) continue;
            List<Integer> _cluster = new ArrayList<>(cluster);

            ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 4);
            int[] choice;

            while ((choice = gen2.next()) != null) {
                int t1 = _cluster.get(choice[0]);
                int t2 = _cluster.get(choice[1]);
                int t3 = _cluster.get(choice[2]);
                int t4 = _cluster.get(choice[3]);

                List<Integer> pentad = cluster(t1, t2, t3, t4);
                pentad.add(o);

                if (scoreFails(pentad)) {
                    continue O;
                }

                if (zeroCorr(pentad, 4)) continue;

                if (!purePentad(pentad)) {
                    continue O;
                }
            }

            log("Extending by " + variables.get(o));
            cluster.add(o);
        }
    }

    private Set<List<Integer>> estimateClustersGAP() {
        List<Integer> _variables = allVariables();

        final Map<List<Integer>, Double> scoredPentads = findPurePentadsGAP(_variables);
        List<List<Integer>> pentads = new ArrayList<>(scoredPentads.keySet());
        Set<List<Integer>> combined = combinePurePentadsGAP(pentads, _variables);
        Set<List<Integer>> _combined = new HashSet<>();

        for (List<Integer> c : combined) {
            List<Integer> a = new ArrayList<>(c);
            Collections.sort(a);
            _combined.add(a);
        }

        return _combined;
    }

    private Map<List<Integer>, Double> findPurePentadsGAP(List<Integer> variables) {
        if (variables.size() < 5) {
            return new HashMap<>();
        }

        log("Finding pure pentads.");

        ChoiceGenerator gen = new ChoiceGenerator(variables.size(), 5);
        int[] choice;
        Map<List<Integer>, Double> purePentads = new HashMap<>();

        while ((choice = gen.next()) != null) {
            int n1 = variables.get(choice[0]);
            int n2 = variables.get(choice[1]);
            int n3 = variables.get(choice[2]);
            int n4 = variables.get(choice[3]);
            int n5 = variables.get(choice[4]);

            List<Integer> pentad = cluster(n1, n2, n3, n4, n5);

            if (scoreFails(pentad)) continue;
            if (zeroCorr(pentad, 4)) continue;

            if (!purePentad(pentad)) {
                continue;
            }

            List<Integer> _cluster = new ArrayList<>(pentad);

            if (verbose) {
                log("++ " + variablesForIndices(pentad));
            }


            double score = score(pentad);
            purePentads.put(_cluster, score);
        }

        return purePentads;
    }

    private Set<List<Integer>> combinePurePentadsGAP(List<List<Integer>> purePentads, List<Integer> _variables) {
        log("Growing pure pentads.");
        Set<List<Integer>> grown = new HashSet<>();
        List<Integer> t = new ArrayList<>();
        laxGrowPhaseGAP(purePentads, _variables, grown, t);
        return optimizedPickPhase(grown);
    }

    private void laxGrowPhaseGAP(List<List<Integer>> purePentads, List<Integer> _variables, Set<List<Integer>> grown, List<Integer> t) {
        int count = 0;
        int total = purePentads.size();

        do {
            if (!purePentads.iterator().hasNext()) {
                break;
            }

            List<Integer> cluster = purePentads.iterator().next();
            List<Integer> _cluster = new ArrayList<>(cluster);

            for (int o : _variables) {
                if (_cluster.contains(o)) continue;

                List<Integer> _cluster2 = new ArrayList<>(_cluster);
                int rejected = 0;
                int accepted = 0;

                ChoiceGenerator gen = new ChoiceGenerator(_cluster2.size(), 4);
                int[] choice;

                while ((choice = gen.next()) != null) {
                    t.clear();
                    t.add(_cluster2.get(choice[0]));
                    t.add(_cluster2.get(choice[1]));
                    t.add(_cluster2.get(choice[2]));
                    t.add(_cluster2.get(choice[3]));
                    t.add(o);

                    if (!purePentads.contains(t)) {
                        rejected++;
                    } else {
                        accepted++;
                    }
                }

                if (rejected > accepted) {
                    continue;
                }

                _cluster.add(o);
            }

            removeSubsets(purePentads, _cluster);

            if (verbose) {
                log("Grown " + (++count) + " of " + total + ": " +
                        variablesForIndices(new ArrayList<Integer>(_cluster)));
            }
            grown.add(_cluster);
        } while (!purePentads.isEmpty());
    }

//    private void strictGrowPhaseGAP(List<List<Integer>> purePentads, List<Integer> _variables, Set<List<Integer>> grown, List<Integer> t) {
//        int count = 0;
//        int total = purePentads.size();
//
//        do {
//            if (!purePentads.iterator().hasNext()) {
//                break;
//            }
//
//            List<Integer> cluster = purePentads.iterator().next();
//            List<Integer> _cluster = new ArrayList<>(cluster);
//
//            VARIABLES:
//            for (int o : _variables) {
//                if (_cluster.contains(o)) continue;
//
//                List<Integer> _cluster2 = new ArrayList<>(_cluster);
//
//                ChoiceGenerator gen = new ChoiceGenerator(_cluster2.size(), 4);
//                int[] choice;
//
//                while ((choice = gen.next()) != null) {
//                    int n1 = _cluster2.get(choice[0]);
//                    int n2 = _cluster2.get(choice[1]);
//                    int n3 = _cluster2.get(choice[2]);
//                    int n4 = _cluster2.get(choice[3]);
//
//                    t.clear();
//                    t.add(n1);
//                    t.add(n2);
//                    t.add(n3);
//                    t.add(n4);
//                    t.add(o);
//
//                    if (!purePentads.contains(t)) {
//                        continue VARIABLES;
//                    }
//                }
//
//                _cluster.add(o);
//            }
//
//            removeSubsets(purePentads, _cluster);
//
//            if (verbose) {
//                log("Grown " + (++count) + " of " + total + ": " + _cluster);
//            }
//
//            grown.add(_cluster);
//        } while (!purePentads.isEmpty());
//    }

//    private void laxGrowPhaseWithoutSpeedupGAP(List<List<Integer>> purePentads, List<Integer> _variables, Set<List<Integer>> grown) {
//        int count = 0;
//        int total = purePentads.size();
//
//        // Optimized lax version of grow phase.
//        for (List<Integer> cluster : new HashSet<>(purePentads)) {
//            List<Integer> _cluster = new ArrayList<>(cluster);
//
//            for (int o : _variables) {
//                if (_cluster.contains(o)) continue;
//
//                List<Integer> _cluster2 = new ArrayList<>(_cluster);
//                int rejected = 0;
//                int accepted = 0;
//
//                ChoiceGenerator gen = new ChoiceGenerator(_cluster2.size(), 6);
//                int[] choice;
//
//                while ((choice = gen.next()) != null) {
//                    int n1 = _cluster2.get(choice[0]);
//                    int n2 = _cluster2.get(choice[1]);
//                    int n3 = _cluster2.get(choice[2]);
//                    int n4 = _cluster2.get(choice[3]);
//
//                    List<Integer> pentad = cluster(n1, n2, n3, n4, o);
//
//                    if (!purePentads.contains(pentad)) {
//                        rejected++;
//                    } else {
//                        accepted++;
//                    }
//                }
//
//                if (rejected > accepted) {
//                    continue;
//                }
//
//                _cluster.add(o);
//            }
//
//            for (List<Integer> c : new HashSet<>(purePentads)) {
//                if (_cluster.containsAll(c)) {
//                    purePentads.remove(c);
//                }
//            }
//
//            if (verbose) {
//                log("Grown " + (++count) + " of " + total + ": " + variablesForIndices(_cluster));
//            }
//
//            grown.add(_cluster);
//        }
//    }

    private Set<List<Integer>> optimizedPickPhase(Set<List<Integer>> grown) {
        log("Choosing among grown clusters.");

        List<List<Integer>> list = new ArrayList<>(grown);

        Collections.sort(list, new Comparator<List<Integer>>() {
            @Override
            public int compare(List<Integer> o1, List<Integer> o2) {
                return o2.size() - o1.size();
            }
        });

        List<Integer> all = new ArrayList<>();
        Set<List<Integer>> out = new HashSet<>();

        while (!(new ArrayList<>(list).isEmpty())) {
            List<List<Integer>> tier = getTier(list, list.get(0).size());
            list.removeAll(tier);

            C:
            for (List<Integer> cluster : new ArrayList<>(tier)) {
                for (Integer i : cluster) {
                    if (all.contains(i)) {
                        continue C;
                    }
                }

                out.add(cluster);
                all.addAll(cluster);
            }
        }

        return out;
    }

    // Comparable only for same-size clusters.
    private double score(List<Integer> cluster) {
        double score = 0.0;
        int n = 0;

        for (int i = 0; i < cluster.size(); i++) {
            for (int j = i + 1; j < cluster.size(); j++) {
                double score0 = Math.log(abs(corr.getValue(cluster.get(i), cluster.get(j))));
                score += score0;
                n++;
            }
        }

        score /= n;
        return score;
    }

    private boolean purePentad(List<Integer> pentad) {
        for (int o : allVariables()) {
            if (pentad.contains(o)) {
                continue;
            }

            List<Integer> sextad = new ArrayList<>(pentad);
            sextad.add(o);

            if (!vanishes(sextad)) {
                return false;
            }
        }

        return true;
    }

    private List<List<Integer>> getTier(List<List<Integer>> list, int n) {
        List<List<Integer>> tier = new ArrayList<>();

        for (List<Integer> aList : list) {
            if (aList.size() == n) tier.add(aList);
        }

        return tier;
    }

    private void removeSubsets(List<List<Integer>> purePentads, List<Integer> _cluster) {
        List<Integer> t = new ArrayList<>();

        // This takes out all pure clusters that are subsets of _cluster.
        ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 5);
        int[] choice2;
        List<Integer> _cluster3 = new ArrayList<>(_cluster);

        while ((choice2 = gen2.next()) != null) {
            int n1 = _cluster3.get(choice2[0]);
            int n2 = _cluster3.get(choice2[1]);
            int n3 = _cluster3.get(choice2[2]);
            int n4 = _cluster3.get(choice2[3]);
            int n5 = _cluster3.get(choice2[4]);

            t.clear();
            t.add(n1);
            t.add(n2);
            t.add(n3);
            t.add(n4);
            t.add(n5);

            purePentads.remove(t);
        }
    }

    private List<Node> variablesForIndices(List<Integer> cluster) {
        List<Node> _cluster = new ArrayList<>();

        for (int c : cluster) {
            _cluster.add(variables.get(c));
        }

        return _cluster;
    }

    private List<List<Node>> variablesForIndices(Set<List<Integer>> clusters) {
        List<List<Node>> variables = new ArrayList<>();

        for (List<Integer> cluster : clusters) {
            variables.add(variablesForIndices(cluster));
        }

        return variables;
    }

    private List<Integer> cluster(int... vars) {
        List<Integer> cluster = new ArrayList<>();

        for (int i : vars) {
            cluster.add(i);
        }

        if (new HashSet<>(cluster).size() < vars.length)
            throw new IllegalArgumentException("cluster elements must be unique: <" + cluster + ">");

        return cluster;
    }

    private boolean scoreFails(List<Integer> cluster) {
        if (!scoreChecked) return false;
        return score(cluster) < scoreCutoff;
    }

    private boolean zeroCorr(List<Integer> cluster, int n) {
        if (!zeroCorrChecked) return false;

        int count = 0;

        for (int i = 0; i < cluster.size(); i++) {
            for (int j = i + 1; j < cluster.size(); j++) {
                double r = this.corr.getValue(cluster.get(i), cluster.get(j));
                double p = getCorrelationP(r);
                if (p > alpha) count++;
            }
        }

        return count >= n;
    }

    private double getCorrelationP(double r) {
        int N = this.corr.getSampleSize();
        double f = sqrt(N) * Math.log((1. + r) / (1. - r));
        return 2.0 * (1.0 - RandomUtil.getInstance().normalCdf(0, 1, abs(f)));
    }

    /**
     * The clusters output by the algorithm from the last call to search().
     */
    public List<List<Node>> getClusters() {
        return clusters;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    private boolean vanishes(List<Integer> sextad) {
        return pValue(sextad) > alpha;
    }

    private double pValue(List<Integer> sextad) {
        List<IntSextad> sextads = getIntSextads(sextad.get(0), sextad.get(1), sextad.get(2),
                sextad.get(3), sextad.get(4), sextad.get(5));
        return test.getPValue(sextads);
    }

    private List<IntSextad> getIntSextads(int... n) {
        if (n.length != 6) throw new IllegalArgumentException();

//        IntSextad t1 = new IntSextad(n[0], n[1], n[2], n[3], n[4], n[5]);
//        IntSextad t2 = new IntSextad(n[0], n[4], n[5], n[1], n[2], n[3]);
//        IntSextad t3 = new IntSextad(n[0], n[3], n[5], n[1], n[2], n[4]);
//        IntSextad t4 = new IntSextad(n[0], n[3], n[4], n[1], n[2], n[5]);
//        IntSextad t5 = new IntSextad(n[0], n[2], n[3], n[1], n[4], n[5]);
        IntSextad t6 = new IntSextad(n[0], n[2], n[4], n[1], n[3], n[5]);
        IntSextad t7 = new IntSextad(n[0], n[2], n[5], n[1], n[3], n[4]);
        IntSextad t8 = new IntSextad(n[0], n[1], n[3], n[2], n[4], n[5]);
        IntSextad t9 = new IntSextad(n[0], n[1], n[4], n[2], n[3], n[5]);
        IntSextad t10 = new IntSextad(n[0], n[1], n[5], n[2], n[3], n[4]);

        // The four sextads implied by equation 5.17 in Harmann.
//        return sextadList(t3, t7, t8, t9);

//        return sextadList(t1, t2, t3, t5, t6);
//        return sextadList(t1, t2, t3, t9, t10);
//        return sextadList(t6, t7, t8, t9, t10);
//        return sextadList(t1, t2, t4, t5, t9);
//        return sextadList(t1, t3, t4, t6, t10);
        return sextadList(t6, t7, t8, t9, t10);
    }

    private List<IntSextad> sextadList(IntSextad... t) {
        List<IntSextad> list = new ArrayList<>();
        Collections.addAll(list, t);
        return list;
    }

    private Graph convertSearchGraphNodes(Set<Set<Node>> clusters) {
        Graph graph = new EdgeListGraph(variables);

        List<Node> latents = new ArrayList<>();
        for (int i = 0; i < clusters.size(); i++) {
            Node latent = new GraphNode(MimBuild.LATENT_PREFIX + (i + 1));
            latent.setNodeType(NodeType.LATENT);
            latents.add(latent);
            graph.addNode(latent);
        }

        List<Set<Node>> _clusters = new ArrayList<>(clusters);

        for (int i = 0; i < latents.size(); i++) {
            for (Node node : _clusters.get(i)) {
                if (!graph.containsNode(node)) graph.addNode(node);
                graph.addDirectedEdge(latents.get(i), node);
            }
        }

        return graph;
    }

    private Graph convertToGraph(Set<List<Integer>> allClusters) {
        Set<Set<Node>> _clustering = new HashSet<>();

        for (List<Integer> cluster : allClusters) {
            Set<Node> nodes = new HashSet<>();

            for (int i : cluster) {
                nodes.add(variables.get(i));
            }

            _clustering.add(nodes);
        }

        return convertSearchGraphNodes(_clustering);
    }

    private void log(String s) {
        if (verbose) {
            TetradLogger.getInstance().log("info", s);
            System.out.println(s);
        }
    }

    private List<Integer> allVariables() {
        if (this.variableIndices != null) {
            return this.variableIndices;
        }

        List<Integer> _variables = new ArrayList<>();
        for (int i = 0; i < variables.size(); i++) _variables.add(i);

        this.variableIndices = _variables;
        return variableIndices;
    }
}




