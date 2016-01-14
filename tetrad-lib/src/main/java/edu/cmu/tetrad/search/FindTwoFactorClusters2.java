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
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.ChoiceGenerator;
import edu.cmu.tetrad.util.MathUtils;
import edu.cmu.tetrad.util.RandomUtil;
import edu.cmu.tetrad.util.TetradLogger;

import java.util.*;

import static java.lang.Math.*;


/**
 * Implements FindOneFactorCluster by Erich Kummerfeld (adaptation of a two factor
 * sextad algorithm to a one factor IntSextad algorithm).
 *
 * @author Joseph Ramsey
 */
public class FindTwoFactorClusters2 {

    private List<Integer> allVariables;

    public FindTwoFactorClusters.Algorithm getAlgorithm() {
        return algorithm;
    }

    public void setAlgorithm(FindTwoFactorClusters.Algorithm algorithm) {
        this.algorithm = algorithm;
    }

    public enum Algorithm {SAG, GAP}

    private CorrelationMatrix corr;
    // The list of all variables.
    private List<Node> variables;

    // The significance level.
    private double alpha;

    // The Delta test. Testing two sextads simultaneously.
    private DeltaSextadTest2 test;

    // The data.
    private transient DataModel dataModel;

    private List<List<Node>> clusters;

    private boolean verbose = true;
    private FindTwoFactorClusters.Algorithm algorithm = FindTwoFactorClusters.Algorithm.GAP;
    double scoreCutoff = Double.NEGATIVE_INFINITY;


    //========================================PUBLIC METHODS====================================//

    public FindTwoFactorClusters2(ICovarianceMatrix cov, FindTwoFactorClusters.Algorithm algorithm, double alpha) {
//        CovarianceMatrix covarianceMatrix = new CovarianceMatrix(cov);
        this.variables = cov.getVariables();
        this.alpha = alpha;
        this.test = new DeltaSextadTest2(cov);
        this.dataModel = cov;
        this.algorithm = algorithm;

        this.corr = new CorrelationMatrix(cov);

//        List<Double> corrList = new ArrayList<Double>();
//
//        for (int i = 0; i < covarianceMatrix.getDimension(); i++) {
//            for (int j = i + 1; j < covarianceMatrix.getDimension(); j++) {
//                double r = corr.getValue(i, j);
//                corrList.add(abs(r));
//            }
//        }

//        Collections.sort(corrList);
//        double threshold = corrList.get((int) (.9 * corrList.size()) - 1);
//
//        for (int i = 0; i < covarianceMatrix.getDimension(); i++) {
//            for (int j = 0; j < covarianceMatrix.getDimension(); j++) {
//                if (i == j) ;//covarianceMatrix.setValue(i, j, 0);
//                else {
//                    double r = corr.getValue(i, j);
//                    if (abs(r) > threshold) {//getCorrelationP(r) > alpha) {
//                        covarianceMatrix.setValue(i, j, 0); //signum(r) * covarianceMatrix.getValue(i, j));
//                    }
//                }
//            }
//        }

//        System.out.println(covarianceMatrix);
    }

    public FindTwoFactorClusters2(DataSet dataSet, FindTwoFactorClusters.Algorithm algorithm, double alpha) {
        this.variables = dataSet.getVariables();
        this.alpha = alpha;
        this.test = new DeltaSextadTest2(dataSet);

//        dataSet = dataSet.copy();
//
//        for (int i = 0; i < dataSet.getNumRows(); i++) {
//            for (int j = 0; j < dataSet.getNumColumns(); j++) {
//                dataSet.setDouble(i, j, Math.log(dataSet.getDouble(i, j) + 1));
//            }
//        }

//        this.dataModel = new CovarianceMatrix(dataSet);

        this.dataModel = dataSet;
        this.algorithm = algorithm;

        this.corr = new CorrelationMatrix(dataSet);
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

    //========================================PRIVATE METHODS====================================//

    // This is the main algorithm.
    private Set<List<Integer>> estimateClustersGAP() {
        List<Integer> _variables = allVariables();

        final Map<List<Integer>, Double> scoredPentads = findPurePentads(_variables);
        List<List<Integer>> pentads = new ArrayList<>(scoredPentads.keySet());

//        Collections.sort(pentads, new Comparator<List<Integer>>() {
//            @Override
//            public int compare(List<Integer> o1, List<Integer> o2) {
//                return -Double.compare(scoredPentads.get(o1), scoredPentads.get(o2));
//            }
//        });
//
//        List<List<Integer>> sublist = new ArrayList<>();
//
//        for (int i = 0; i < (int) (pentads.size() * 0.8); i++) {
//            sublist.add(pentads.get(i));
//        }
//
//        pentads = sublist;

        Set<List<Integer>> combined = combinePurePentads(pentads, _variables);

        Set<List<Integer>> _combined = new HashSet<>();

        for (List<Integer> c : combined) {
            List<Integer> a = new ArrayList<>(c);
            Collections.sort(a);
            _combined.add(a);
        }

        return _combined;
    }

    private List<Integer> allVariables() {
        if (this.allVariables != null) {
            return this.allVariables;
        }

        List<Integer> _variables = new ArrayList<>();
        for (int i = 0; i < variables.size(); i++) _variables.add(i);

        this.allVariables = _variables;
        return allVariables;
    }

    private Set<List<Integer>> estimateClustersSAG() {
        List<Integer> _variables = new ArrayList<>(allVariables());
        Set<List<Integer>> clusters1 = new HashSet<>();

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

//                    if (score(cluster) < scoreCutoff) continue;

//                    if (zeroCorr(cluster, 4)) continue;

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

                        clusters1.add(nodes);
                        _variables.removeAll(nodes);

                        continue VARIABLES;
                    }
                }

                break;
            }

        }

        return clusters1;

    }

    private Map<List<Integer>, Double> findPurePentads(List<Integer> variables) {
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

            if (score(pentad) < scoreCutoff) {
                continue;
            }

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

    // Comparably only for same-size clusters.
    private double score(List<Integer> cluster) {
//        if (true) return 10000;

        double score1 = 0.0;

//        for (int i = 0; i < cluster.size(); i++) {
//            for (int j = i + 1; j < cluster.size(); j++) {
//                double score0 = Math.log(abs(corr.getValue(cluster.get(i), cluster.get(j))));
//                score1 += score0;
//            }
//        }

        double score2 = 0.0;

        ChoiceGenerator gen = new ChoiceGenerator(cluster.size(), 5);
        int[] choice;

        while ((choice = gen.next()) != null) {
            int n1 = cluster.get(choice[0]);
            int n2 = cluster.get(choice[1]);
            int n3 = cluster.get(choice[2]);
            int n4 = cluster.get(choice[3]);
            int n5 = cluster.get(choice[4]);

            List<Integer> pentad = cluster(n1, n2, n3, n4, n5);

            for (int o : allVariables()) {
                if (pentad.contains(o)) {
                    continue;
                }

                List<Integer> sextad = new ArrayList<>(pentad);
                sextad.add(o);

                IntSextad intSextad = new IntSextad(sextad.get(0), sextad.get(1), sextad.get(2),
                        sextad.get(3), sextad.get(4), sextad.get(5));

                double pValue = test.getPValue(Collections.singletonList(intSextad));
                    score2 += Math.log(pValue);
            }
        }

        return score2;

//        return (score1 + score2) / 2;
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

    private Set<List<Integer>> combinePurePentads(List<List<Integer>> purePentads, List<Integer> _variables) {
        log("Growing pure pentads.");
        Set<List<Integer>> grown = new HashSet<>();
        List<Integer> t = new ArrayList<>();

        // Lax grow phase with speedup.
        if (true) {
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

        // Lax grow phase without speedup.
        if (false) {
            int count = 0;
            int total = purePentads.size();

            // Optimized lax version of grow phase.
            for (List<Integer> cluster : new HashSet<>(purePentads)) {
                List<Integer> _cluster = new ArrayList<>(cluster);

                for (int o : _variables) {
                    if (_cluster.contains(o)) continue;

                    List<Integer> _cluster2 = new ArrayList<>(_cluster);
                    int rejected = 0;
                    int accepted = 0;

                    ChoiceGenerator gen = new ChoiceGenerator(_cluster2.size(), 6);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        int n1 = _cluster2.get(choice[0]);
                        int n2 = _cluster2.get(choice[1]);
                        int n3 = _cluster2.get(choice[2]);
                        int n4 = _cluster2.get(choice[3]);

                        List<Integer> pentad = cluster(n1, n2, n3, n4, o);

                        if (!purePentads.contains(pentad)) {
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

                for (List<Integer> c : new HashSet<>(purePentads)) {
                    if (_cluster.containsAll(c)) {
                        purePentads.remove(c);
                    }
                }

                if (verbose) {
                    log("Grown " + (++count) + " of " + total + ": " + variablesForIndices(_cluster));
                }

                grown.add(_cluster);
            }
        }

        // Strict grow phase.
        if (false) {
            int count = 0;
            int total = purePentads.size();

            do {
                if (!purePentads.iterator().hasNext()) {
                    break;
                }

                List<Integer> cluster = purePentads.iterator().next();
                List<Integer> _cluster = new ArrayList<>(cluster);

                VARIABLES:
                for (int o : _variables) {
                    if (_cluster.contains(o)) continue;

                    List<Integer> _cluster2 = new ArrayList<>(_cluster);

                    ChoiceGenerator gen = new ChoiceGenerator(_cluster2.size(), 4);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        int n1 = _cluster2.get(choice[0]);
                        int n2 = _cluster2.get(choice[1]);
                        int n3 = _cluster2.get(choice[2]);
                        int n4 = _cluster2.get(choice[3]);

                        t.clear();
                        t.add(n1);
                        t.add(n2);
                        t.add(n3);
                        t.add(n4);
                        t.add(o);

                        if (!purePentads.contains(t)) {
                            continue VARIABLES;
                        }
                    }

                    _cluster.add(o);
                }

//                for (Set<Integer> c : new HashSet<>(purePentads)) {
////                    for (Integer d : c) {
////                        if (_cluster.contains(d)) {
////                            purePentads.remove(c);
////                        }
////                    }
//
//                    if (_cluster.containsAll(c)) {
//                        purePentads.remove(c);
//                    }
//                }
                removeSubsets(purePentads, _cluster);

                if (verbose) {
                    log("Grown " + (++count) + " of " + total + ": " + _cluster);
                }

                grown.add(_cluster);
            } while (!purePentads.isEmpty());
        }

        // Optimized pick phase.
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

//            if (tier.get(0).size() == 5) {
//                Collections.sort(tier, new Comparator<List<Integer>>() {
//                    @Override
//                    public int compare(List<Integer> o1, List<Integer> o2) {
//                        return -Double.compare(score(o1), score(o2));
//                    }
//                });
//            }

            C2:
            for (List<Integer> cluster : new ArrayList<>(tier)) {
                for (Integer i : cluster) {
                    if (all.contains(i)) {
                        continue C2;
                    }
                }

                out.add(cluster);
                all.addAll(cluster);
            }
        }

        return out;
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

                if (score(pentad) < scoreCutoff) {
                    continue O;
                }

//                if (zeroCorr(pentad, 4)) continue;

                if (!purePentad(pentad)) {
                    continue O;
                }
            }

            log("Extending by " + variables.get(o));
            cluster.add(o);
        }
    }

    private void addOtherVariablesLax(List<Integer> _variables, List<Integer> cluster) {

        O:
        for (int o : _variables) {
            if (cluster.contains(o)) continue;
            List<Integer> _cluster = new ArrayList<>(cluster);

            int accepted = 0;

            ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 4);
            int threshold = MathUtils.choose(_cluster.size(), 4);
            int[] choice;

            while ((choice = gen2.next()) != null) {
                int t1 = _cluster.get(choice[0]);
                int t2 = _cluster.get(choice[1]);
                int t3 = _cluster.get(choice[2]);
                int t4 = _cluster.get(choice[3]);

                List<Integer> sextad = cluster(t1, t2, t3, t4);
                sextad.add(o);

//                if (score(sextad) < scoreCutoff) continue;
//                if (zeroCorr(sextad, 4)) continue;

                if (!purePentad(sextad)) {
//                    rejected++;
                    continue;
                }

                accepted++;

                if (accepted > threshold) {
                    break;
                }
            }

            if (accepted >= (int)round(threshold * 0.5)) {
                log("Extending by " + variables.get(o));
                cluster.add(o);
            }
        }
    }

    private void addOtherVariablesLax2(List<Integer> _variables, List<Integer> cluster) {
        for (int o : _variables) {
            if (cluster.contains(o)) continue;
            List<Integer> _cluster = new ArrayList<>(cluster);

            int[] missed = new int[_cluster.size()];

            ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 4);
            int[] choice;

            while ((choice = gen2.next()) != null) {
                int t1 = _cluster.get(choice[0]);
                int t2 = _cluster.get(choice[1]);
                int t3 = _cluster.get(choice[2]);
                int t4 = _cluster.get(choice[3]);

                List<Integer> pentad = cluster(t1, t2, t3, t4);
                pentad.add(o);

                if (!purePentad(pentad)) {
                    for (int h : choice) {
                        missed[h]++;
                    }
                }
            }

            System.out.println("Missed " + Arrays.toString(missed) + " " + variablesForIndices(_cluster) + " " + variables.get(o));

            boolean allZero = getSum(missed) == 0;

            List<Integer> keep = new ArrayList<>();
            List<Integer> lose = new ArrayList<>();

            if (!allZero) {
                for (int i = 0; i < missed.length; i++) {
                    if (missed[i] == 0) {
                        keep.add(_cluster.get(i));
                    }
                    else if (missed[i] == MathUtils.choose(_cluster.size() - 1, 4)) {
                        lose.add(_cluster.get(i));
                    }
                }
            }

            if (!keep.isEmpty()) {
                for (int i = 0; i < keep.size(); i++) {
                    if (!cluster.contains(keep.get(i))) {
                        cluster.add(keep.get(i));
                        log("Keep Extending by " + variables.get(keep.get(i)));
                    }
                }
//                for (int i = 0; i < lose.size(); i++) {
//                    if (cluster.contains(lose.get(i))) {
//                        cluster.remove(lose.get(i));
//                    }
//                }
            }
            else if (allZero) {
                log("Extending by " + variables.get(o));
                cluster.add(o);
            }
        }
    }

    private int getSum(int[] missed) {
        int sum = 0;

        for (int m : missed) {
            sum += m;
        }
        return sum;
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

    private boolean zeroCorr(List<Integer> cluster, int n) {
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
        double p = pValue(sextad);
//        System.out.println(p);
        return p > alpha;
    }

    private double pValue(List<Integer> sextad) {
        List<IntSextad> sextads = getIntSextads(sextad.get(0), sextad.get(1), sextad.get(2),
                sextad.get(3), sextad.get(4), sextad.get(5));
        return test.getPValue(sextads);
    }

    private List<IntSextad> getIntSextads(int... n) {
        if (n.length != 6) throw new IllegalArgumentException();

        IntSextad t1 = new IntSextad(n[0], n[1], n[2], n[3], n[4], n[5]);
        IntSextad t2 = new IntSextad(n[0], n[4], n[5], n[1], n[2], n[3]);
        IntSextad t3 = new IntSextad(n[0], n[3], n[5], n[1], n[2], n[4]);
        IntSextad t4 = new IntSextad(n[0], n[3], n[4], n[1], n[2], n[5]);
        IntSextad t5 = new IntSextad(n[0], n[2], n[3], n[1], n[4], n[5]);
        IntSextad t6 = new IntSextad(n[0], n[2], n[4], n[1], n[3], n[5]);
        IntSextad t7 = new IntSextad(n[0], n[2], n[5], n[1], n[3], n[4]);
        IntSextad t8 = new IntSextad(n[0], n[1], n[3], n[2], n[4], n[5]);
        IntSextad t9 = new IntSextad(n[0], n[1], n[4], n[2], n[3], n[5]);
        IntSextad t10 = new IntSextad(n[0], n[1], n[5], n[2], n[3], n[4]);

        // The four sextads implied by equation 5.17 in Harmann.
//        return sextadList(t3, t7, t8, t9);

//        return sextadList(t1, t8, t9, t10, t5);

//            IntSextad[] independents = {t2, t5, t10, t3, t6};

//        return sextadList(t1, t2, t3, t5, t6);
//        return sextadList(t1, t2, t3, t9, t10);
//        return sextadList(t6, t7, t8, t9, t10);
//        return sextadList(t1, t2, t4, t5, t9);
//        return sextadList(t1, t3, t4, t6, t10);
        return sextadList(t6, t7, t8, t9, t10);

//        List<IntSextad> intSextads = sextadList(t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);

//        List<IntSextad> sextads = new ArrayList<>(choice.length);
//
//        for (int i = 0; i < choice.length; i++) {
//            sextads.add(intSextads.get(choice[i]));
//        }
//
//        return sextads;
//
//        return intSextads;

    }

    private List<IntSextad> sextadList(IntSextad... t) {
        List<IntSextad> list = new ArrayList<>();

        for (int i = 0; i < t.length; i++) {
            list.add(t[i]);
        }

        return list;
    }

    private double logsum(List<Integer> cluster, double pValue) {
        double logsum = 0.0;

        for (int i = 0; i < cluster.size(); i++) {
            for (int j = i + 1; j < cluster.size(); j++) {
                double r = this.corr.getValue(cluster.get(i), cluster.get(j));
                if (abs(r) > 0) {
                    logsum += Math.log(abs(r));
                }
            }
        }

        if (pValue > 0) {
            logsum += Math.log(pValue);
        }

        return logsum;
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
}




