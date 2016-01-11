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
import edu.cmu.tetrad.sem.*;
import edu.cmu.tetrad.util.*;

import java.util.*;

import static java.lang.Math.abs;
import static java.lang.Math.signum;
import static java.lang.Math.sqrt;


/**
 * Implements FindOneFactorCluster by Erich Kummerfeld (adaptation of a two factor
 * sextet algorithm to a one factor IntSextad algorithm).
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

    private boolean verbose = false;
    private FindTwoFactorClusters.Algorithm algorithm = FindTwoFactorClusters.Algorithm.GAP;
    double scoreCutoff = Double.NEGATIVE_INFINITY;


    //========================================PUBLIC METHODS====================================//

    public FindTwoFactorClusters2(ICovarianceMatrix cov, FindTwoFactorClusters.Algorithm algorithm, double alpha) {
        CovarianceMatrix covarianceMatrix = new CovarianceMatrix(cov);
        cov = covarianceMatrix;
        this.variables = cov.getVariables();
        this.alpha = alpha;
        this.test = new DeltaSextadTest2(cov);
        this.dataModel = covarianceMatrix;
        this.algorithm = algorithm;

        this.corr = new CorrelationMatrix(cov);

        List<Double> corrList = new ArrayList<Double>();

        for (int i = 0; i < covarianceMatrix.getDimension(); i++) {
            for (int j = 0; j < covarianceMatrix.getDimension(); j++) {
                double r = corr.getValue(i, j);
                corrList.add(abs(r));
            }
        }

        Collections.sort(corrList);
        double threshold = corrList.get((int)(.9 * corrList.size()) - 1);

        for (int i = 0; i < covarianceMatrix.getDimension(); i++) {
            for (int j = 0; j < covarianceMatrix.getDimension(); j++) {
                if (i == j) covarianceMatrix.setValue(i, j, 0);
                else {
                    double r = corr.getValue(i, j);
                    if (abs(r) > threshold) {//getCorrelationP(r) > alpha) {
                        covarianceMatrix.setValue(i, j, 0);//i, j, signum(r) * threshold);
                    }
                }
            }
        }

        System.out.println(covarianceMatrix);
    }

    public FindTwoFactorClusters2(DataSet dataSet, FindTwoFactorClusters.Algorithm algorithm, double alpha) {
        this.variables = dataSet.getVariables();
        this.alpha = alpha;
        this.test = new DeltaSextadTest2(dataSet);
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

//        List<List<Integer>> sublist = new ArrayList<>();
//
//        for (int i = 0; i < (int) (pentads.size() * 0.8); i++) {
//            sublist.add(pentads.get(i));
//        }

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
        List<Integer> _variables = allVariables();

        Set<List<Integer>> pureClusters = findPureClusters(_variables);
//        Set<List<Integer>> mixedClusters = findMixedClusters(pureClusters, _variables, unionPure(pureClusters));
        Set<List<Integer>> allClusters = new HashSet<>(pureClusters);
//        allClusters.addAll(mixedClusters);
        return allClusters;
    }

    private Map<List<Integer>, Double> findPurePentads(List<Integer> variables) {
        if (variables.size() < 5) {
            return new HashMap<>();
        }

        log("Finding pure pentads.", true);

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

//            if (score(pentad) < scoreCutoff) {
//                continue;
//            }

            if (!purePentad(pentad)) {
                continue;
            }

            List<Integer> _cluster = new ArrayList<>(pentad);

//            if (verbose) {
                System.out.println("++ " + variablesForIndices(pentad));
                log("++ " + variablesForIndices(pentad), false);
//            }


            double score = score(pentad);

            purePentads.put(_cluster, score);
        }

        return purePentads;
    }

    // Comparably only for same-size clusters.
    private double score(List<Integer> cluster) {
//        if (true) return 100;

        double score1 = 0.0;

        for (int i = 0; i < cluster.size(); i++) {
            for (int j = i + 1; j < cluster.size(); j++) {
                double score0 = Math.log(abs(corr.getValue(cluster.get(i), cluster.get(j))));
                score1 += score0;
            }
        }

//        double score2 = 0.0;
//
//        ChoiceGenerator gen = new ChoiceGenerator(cluster.size(), 5);
//        int[] choice;
//
//        while ((choice = gen.next()) != null) {
//            int n1 = cluster.get(choice[0]);
//            int n2 = cluster.get(choice[1]);
//            int n3 = cluster.get(choice[2]);
//            int n4 = cluster.get(choice[3]);
//            int n5 = cluster.get(choice[4]);
//
//            List<Integer> pentad = cluster(n1, n2, n3, n4, n5);
//
//            for (int o : allVariables()) {
//                if (pentad.contains(o)) {
//                    continue;
//                }
//
//                List<Integer> sextet = new ArrayList<>(pentad);
//                sextet.add(o);
//
//                Collections.sort(sextet);
//                IntSextad sextad = new IntSextad(sextet.get(0), sextet.get(1), sextet.get(2),
//                        sextet.get(3), sextet.get(4), sextet.get(5));
//                double pValue = test.getPValue(Collections.singletonList(sextad));
//                score2 += Math.log(pValue);
//            }
//        }
//
//        double score = (score1 + score2) / 2;

//        System.out.println(score);

        return score1;
    }

    private boolean purePentad(List<Integer> pentad) {
//        if (true) return true;

        for (int o : allVariables()) {
            if (pentad.contains(o)) {
                continue;
            }

            List<Integer> sextet = new ArrayList<>(pentad);
            sextet.add(o);

            Collections.sort(sextet);

            if (!vanishes(sextet)) {
//                System.out.println("Doesnn't vanish: " + pentad);
                return false;
            }
        }

//        System.out.println("PPP " + pentad);

        return true;


    }

//    private boolean existsImpurePentad(List<Integer> variables) {
//        log("exists impure pentad.", true);
//
//        ChoiceGenerator gen = new ChoiceGenerator(variables.size(), 5);
//        int[] choice;
//
//        while ((choice = gen.next()) != null) {
//            int n1 = variables.get(choice[0]);
//            int n2 = variables.get(choice[1]);
//            int n3 = variables.get(choice[2]);
//            int n4 = variables.get(choice[3]);
//            int n5 = variables.get(choice[4]);
//
//            List<Integer> pentad = pentad(n1, n2, n3, n4, n5);
//
//            if (zeroCorr(pentad, 4)) continue;
//
//            for (int o : allVariables()) {
//                if (pentad.contains(o)) {
//                    continue;
//                }
//
//                List<Integer> sextet = sextet(n1, n2, n3, n4, n5, o);
//
//                Collections.sort(sextet);
//
//                boolean vanishes = vanishes(sextet);
//
//                if (!vanishes) {
//                    return true;
//                }
//            }
//        }
//
//        return false;
//    }


    private Set<List<Integer>> combinePurePentads(List<List<Integer>> purePentads, List<Integer> _variables) {
        log("Growing pure pentads.", true);
        Set<List<Integer>> grown = new HashSet<>();

        // Lax grow phase with speedup.
        if (false) {
            List<Integer> t = new ArrayList<>();
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

//                    if (!(avgSumLnP(new ArrayList<Integer>(_cluster)) > -10)) {
//                        _cluster.remove(o);
//                    }
                }

                // This takes out all pure clusters that are subsets of _cluster.
                ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 5);
                int[] choice2;
                List<Integer> _cluster3 = new ArrayList<Integer>(_cluster);

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

                if (verbose) {
                    System.out.println("Grown " + (++count) + " of " + total + ": " + variablesForIndices(new ArrayList<Integer>(_cluster)));
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

                        List<Integer> t = new ArrayList<>(pentad);

                        Collections.sort(t);

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

                for (List<Integer> c : new HashSet<>(purePentads)) {
                    if (_cluster.containsAll(c)) {
                        purePentads.remove(c);
                    }
                }

                if (verbose) {
                    System.out.println("Grown " + (++count) + " of " + total + ": " + _cluster);
                }

                grown.add(_cluster);
            }
        }

        // Strict grow phase.
        if (true) {
            List<Integer> t = new ArrayList<>();
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

                        Collections.sort(t);

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
                removeSubsets(purePentads, t, _cluster);

                if (verbose) {
                    System.out.println("Grown " + (++count) + " of " + total + ": " + _cluster);
                }

                grown.add(_cluster);
            } while (!purePentads.isEmpty());
        }

        // Optimized pick phase.
        log("Choosing among grown clusters.", true);

        List<List<Integer>> list = new ArrayList<>(grown);

        Collections.sort(list, new Comparator<List<Integer>>() {
            @Override
            public int compare(List<Integer> o1, List<Integer> o2) {
                return o2.size() - o1.size();
            }
        });

        System.out.println("A");

        List<Integer> all = new ArrayList<>();
        Set<List<Integer>> out = new HashSet<>();

        while (!(new ArrayList<>(list).isEmpty())) {
            List<List<Integer>> tier = getTier(list, list.get(0).size());
            list.removeAll(tier);


            System.out.println("B");

//            if (tier.get(0).size() == 5) {
//                Collections.sort(tier, new Comparator<List<Integer>>() {
//                    @Override
//                    public int compare(List<Integer> o1, List<Integer> o2) {
//                        return -Double.compare(score(o1), score(o2));
//                    }
//                });
//            }

            System.out.println("C");

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

            System.out.println("D");
        }

        System.out.println("OOUUTT = " + out);

        return out;
    }

    private List<List<Integer>> getTier(List<List<Integer>> list, int n) {
        List<List<Integer>> tier = new ArrayList<>();

        for (int i = 0; i < list.size(); i++) {
            if (list.get(i).size() == n) tier.add(list.get(i));
        }

        return tier;
    }

    private void removeSubsets(List<List<Integer>> purePentads, List<Integer> t, List<Integer> _cluster) {
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

            Collections.sort(t);

            purePentads.remove(t);
        }
    }

    // Finds clusters of size 5 or higher for the SAG algorithm.
    private Set<List<Integer>> findPureClusters(List<Integer> _variables) {
        Set<List<Integer>> clusters = new HashSet<>();

        int clusterSize = 5;

        for (int k = 6; k >= clusterSize; k--) {
            VARIABLES:
            while (!_variables.isEmpty()) {
                if (verbose) {
                    System.out.println(_variables);
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
//                    int n6 = _variables.get(choice[5]);

                    List<Integer> cluster = cluster(n1, n2, n3, n4, n5);

//                    if (score(cluster) < scoreCutoff) continue;

//                    if (zeroCorr(cluster, 4)) continue;

                    // Note that purity needs to be assessed with respect to all of the variables in order to
                    // remove all latent-measure impurities between pairs of latents.
                    if (purePentad(cluster)) {
                        if (verbose) {
                            log("Found a pure: " + variablesForIndices(cluster) + " score = " + score(cluster), true);
                        }

                        addOtherVariables1(_variables, cluster);

                        if (cluster.size() < k) continue;

                        if (verbose) {
                            log("Cluster found: " + variablesForIndices(cluster), true);
//                            System.out.println("Indices for cluster = " + cluster);
                        }

                        clusters.add(cluster);
                        _variables.removeAll(cluster);

                        continue VARIABLES;
                    }
                }

                break;
            }

        }

        return clusters;
    }

    private double logsum2(List<Integer> cluster) {
        List<List<IntSextad>> intSextads = getIntSextads(cluster.get(0), cluster.get(1), cluster.get(2),
                cluster.get(3), cluster.get(4), cluster.get(5));
        double pValue = test.getPValue(intSextads.get(0));
        double logsum = Math.log(pValue);
        logsum += test.getLogProductCorrelations(intSextads.get(0));
//        System.out.println(logsum);
        return logsum;
    }

    private void addOtherVariables1(List<Integer> _variables, List<Integer> cluster) {

        O:
        for (int o : _variables) {
            if (cluster.contains(o)) continue;
            List<Integer> _cluster = new ArrayList<>(cluster);

            int accepted = 0;
            int rejected = 0;

            ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 4);
            int threshold = ChoiceGenerator.getNumCombinations(_cluster.size(), 4) / 3;
            int[] choice;

            while ((choice = gen2.next()) != null) {
                int t1 = _cluster.get(choice[0]);
                int t2 = _cluster.get(choice[1]);
                int t3 = _cluster.get(choice[2]);
                int t4 = _cluster.get(choice[3]);

                List<Integer> sextad = cluster(t1, t2, t3, t4);
                sextad.add(o);

//                if (score(sextad) < scoreCutoff) {
//                    rejected++;
//                    continue O;
//                }
//                if (zeroCorr(sextad, 4)) continue;

                if (!purePentad(sextad)) {
//                    rejected++;
                    continue O;
                }

                accepted++;

//                if (accepted > threshold) {
//                    break;
//                }
            }

            log("Extending by " + variables.get(o), true);
            cluster.add(o);
        }
    }


//    private void addOtherVariables(List<Integer> _variables, List<Integer> cluster) {
//
//        List<Integer> clusterSaved = new ArrayList<>(cluster);
//        List<Integer> firstPass = new ArrayList<>(cluster);
//
//        O:
//        for (int o : _variables) {
//            if (clusterSaved.contains(o)) continue;
//            List<Integer> _cluster = new ArrayList<>(clusterSaved);
//
//            ChoiceGenerator gen2 = new ChoiceGenerator(clusterSaved.size(), 5);
//            int[] choice;
////            int[] counts = new int[_cluster.size()];
//
//            while ((choice = gen2.next()) != null) {
//                int t1 = clusterSaved.get(choice[0]);
//                int t2 = clusterSaved.get(choice[1]);
//                int t3 = clusterSaved.get(choice[2]);
//                int t4 = clusterSaved.get(choice[3]);
//                int t5 = clusterSaved.get(choice[4]);
//
//                List<Integer> _cluster = pentad(t1, t2, t3, t4, t5);
//                _cluster.add(o);
//
//                if (zeroCorr(_cluster, 4)) continue;
//
//                if (!pure(_cluster)) {
//                    continue O;
//                }
//            }
//
//            log("Extending by " + variables.get(o), true);
//            firstPass.add(o);
//        }
//
//        O:
//        for (int o : _variables) {
//            List<Integer> _cluster = new ArrayList<>(firstPass);
//
//            ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 4);
//            int[] choice;
////            int[] counts = new int[_cluster.size()];
//
//            while ((choice = gen2.next()) != null) {
//                int t1 = _cluster.get(choice[0]);
//                int t2 = _cluster.get(choice[1]);
//                int t3 = _cluster.get(choice[2]);
//                int t4 = _cluster.get(choice[3]);
//
//                List<Integer> pentad = quartet(t1, t2, t3, t4);
//                if (pentad.contains(o)) continue;
//                pentad.add(o);
//
//                if (!purePentad(pentad)) {
//                    continue O;
//                }
//            }
//
//            log("Extending by " + variables.get(o), true);
//            cluster.add(o);
//        }
//    }

    //  Finds clusters of size 5 for the sextet first algorithm.
    private Set<List<Integer>> findMixedClusters(Set<List<Integer>> clusters, List<Integer> remaining, Set<Integer> unionPure) {
        Set<List<Integer>> pentads = new HashSet<>();
        Set<List<Integer>> _clusters = new HashSet<>(clusters);

        if (unionPure.isEmpty()) {
            return new HashSet<>();
        }

        REMAINING:
        while (true) {
            if (remaining.size() < 5) break;

            if (verbose) {
                log("UnionPure = " + variablesForIndices(new ArrayList<>(unionPure)), false);
            }

            ChoiceGenerator gen = new ChoiceGenerator(remaining.size(), 5);
            int[] choice;

            while ((choice = gen.next()) != null) {
                int t2 = remaining.get(choice[0]);
                int t3 = remaining.get(choice[1]);
                int t4 = remaining.get(choice[2]);
                int t5 = remaining.get(choice[3]);
                int t6 = remaining.get(choice[4]);

                List<Integer> cluster = new ArrayList<>();
                cluster.add(t2);
                cluster.add(t3);
                cluster.add(t4);
                cluster.add(t5);
                cluster.add(t6);

                if (zeroCorr(cluster, 3)) {
                    continue;
                }

                // Check all x as a cross check; really only one should be necessary.
                boolean allvanish = true;
                boolean someVanish = false;

                for (int t1 : allVariables()) {
                    if (cluster.contains(t1)) continue;

                    List<Integer> _cluster = new ArrayList<>(cluster);
                    _cluster.add(t1);


                    if (vanishes(_cluster)) {
                        someVanish = true;
                    } else {
                        allvanish = false;
                        break;
                    }
                }

                if (someVanish && allvanish) {
                    pentads.add(cluster);
                    _clusters.add(cluster);
                    unionPure.addAll(cluster);
                    remaining.removeAll(cluster);

                    if (verbose) {
                        log("3-cluster found: " + variablesForIndices(cluster), false);
                    }

                    continue REMAINING;
                }
            }

            break;
        }

        return pentads;
    }

    private double significance(List<Integer> cluster) {
        double chisq = getClusterChiSquare(cluster);

        // From "Algebraic factor analysis: sextads, pentads and beyond" Drton et al.
        int n = cluster.size();
        int dof = dofHarman(n);
        double q = ProbUtils.chisqCdf(chisq, dof);
        return 1.0 - q;
    }

    private int dofDrton(int n) {
        int dof = ((n - 2) * (n - 3)) / 2 - 2;
        if (dof < 0) dof = 0;
        return dof;
    }

    private int dofHarman(int n) {
        int dof = n * (n - 5) / 2 + 1;
        if (dof < 0) dof = 0;
        return dof;
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

    private boolean pure(List<Integer> sextet) {
        if (vanishes(sextet)) {
            for (int o : allVariables()) {
                if (sextet.contains(o)) continue;

                for (int i = 0; i < sextet.size(); i++) {
                    List<Integer> _sextet = new ArrayList<>(sextet);
                    _sextet.remove(sextet.get(i));
                    _sextet.add(i, o);

                    if (!(vanishes(_sextet))) {
                        return false;
                    }
                }
            }

//            if (existsImpurePentad(sextet)) return false;

            System.out.println("PURE: " + variablesForIndices(sextet));

            return true;
        }

        return false;
    }

    private double getClusterChiSquare(List<Integer> cluster) {
        SemIm im = estimateClusterModel(cluster);
        return im.getChiSquare();
    }

    private SemIm estimateClusterModel(List<Integer> sextet) {
        Graph g = new EdgeListGraph();
        Node l1 = new GraphNode("L1");
        l1.setNodeType(NodeType.LATENT);
        Node l2 = new GraphNode("L2");
        l2.setNodeType(NodeType.LATENT);
        g.addNode(l1);
        g.addNode(l2);

        for (Integer aQuartet : sextet) {
            Node n = this.variables.get(aQuartet);
            g.addNode(n);
            g.addDirectedEdge(l1, n);
            g.addDirectedEdge(l2, n);
        }

        SemPm pm = new SemPm(g);

        SemEstimator est;

        if (dataModel instanceof DataSet) {
            est = new SemEstimator((DataSet) dataModel, pm, new SemOptimizerEm());
        } else {
            est = new SemEstimator((CovarianceMatrix) dataModel, pm, new SemOptimizerEm());
        }

        return est.estimate();
    }

    private SemIm estimateModel(List<List<Integer>> clusters) {
        Graph g = new EdgeListGraph();

        List<Node> upperLatents = new ArrayList<>();
        List<Node> lowerLatents = new ArrayList<>();

        for (int i = 0; i < clusters.size(); i++) {
            List<Integer> cluster = clusters.get(i);
            Node l1 = new GraphNode("L1." + (i + 1));
            l1.setNodeType(NodeType.LATENT);

            Node l2 = new GraphNode("L2." + (i + 1));
            l2.setNodeType(NodeType.LATENT);

            upperLatents.add(l1);
            lowerLatents.add(l2);

            g.addNode(l1);
            g.addNode(l2);

            for (Integer aCluster : cluster) {
                Node n = this.variables.get(aCluster);
                g.addNode(n);
                g.addDirectedEdge(l1, n);
                g.addDirectedEdge(l2, n);
            }
        }

        for (int i = 0; i < upperLatents.size(); i++) {
            for (int j = i + 1; j < upperLatents.size(); j++) {
                g.addDirectedEdge(upperLatents.get(i), upperLatents.get(j));
                g.addDirectedEdge(lowerLatents.get(i), lowerLatents.get(j));
            }
        }

        for (int i = 0; i < upperLatents.size(); i++) {
            for (int j = 0; j < lowerLatents.size(); j++) {
                if (i == j) continue;
                g.addDirectedEdge(upperLatents.get(i), lowerLatents.get(j));
            }
        }

        SemPm pm = new SemPm(g);

        for (Node node : upperLatents) {
            Parameter p = pm.getParameter(node, node);
            p.setFixed(true);
            p.setStartingValue(1.0);
        }

        for (Node node : lowerLatents) {
            Parameter p = pm.getParameter(node, node);
            p.setFixed(true);
            p.setStartingValue(1.0);
        }

        SemEstimator est;

        if (dataModel instanceof DataSet) {
            est = new SemEstimator((DataSet) dataModel, pm, new SemOptimizerEm());
        } else {
            est = new SemEstimator((CovarianceMatrix) dataModel, pm, new SemOptimizerEm());
        }

        return est.estimate();
    }

    private List<Integer> cluster(int...vars) {
        List<Integer> cluster = new ArrayList<>();

        for (int i : vars) {
            cluster.add(i);
        }

        if (new HashSet<>(cluster).size() < vars.length)
            throw new IllegalArgumentException("cluster elements must be unique: <" + cluster + ">");

        return cluster;
    }

    private boolean vanishes(List<Integer> sextet) {
//        if (zeroCorr(sextet, 4)) {
//            return false;
//        }

        Collections.sort(sextet);

//        PermutationGenerator gen = new PermutationGenerator(6);
//        int[] perm;

//        while ((perm = gen.next()) != null) {
//            int n1 = sextet.get(perm[0]);
//            int n2 = sextet.get(perm[1]);
//            int n3 = sextet.get(perm[2]);
//            int n4 = sextet.get(perm[3]);
//            int n5 = sextet.get(perm[4]);
//            int n6 = sextet.get(perm[5);
//
//            if (!vanishes(n1, n2, n3, n4, n5, n6)) return false;
//        }
//
//        return true;

        int n1 = sextet.get(0);
        int n2 = sextet.get(1);
        int n3 = sextet.get(2);
        int n4 = sextet.get(3);
        int n5 = sextet.get(4);
        int n6 = sextet.get(5);

        return vanishes(n1, n2, n3, n4, n5, n6);
//                && vanishes(n3, n2, n1, n6, n5, n4)
//                && vanishes(n4, n5, n6, n1, n2, n3)
//                && vanishes(n6, n5, n4, n3, n2, n1);
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

//        System.out.println("num zeroes = " + count);

        return count >= n;
//        return false;
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

    private boolean vanishes(int n1, int n2, int n3, int n4, int n5, int n6) {
        List<List<IntSextad>> independents = getIntSextads(n1, n2, n3, n4, n5, n6);

//        IntSextad[] all = {t1, t2, t3, t4, t5, t6, t7, t8, t9, t10};
//        independents.add(all);

        for (List<IntSextad> sextads : independents) {
            for (IntSextad sextad : sextads) {
                if (zeroCorr(sextad.getNodes(), 5)) continue;
            }

            double p = 0;
            try {
                p = test.getPValue(sextads);
            } catch (Exception e) {
                return true;
            }
//            System.out.println("p = " + p);
            if (p < alpha) return false;
        }


//        IntSextad[] sextads = new IntSextad[]{t1, t2, t3, t4, t5, t6, t7, t8, t9, t10};
//
//        for (IntSextad sextad : sextads) {
//            if (test.getPValue(sextad) < alpha) return false;
//        }

//        for (int i = 0; i < independents.size(); i++) {
//            if (test.getPValue(independents.get(0)) < alpha) {
//                return false;
//            }
//        }
//
//
//        return true;

//        for (List<IntSextad> sextads : independents) {
//            for (IntSextad sextad : sextads) {
//                if (logsum2(sextad.getNodes()) < -30) {
//                    return false;
//                }
//            }
//        }

        return true;
    }

    private List<List<IntSextad>> getIntSextads(int n1, int n2, int n3, int n4, int n5, int n6) {
        IntSextad t1 = new IntSextad(n1, n2, n3, n4, n5, n6);
        IntSextad t2 = new IntSextad(n1, n5, n6, n2, n3, n4);
        IntSextad t3 = new IntSextad(n1, n4, n6, n2, n3, n5);
        IntSextad t4 = new IntSextad(n1, n4, n5, n2, n3, n6);
        IntSextad t5 = new IntSextad(n1, n3, n4, n2, n5, n6);
        IntSextad t6 = new IntSextad(n1, n3, n5, n2, n4, n6);
        IntSextad t7 = new IntSextad(n1, n3, n6, n2, n4, n5);
        IntSextad t8 = new IntSextad(n1, n2, n4, n3, n5, n6);
        IntSextad t9 = new IntSextad(n1, n2, n5, n3, n4, n6);
        IntSextad t10 = new IntSextad(n1, n2, n6, n3, n4, n5);

        // The four sextads implied by equation 5.17 in Harmann.
        // independents.add(new IntSextad[]{t3, t7, t8, t9});

//            IntSextad[] independents = {t2, t5, t10, t3, t6};

        List<List<IntSextad>> independents = new ArrayList<>();
//        independents.add(sextadList(t1, t2, t3, t5, t6));
//        independents.add(sextadList(t1, t2, t3, t9, t10));
//        independents.add(sextadList(t6, t7, t8, t9, t10));
//        independents.add(sextadList(t1, t2, t4, t5, t9));
        independents.add(sextadList(t1, t3, t4, t6, t10));
        return independents;
    }

    private List<IntSextad> sextadList(IntSextad t1, IntSextad t2, IntSextad t3, IntSextad t5, IntSextad t6) {
        List<IntSextad> list = new ArrayList<>();
        list.add(t1);
        list.add(t2);
        list.add(t3);
        list.add(t5);
        list.add(t6);
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

    private Set<Integer> unionPure(Set<List<Integer>> pureClusters) {
        Set<Integer> unionPure = new HashSet<>();

        for (List<Integer> cluster : pureClusters) {
            unionPure.addAll(cluster);
        }

        return unionPure;
    }

    private void log(String s, boolean toLog) {
        if (toLog) {
            TetradLogger.getInstance().log("info", s);
            System.out.println(s);
        }
    }
}




