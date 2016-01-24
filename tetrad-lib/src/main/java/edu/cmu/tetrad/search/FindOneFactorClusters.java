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
import static java.lang.Math.sqrt;


/**
 * Implements FindOneFactorCluster by Erich Kummerfeld (adaptation of a two factor
 * quartet algorithm to a one factor tetrad algorithm).
 *
 * @author Joseph Ramsey
 */
public class FindOneFactorClusters {

    private final ICovarianceMatrix cov;

    public enum Algorithm {SAG, GAP}

    private CorrelationMatrix corr;
    // The list of all variables.
    private List<Node> variables;

    // The significance level.
    private double alpha;

    private TestType testType = TestType.TETRAD_DELTA;

    // The Delta test. Testing two tetrads simultaneously.
    private DeltaTetradTest test;

    // The tetrad test--using Ricardo's. Used only for Wishart.
    private ContinuousTetradTest test2;

    // The data.
    private transient DataModel dataModel;

    private List<List<Node>> clusters;

    private boolean verbose = false;
    private boolean significanceCalculated = false;
    private Algorithm algorithm = Algorithm.GAP;
    private Map<Set<Integer>, Double> avgSumLnPs = new HashMap<>();


    //========================================PUBLIC METHODS====================================//

    public FindOneFactorClusters(ICovarianceMatrix cov, TestType testType, Algorithm algorithm, double alpha) {
        if (testType == null) throw new NullPointerException("Null test type.");
        cov = new CovarianceMatrix(cov);
        this.variables = cov.getVariables();
        this.alpha = alpha;
        this.testType = testType;
        this.test = new DeltaTetradTest(cov);
        this.test2 = new ContinuousTetradTest(cov, testType, alpha);
        this.dataModel = cov;
        this.algorithm = algorithm;
        this.cov = cov;

        this.corr = new CorrelationMatrix(cov);


    }

    public FindOneFactorClusters(DataSet dataSet, TestType testType, Algorithm algorithm, double alpha) {
        if (testType == null) throw new NullPointerException("Null test type.");
        this.variables = dataSet.getVariables();
        this.alpha = alpha;
        this.testType = testType;
        this.test = new DeltaTetradTest(dataSet);
        this.test2 = new ContinuousTetradTest(dataSet, testType, alpha);
        this.dataModel = dataSet;
        this.algorithm = algorithm;

        this.corr = new CorrelationMatrix(dataSet);
        this.cov = new CovarianceMatrix(dataSet);
    }

    public Algorithm getAlgorithm() {
        return algorithm;
    }

    public void setAlgorithm(Algorithm algorithm) {
        this.algorithm = algorithm;
    }

    public Graph search() {
        Set<List<Integer>> allClusters;

        if (algorithm == Algorithm.SAG) {
            allClusters = estimateClustersSAG();
        } else if (algorithm == Algorithm.GAP) {
            allClusters = estimateClustersGAP();
        } else {
            throw new IllegalStateException("Expected SAG or GAP: " + testType);
        }
        this.clusters = variablesForIndices2(allClusters);
        return convertToGraph(allClusters);
    }

    public boolean isSignificanceCalculated() {
        return significanceCalculated;
    }

    public void setSignificanceCalculated(boolean significanceCalculated) {
        this.significanceCalculated = significanceCalculated;
    }

    //========================================PRIVATE METHODS====================================//


    // renjiey
    private int findFrequentestIndex(Integer outliers[]) {
        Map<Integer, Integer> map = new HashMap<>();

        for (int i = 0; i < outliers.length; i++) {
            if (map.containsKey(outliers[i])) {
                map.put(outliers[i], map.get(outliers[i]) + 1);
            } else {
                map.put(outliers[i], 1);
            }
        }

        Set<Map.Entry<Integer, Integer>> set = map.entrySet();
        Iterator<Map.Entry<Integer, Integer>> it = set.iterator();
        int nums = 0;// how many times variable occur
        int key = 0;// the number occur the most times

        while (it.hasNext()) {
            Map.Entry<Integer, Integer> entry = it.next();
            if (entry.getValue() > nums) {
                nums = entry.getValue();
                key = entry.getKey();
            }
        }

        return (key);
    }

    // This is the main function. It remove variables in the data such that the remaining correlation matrix
    // does not contain extreme value
    // Inputs: correlation matrix, upper and lower bound for unacceptable correlations
    // Output: and dynamic array of removed variables
    // renjiey
    private ArrayList<Integer> removeVariables(TetradMatrix correlationMatrix, double lowerBound, double upperBound,
                                               double percentBound) {
        Integer outlier[] = new Integer[correlationMatrix.rows() * (correlationMatrix.rows() - 1)];
        int count = 0;
        for (int i = 2; i < (correlationMatrix.rows() + 1); i++) {
            for (int j = 1; j < i; j++) {

                if ((Math.abs(correlationMatrix.get(i - 1, j - 1)) < lowerBound)
                        || (Math.abs(correlationMatrix.get(i - 1, j - 1)) > upperBound)) {
                    outlier[count * 2] = i;
                    outlier[count * 2 + 1] = j;

                } else {
                    outlier[count * 2] = 0;
                    outlier[count * 2 + 1] = 0;
                }
                count = count + 1;
            }
        }

        //find out the variables that should be deleted
        ArrayList<Integer> removedVariables = new ArrayList<>();

        // Added the percent bound jdramsey
        while (outlier.length > 1 && removedVariables.size() < percentBound * correlationMatrix.rows()) {
            //find out the variable that occurs most frequently in outlier
            int worstVariable = findFrequentestIndex(outlier);
            if (worstVariable > 0) {
                removedVariables.add(worstVariable);
            }

            //remove the correlations having the bad variable (change the relevant variables to 0)
            for (int i = 1; i < outlier.length + 1; i++) {
                if (outlier[i - 1] == worstVariable) {
                    outlier[i - 1] = 0;

                    if (i % 2 != 0) {
                        outlier[i] = 0;
                    } else {
                        outlier[i - 2] = 0;
                    }
                }
            }

            //delete zero elements in outlier
            outlier = removeZeroIndex(outlier);
        }

        log(removedVariables.size() + " variables removed: " + variablesForIndices(removedVariables), true);

        return (removedVariables);
    }

    // renjiey
    private Integer[] removeZeroIndex(Integer outlier[]) {
        List<Integer> list = new ArrayList<>();
        for (int i = 0; i < outlier.length; i++) {
            list.add(outlier[i]);
        }
        for (Integer element : outlier) {
            if (element < 1) {
                list.remove(element);
            }
        }
        return list.toArray(new Integer[1]);
    }

    // This is the main algorithm.
    private Set<List<Integer>> estimateClustersGAP() {
//        List<Integer> _variables = new ArrayList<Integer>();
//        for (int i = 0; i < variables.size(); i++) _variables.add(i);
        List<Integer> _variables = allVariables();

        Set<Set<Integer>> triples = findPuretriplesGAP(_variables);
        Set<Set<Integer>> combined = combinePuretriplesGAP(triples, _variables);

        Set<List<Integer>> _combined = new HashSet<>();

        for (Set<Integer> c : combined) {
            List a = new ArrayList<>(c);
            Collections.sort(a);
            _combined.add(a);
        }

        return _combined;

    }

    private List<Integer> allVariables() {
        List<Integer> _variables = new ArrayList<>();
        for (int i = 0; i < variables.size(); i++) _variables.add(i);
        return _variables;
    }

    private Set<Set<Integer>> findPuretriplesGAP(List<Integer> allVariables) {
        if (allVariables.size() < 4) {
            return new HashSet<>();
        }

        log("Finding pure triples.", true);

        ChoiceGenerator gen = new ChoiceGenerator(allVariables.size(), 3);
        int[] choice;
        Set<Set<Integer>> puretriples = new HashSet<>();

        while ((choice = gen.next()) != null) {
            int n1 = allVariables.get(choice[0]);
            int n2 = allVariables.get(choice[1]);
            int n3 = allVariables.get(choice[2]);

            List<Integer> triple = triple(n1, n2, n3);

            if (zeroCorr(triple)) continue;

            if (!pureTriple(triple)) continue;

            if (verbose) {
                log("++" + variablesForIndices(triple), false);
            }

            puretriples.add(new HashSet<>(triple));
        }

        return puretriples;
    }

    private Set<List<Integer>> estimateClustersSAG() {
        List<Integer> _variables = new ArrayList<>(allVariables());
        Set<List<Integer>> clusters = new HashSet<>();

        VARIABLES:
        while (!_variables.isEmpty()) {
            if (_variables.size() < 3) break;

            ChoiceGenerator gen = new ChoiceGenerator(_variables.size(), 3);
            int[] choice;

            while ((choice = gen.next()) != null) {
                int n1 = _variables.get(choice[0]);
                int n2 = _variables.get(choice[1]);
                int n3 = _variables.get(choice[2]);

                List<Integer> cluster = triple(n1, n2, n3);

                if (zeroCorr(cluster)) continue;

                // Note that purity needs to be assessed with respect to all of the variables in order to
                // remove all latent-measure impurities between pairs of latents.
                if (pureTriple(cluster)) {
                    if (verbose) {
                        log("Found a pure: " + variablesForIndices(cluster), true);
                    }

                    List<Integer> nodes = new ArrayList<>(cluster);

                    addOtherVariablesStrict(_variables, nodes);
                    nodes = new ArrayList<>(nodes);
                    Set<Integer> _nodes = pickSignificantSubcluster(new HashSet<>(nodes));
                    if (_nodes != null) nodes = new ArrayList<>(_nodes);

                    if (nodes.size() < 5) continue;

                    if (verbose) {
                        log("Cluster found: " + variablesForIndices(cluster), true);
                    }

                    clusters.add(nodes);
                    _variables.removeAll(nodes);

                    continue VARIABLES;
                }
            }

            break;
        }

        return clusters;
    }

    private void addOtherVariablesStrict(List<Integer> _variables, List<Integer> cluster) {

        O:
        for (int o : _variables) {
            if (cluster.contains(o)) continue;
            List<Integer> _cluster = new ArrayList<>(cluster);

            ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 2);
            int[] choice;

            while ((choice = gen2.next()) != null) {
                int t1 = _cluster.get(choice[0]);
                int t2 = _cluster.get(choice[1]);

                List<Integer> triple = new ArrayList<>();
                triple.add(t1);
                triple.add(t2);
                triple.add(o);

                if (!pureTriple(triple)) {
                    continue O;
                }
            }

            Integer o1 = o;
            cluster.add(o1);

//            if (getPValue(cluster) < alpha) {
//                cluster.remove(o1);
//            }

            log("Extending by " + variables.get(o1), true);
        }
    }

    private boolean pureTriple(List<Integer> triple) {
        for (int o : allVariables()) {
            if (triple.contains(o)) {
                continue;
            }

            List<Integer> quartet = new ArrayList<>(triple);
            quartet.add(o);

            if (!vanishes(quartet)) {
                return false;
            }
        }

        return true;
    }

    private Set<Set<Integer>> combinePuretriplesGAP(Set<Set<Integer>> puretriples, List<Integer> _variables) {
        log("Growing pure triples.", true);
        Set<Set<Integer>> grown = new HashSet<>();

        // Lax grow phase with speedup.
        if (true) {
            Set<Integer> t = new HashSet<>();
            int count = 0;
            int total = puretriples.size();

            do {
                if (!puretriples.iterator().hasNext()) {
                    break;
                }

                Set<Integer> cluster = puretriples.iterator().next();
                Set<Integer> _cluster = new HashSet<>(cluster);

                for (int o : _variables) {
                    if (_cluster.contains(o)) continue;

                    List<Integer> _cluster2 = new ArrayList<>(_cluster);
                    int possible = MathUtils.choose(_cluster2.size(), 2);
                    int accepted = 0;

                    ChoiceGenerator gen = new ChoiceGenerator(_cluster2.size(), 2);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        t.clear();
                        t.add(_cluster2.get(choice[0]));
                        t.add(_cluster2.get(choice[1]));
                        t.add(o);

                        if (!puretriples.contains(t)) {
                            continue;
                        }

                        accepted++;
                    }

                    if (accepted < .8 * possible) {
                        continue;
                    }

//                    if (_cluster.size() > 3 && getPValue(new ArrayList<>(_cluster)) > alpha) {
//                        _cluster.add(o);
//                    } else if (_cluster.size() == 3) {
//                        _cluster.add(o);
//                    }

                    _cluster.add(o);

//                    if (!(avgSumLnP(new ArrayList<Integer>(_cluster)) > -10)) {
//                        _cluster.remove(o);
//                    }
                }

//                Set<Integer> bestCluster = pickBestSubcluster(_cluster);
                Set<Integer> bestCluster = pickSignificantSubcluster(_cluster);
//                Set<Integer> bestCluster = _cluster;

                // This takes out all pure clusters that are subsets of _cluster.
                ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 3);
                int[] choice2;
                List<Integer> _cluster3 = new ArrayList<>(_cluster);

                while ((choice2 = gen2.next()) != null) {
                    int n1 = _cluster3.get(choice2[0]);
                    int n2 = _cluster3.get(choice2[1]);
                    int n3 = _cluster3.get(choice2[2]);

                    t.clear();
                    t.add(n1);
                    t.add(n2);
                    t.add(n3);

                    puretriples.remove(t);
                }

                if (verbose) {
                    System.out.println("Grown " + (++count) + " of " + total + ": " + variablesForIndices(new ArrayList<>(_cluster)));
                }

                if (bestCluster != null && bestCluster.size() >= 5) {
                    grown.add(bestCluster);
                }
            } while (!puretriples.isEmpty());
        }

        // Lax grow phase without speedup.
        if (false) {
            int count = 0;
            int total = puretriples.size();

            // Optimized lax version of grow phase.
            for (Set<Integer> cluster : new HashSet<>(puretriples)) {
                Set<Integer> _cluster = new HashSet<>(cluster);

                for (int o : _variables) {
                    if (_cluster.contains(o)) continue;

                    List<Integer> _cluster2 = new ArrayList<>(_cluster);
                    int rejected = 0;
                    int accepted = 0;

                    ChoiceGenerator gen = new ChoiceGenerator(_cluster2.size(), 4);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        int n1 = _cluster2.get(choice[0]);
                        int n2 = _cluster2.get(choice[1]);

                        List<Integer> triple = triple(n1, n2, o);

                        Set<Integer> t = new HashSet<>(triple);

                        if (!puretriples.contains(t)) {
                            rejected++;
                        } else {
                            accepted++;
                        }

//                        if (avgSumLnP(triple) < -10) continue CLUSTER;
                    }

                    if (rejected > accepted) {
                        continue;
                    }

                    _cluster.add(o);
                }

                for (Set<Integer> c : new HashSet<>(puretriples)) {
                    if (_cluster.containsAll(c)) {
                        puretriples.remove(c);
                    }
                }

                if (verbose) {
                    System.out.println("Grown " + (++count) + " of " + total + ": " + _cluster);
                }

                grown.add(_cluster);
            }
        }

        // Strict grow phase.
        if (false) {
            Set<Integer> t = new HashSet<>();
            int count = 0;
            int total = puretriples.size();

            do {
                if (!puretriples.iterator().hasNext()) {
                    break;
                }

                Set<Integer> cluster = puretriples.iterator().next();
                Set<Integer> _cluster = new HashSet<>(cluster);

                VARIABLES:
                for (int o : _variables) {
                    if (_cluster.contains(o)) continue;

                    List<Integer> _cluster2 = new ArrayList<>(_cluster);

                    ChoiceGenerator gen = new ChoiceGenerator(_cluster2.size(), 2);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        int n1 = _cluster2.get(choice[0]);
                        int n2 = _cluster2.get(choice[1]);

                        t.clear();
                        t.add(n1);
                        t.add(n2);
                        t.add(o);

                        if (!puretriples.contains(t)) {
                            continue VARIABLES;
                        }

//                        if (avgSumLnP(new ArrayList<Integer>(t)) < -10) continue CLUSTER;
                    }

                    _cluster.add(o);
                }

                // This takes out all pure clusters that are subsets of _cluster.
                ChoiceGenerator gen2 = new ChoiceGenerator(_cluster.size(), 3);
                int[] choice2;
                List<Integer> _cluster3 = new ArrayList<>(_cluster);

                while ((choice2 = gen2.next()) != null) {
                    int n1 = _cluster3.get(choice2[0]);
                    int n2 = _cluster3.get(choice2[1]);
                    int n3 = _cluster3.get(choice2[2]);

                    t.clear();
                    t.add(n1);
                    t.add(n2);
                    t.add(n3);

                    puretriples.remove(t);
                }

                if (verbose) {
                    System.out.println("Grown " + (++count) + " of " + total + ": " + _cluster);
                }
                grown.add(_cluster);
            } while (!puretriples.isEmpty());
        }

        if (false) {
            System.out.println("# pure triples = " + puretriples.size());

            List<Set<Integer>> clusters = new LinkedList<>(puretriples);
            Set<Integer> t = new HashSet<>();

            I:
            for (int i = 0; i < clusters.size(); i++) {
                System.out.println("I = " + i);

//                // remove "i" clusters that intersect with previous clusters.
//                for (int k = 0; k < i - 1; k++) {
//                    Set<Integer> ck = clusters.get(k);
//                    Set<Integer> ci = clusters.get(i);
//
//                    if (ck == null) continue;
//                    if (ci == null) continue;
//
//                    Set<Integer> cm = new HashSet<Integer>(ck);
//                    cm.retainAll(ci);
//
//                    if (!cm.isEmpty()) {
//                        clusters.remove(i);
//                        i--;
//                        continue I;
//                    }
//                }

                J:
                for (int j = i + 1; j < clusters.size(); j++) {
                    Set<Integer> ci = clusters.get(i);
                    Set<Integer> cj = clusters.get(j);

                    if (ci == null) continue;
                    if (cj == null) continue;

                    Set<Integer> ck = new HashSet<>(ci);
                    ck.addAll(cj);

                    List<Integer> cm = new ArrayList<>(ck);

                    ChoiceGenerator gen = new ChoiceGenerator(cm.size(), 3);
                    int[] choice;

                    while ((choice = gen.next()) != null) {
                        t.clear();
                        t.add(cm.get(choice[0]));
                        t.add(cm.get(choice[1]));
                        t.add(cm.get(choice[2]));

                        if (!puretriples.contains(t)) {
                            continue J;
                        }
                    }

                    clusters.set(i, ck);
                    clusters.remove(j);
                    j--;
                    System.out.println("Removing " + ci + ", " + cj + ", adding " + ck);
                }
            }

            grown = new HashSet<>(clusters);
        }

        // Optimized pick phase.
        log("Choosing among grown clusters.", true);

//        for (Set<Integer> l : grown) {
//            ArrayList<Integer> _l = new ArrayList<>(l);
//            Collections.sort(_l);
//            if (verbose) {
//                log("Grown: " + variablesForIndices(_l), false);
//            }
//        }

        Set<Set<Integer>> out = new HashSet<>();

        List<Set<Integer>> list = new ArrayList<>(grown);

        Collections.sort(list, new Comparator<Set<Integer>>() {
            @Override
            public int compare(Set<Integer> o1, Set<Integer> o2) {
                return o2.size() - o1.size();
            }
        });

//        final Map<Set<Integer>, Double> significances = new HashMap<Set<Integer>, Double>();

//        Collections.sort(list, new Comparator<Set<Integer>>() {
//            @Override
//            public int compare(Set<Integer> cluster1, Set<Integer> cluster2) {
////                Double sum1 = significances.get(cluster1);
////                if (sum1 == null) {
////                    double sig = significance(new ArrayList<Integer>(cluster1));
////                    significances.put(cluster1, sig);
////                    sum1 = sig;
////                }
////                Double sum2 = significances.get(cluster2);
////                if (sum2 == null) {
////                    double sig = significance(new ArrayList<Integer>(cluster2));
////                    significances.put(cluster2, sig);
////                    sum2 = sig;
////                }
//
//                double avg1 = avgSumLnP(new ArrayList<Integer>(cluster1));
//                double avg2 = avgSumLnP(new ArrayList<Integer>(cluster2));
//
//                return Double.compare(avg2, avg1);
//            }
//        });

        Set<Integer> all = new HashSet<>();

        CLUSTER:
        for (Set<Integer> cluster : list) {
            for (Integer i : cluster) {
                if (all.contains(i)) continue CLUSTER;
            }

            if (cluster.size() > 3) {
                out.add(cluster);
                all.addAll(cluster);
            }

            System.out.println("P value for cluster " + variablesForIndices(new ArrayList<>(cluster)) + " = " + getPValue(new ArrayList<>(cluster)));
        }

        for (Set<Integer> _out : out) {
            log("OUT: " + variablesForIndices(new ArrayList<>(_out)), true);
        }

        return out;
    }

    private Set<Integer> pickBestSubcluster(Set<Integer> _cluster) {
        DepthChoiceGenerator gen = new DepthChoiceGenerator(_cluster.size(), Math.min(10, _cluster.size()));
        int[] choice;
        List<Integer> __cluster = new ArrayList<>(_cluster);
        double maxP = 0.0;
        List<Node> maxCluster = null;
        int minSize = 4;

        while ((choice = gen.next()) != null) {
            if (choice.length < minSize) continue;
            List<Node> nodes = new ArrayList<>();
            for (int aChoice : choice) nodes.add(variables.get(__cluster.get(aChoice)));
            Mimbuild2 mimbuild = new Mimbuild2();
            mimbuild.setMethod(Mimbuild2.Method.STRAIGHT);
            mimbuild.search(Collections.singletonList(nodes), Collections.singletonList("L"), cov);
            double p = mimbuild.getPValue();
            if (p > alpha && p > maxP) {
                maxP = p;
                maxCluster = nodes;
            }
        }

        if (maxCluster == null) return null;

        Set<Integer> cluster = new HashSet<>();
        for (Node node : maxCluster) cluster.add(variables.indexOf(node));

        return cluster;
    }

    private Set<Integer> pickSignificantSubcluster(Set<Integer> _cluster) {
        List<Integer> __cluster = new ArrayList<>(_cluster);
        int minSize = 4;

        for (int i = _cluster.size(); i >= minSize; i--) {
            ChoiceGenerator gen = new ChoiceGenerator(__cluster.size(), i);
            int[] choice;

            while ((choice = gen.next()) != null) {
                List<Node> nodes = new ArrayList<>();
                for (int aChoice : choice) nodes.add(variables.get(__cluster.get(aChoice)));
                Mimbuild2 mimbuild = new Mimbuild2();
                mimbuild.setMethod(Mimbuild2.Method.STRAIGHT);
                mimbuild.search(Collections.singletonList(nodes), Collections.singletonList("L"), cov);
                double p = mimbuild.getPValue();
                if (p > alpha) {
                    Set<Integer> indices = new HashSet<>();
                    for (int aChoice : choice) indices.add(__cluster.get(aChoice));
                    return indices;
                }
            }
        }

        return null;
    }

//    private int dofDrton(int n) {
//        int dof = ((n - 2) * (n - 3)) / 2 - 2;
//        if (dof < 0) dof = 0;
//        return dof;
//    }
//
//    private int dofHarman(int n) {
//        int dof = n * (n - 5) / 2 + 1;
//        if (dof < 0) dof = 0;
//        return dof;
//    }

    private List<Node> variablesForIndices(List<Integer> cluster) {
        List<Node> _cluster = new ArrayList<>();

        for (int c : cluster) {
            _cluster.add(variables.get(c));
        }

//        Collections.sort(_cluster);

        return _cluster;
    }

    private List<List<Node>> variablesForIndices2(Set<List<Integer>> clusters) {
        List<List<Node>> variables = new ArrayList<>();

        for (List<Integer> cluster : clusters) {
            variables.add(variablesForIndices(cluster));
        }

        return variables;
    }

    private boolean pure(List<Integer> quartet) {
        if (zeroCorr(quartet)) {
            return false;
        }

        if (vanishes(quartet)) {
            for (int o : allVariables()) {
                if (quartet.contains(o)) continue;

                for (int i = 0; i < quartet.size(); i++) {
                    List<Integer> _quartet = new ArrayList<>(quartet);
                    _quartet.remove(quartet.get(i));
                    _quartet.add(o);

//                    if (zeroCorr(_quartet)) {
//                        continue;
//                    }

                    if (!(vanishes(_quartet))) {
                        return false;
                    }
                }
            }

            return true;
        }

        return false;
    }

    public double getPValue(List<Integer> cluster) {
        Mimbuild2 mimbuild = new Mimbuild2();
        List<List<Node>> clustering = new ArrayList<>();
        clustering.add(variablesForIndices(cluster));
        mimbuild.search(clustering, Collections.singletonList("L"), cov);
        return mimbuild.getPValue();
    }

    private SemIm estimateClusterModel(List<Integer> quartet) {
        Graph g = new EdgeListGraph();
        Node l1 = new GraphNode("L1");
        l1.setNodeType(NodeType.LATENT);
        Node l2 = new GraphNode("L2");
        l2.setNodeType(NodeType.LATENT);
        g.addNode(l1);
        g.addNode(l2);

        for (int i = 0; i < quartet.size(); i++) {
            Node n = this.variables.get(quartet.get(i));
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

    private List<Integer> quartet(int n1, int n2, int n3, int n4) {
        List<Integer> quartet = new ArrayList<>();
        quartet.add(n1);
        quartet.add(n2);
        quartet.add(n3);
        quartet.add(n4);

        if (new HashSet<>(quartet).size() < 4)
            throw new IllegalArgumentException("quartet elements must be unique: <" + n1 + ", " + n2 + ", " + n3 + ", " + n4 + ">");

        return quartet;
    }

    private List<Integer> triple(int n1, int n2, int n3) {
        List<Integer> triple = new ArrayList<>();
        triple.add(n1);
        triple.add(n2);
        triple.add(n3);

        if (new HashSet<>(triple).size() < 3)
            throw new IllegalArgumentException("triple elements must be unique: <" + n1 + ", " + n2 + ", " + n3 + ">");

        return triple;
    }

    private boolean vanishes(List<Integer> quartet) {
        int n1 = quartet.get(0);
        int n2 = quartet.get(1);
        int n3 = quartet.get(2);
        int n4 = quartet.get(3);

        return vanishes(n1, n2, n3, n4);
    }

    private boolean zeroCorr(List<Integer> cluster) {
        int count = 0;

        for (int i = 0; i < cluster.size(); i++) {
            for (int j = i + 1; j < cluster.size(); j++) {
                double r = this.corr.getValue(cluster.get(i), cluster.get(j));
                int N = this.corr.getSampleSize();
                double f = sqrt(N) * Math.log((1. + r) / (1. - r));
                double p = 2.0 * (1.0 - RandomUtil.getInstance().normalCdf(0, 1, abs(f)));
                if (p > alpha) count++;
            }
        }

        return count >= 1;
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

    private boolean vanishes(int x, int y, int z, int w) {
        if (testType == TestType.TETRAD_DELTA) {
            Tetrad t1 = new Tetrad(variables.get(x), variables.get(y), variables.get(z), variables.get(w));
            Tetrad t2 = new Tetrad(variables.get(x), variables.get(y), variables.get(w), variables.get(z));

            return test.getPValue(t1, t2) > alpha;
        } else if (testType == TestType.TETRAD_WISHART) {
            return test2.tetradPValue(x, y, z, w) > alpha && test2.tetradPValue(x, y, w, z) > alpha;
        }

        throw new IllegalArgumentException("Only the delta and wishart tests are being used: " + testType);
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
        }

        System.out.println(s);
    }
}




