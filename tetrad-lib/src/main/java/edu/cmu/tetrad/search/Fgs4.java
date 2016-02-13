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
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.TetradLogger;

import java.io.PrintStream;
import java.util.*;
import java.util.concurrent.*;


/**
 * GesSearch is an implementation of the GES algorithm, as specified in Chickering (2002) "Optimal structure
 * identification with greedy search" Journal of Machine Learning Research. It works for both BayesNets and SEMs.
 * <p>
 * Some code optimization could be done for the scoring part of the graph for discrete models (method scoreGraphChange).
 * Some of Andrew Moore's approaches for caching sufficient statistics, for instance.
 * <p>
 * To speed things up, it has been assumed that variables X and Y with zero correlation do not correspond to edges in
 * the graph. This is a restricted form of the faithfulness assumption, something GES does not assume. This
 * faithfulness assumption needs to be explicitly turned on using setFaithfulnessAssumed(true).
 * <p>
 * A number of other optimizations were added 5/2015. See code for details.
 *
 * @author Ricardo Silva, Summer 2003
 * @author Joseph Ramsey, Revisions 5/2015
 */
public final class Fgs4 implements GraphSearch, GraphScorer {

    /**
     * Sample size, either from the data set or from the covariances.
     */
    private int sampleSize;

    /**
     * Specification of forbidden and required edges.
     */
    private IKnowledge knowledge = new Knowledge2();

    /**
     * List of variables in the data set, in order.
     */
    private List<Node> variables;

    /**
     * True iff the data set is discrete.
     */
    private boolean discrete;

    /**
     * The true graph, if known. If this is provided, asterisks will be printed out next to false positive added edges
     * (that is, edges added that aren't adjacencies in the true graph).
     */
    private Graph trueGraph;

    /**
     * An initial graph to start from.
     */
    private Graph initialGraph;

    /**
     * If non-null, edges not adjacent in this graph will not be added.
     */
    private Graph boundGraph = null;

    /**
     * Elapsed time of the most recent search.
     */
    private long elapsedTime;

    /**
     * The depth of search for the forward reevaluation step.
     */
    private int depth = -1;

    /**
     * A bound on cycle length.
     */
    private int cycleBound = -1;

    /**
     * The score for discrete searches.
     */
    private GesScore gesScore;

    /**
     * The logger for this class. The config needs to be set.
     */
    private TetradLogger logger = TetradLogger.getInstance();

    private LinkedList<ScoredGraph> topGraphs = new LinkedList<>();

    /**
     * The number of top patterns to store.
     */
    private int numPatternsToStore = 0;

    /**
     * True if logs should be output.
     */
    private boolean log = true;

    /**
     * True if verbose output should be printed.
     */
    private boolean verbose = false;

    // Potential arrows sorted by bump high to low. The first one is a candidate for adding to the graph.
    private SortedSet<Arrow> sortedArrows = new ConcurrentSkipListSet<>();

    // Arrows added to sortedArrows for each <i, j>.
    private Map<OrderedPair<Node>, Set<Arrow>> lookupArrows;

    // Map from variables to their column indices in the data set.
    private ConcurrentMap<Node, Integer> hashIndices;

    // The static ForkJoinPool instance.
    private ForkJoinPool pool = ForkJoinPoolInstance.getInstance().getPool();

    // A running tally of the total BIC score.
    private double score;

    // A graph where X--Y means that X and Y have non-zero total effect on one another.
    private Graph effectEdgesGraph;

    // The minimum number of operations to do before parallelizing.
    private final int minChunk = 100;

    // Where printed output is sent.
    private PrintStream out = System.out;

    // A initial adjacencies graph.
    private Graph adjacencies = null;

    // True if it is assumed that zero effect adjacencies are not in the graph.
    private boolean faithfulnessAssumed = true;

    // A utility map to help with orientation.
    private WeakHashMap<Node, Set<Node>> neighbors = new WeakHashMap<>();

    //===========================CONSTRUCTORS=============================//

    /**
     * The data set must either be all continuous or all discrete.
     */
    public Fgs4(DataSet dataSet) {
        out.println("GES constructor");

        if (dataSet.isDiscrete()) {
            setGesScore(new BDeuScore(dataSet));
        } else {
            setGesScore(new SemBicScore(new CovarianceMatrixOnTheFly(dataSet)));
        }

        out.println("GES constructor done");
    }

    /**
     * Continuous case--where a covariance matrix is already available.
     */
    public Fgs4(ICovarianceMatrix covMatrix) {
        out.println("GES constructor");

        setGesScore(new SemBicScore(covMatrix));

        out.println("GES constructor done");
    }


    public Fgs4(GesScore gesScore) {
        if (gesScore == null) throw new NullPointerException();
        setGesScore(gesScore);
    }

    //==========================PUBLIC METHODS==========================//

    /**
     * Set to true if it is assumed that all path pairs with one length 1 path do not cancel.
     */
    public void setFaithfulnessAssumed(boolean faithfulness) {
        this.faithfulnessAssumed = faithfulness;
    }

    /**
     * Returns true if it is assumed that all path pairs with one length 1 path do not cancel.
     */
    public boolean isFaithfulnessAssumed() {
        return this.faithfulnessAssumed;
    }

    /**
     * Greedy equivalence search: Start from the empty graph, add edges till model is significant. Then start deleting
     * edges till a minimum is achieved.
     *
     * @return the resulting Pattern.
     */
    public Graph search() {
        lookupArrows = new ConcurrentHashMap<>();
        final List<Node> nodes = new ArrayList<>(variables);
        this.effectEdgesGraph = getEffectEdges(nodes);

        if (adjacencies != null) {
            adjacencies = GraphUtils.replaceNodes(adjacencies, nodes);
        }

        Graph graph;

        if (initialGraph == null) {
            graph = new EdgeListGraphSingleConnections(getVariables());
        } else {
            graph = new EdgeListGraphSingleConnections(initialGraph);

            for (Edge edge : initialGraph.getEdges()) {
                if (!effectEdgesGraph.isAdjacentTo(edge.getNode1(), edge.getNode2())) {
                    effectEdgesGraph.addUndirectedEdge(edge.getNode1(), edge.getNode2());
                }
            }
        }

        addRequiredEdges(graph);

        topGraphs.clear();

        storeGraph(graph);

        long start = System.currentTimeMillis();
        score = 0.0;

        // Do forward search.
        fes(graph);

        // Do backward search.
        bes(graph);

        long endTime = System.currentTimeMillis();
        this.elapsedTime = endTime - start;
        this.logger.log("graph", "\nReturning this graph: " + graph);

        this.logger.log("info", "Elapsed time = " + (elapsedTime) / 1000. + " s");
        this.logger.flush();

        return graph;

    }

    /**
     * Returns the background knowledge.
     */
    public IKnowledge getKnowledge() {
        return knowledge;
    }

    /**
     * Sets the background knowledge.
     *
     * @param knowledge the knowledge object, specifying forbidden and required edges.
     */
    public void setKnowledge(IKnowledge knowledge) {
        if (knowledge == null) throw new NullPointerException();
        this.knowledge = knowledge;
    }

    public long getElapsedTime() {
        return elapsedTime;
    }

    /**
     * For BIC score, a multiplier on the penalty term. For continuous searches.
     */
    public double getPenaltyDiscount() {
        if (gesScore instanceof SemBicScore) {
            return ((SemBicScore) gesScore).getPenaltyDiscount();
        }

        return Double.NaN;
    }

    /**
     * For BIC score, a multiplier on the penalty term. For continuous searches.
     */
    public void setPenaltyDiscount(double penaltyDiscount) {
        if (penaltyDiscount < 0) {
            throw new IllegalArgumentException("Penalty discount must be >= 0: "
                    + penaltyDiscount);
        }

        if (gesScore instanceof SemBicScore) {
            ((SemBicScore) gesScore).setPenaltyDiscount(penaltyDiscount);
        }
    }

    /**
     * If the true graph is set, askterisks will be printed in log output for the true edges.
     */
    public void setTrueGraph(Graph trueGraph) {
        this.trueGraph = trueGraph;
    }

    /**
     * Returns the score of the given DAG, up to a constant.
     */
    public double getScore(Graph dag) {
        return scoreDag(dag);
    }

    /**
     * Returns the number of patterns to store.
     */
    public int getNumPatternsToStore() {
        return numPatternsToStore;
    }

    /**
     * Sets the number of patterns to store. This should be set to zero for fast search.
     */
    public void setNumPatternsToStore(int numPatternsToStore) {
        if (numPatternsToStore < 0) {
            throw new IllegalArgumentException("# graphs to store must at least 0: " + numPatternsToStore);
        }

        this.numPatternsToStore = numPatternsToStore;
    }

    /**
     * True iff log output should be produced.
     */
    public boolean isLog() {
        return log;
    }

    /**
     * Sets whether log output should be produced. Set to false a faster search.
     */
    public void setLog(boolean log) {
        this.log = log;
    }

    /**
     * Returns the initial graph for the search. The search is initialized to this graph and
     * proceeds from there.
     */
    public Graph getInitialGraph() {
        return initialGraph;
    }

    /**
     * Sets the initial graph.
     */
    public void setInitialGraph(Graph initialGraph) {
        if (initialGraph != null) {
            initialGraph = GraphUtils.replaceNodes(initialGraph, variables);

            out.println("Initial graph variables: " + initialGraph.getNodes());
            out.println("Data set variables: " + variables);

            if (!new HashSet<>(initialGraph.getNodes()).equals(new HashSet<>(variables))) {
                throw new IllegalArgumentException("Variables aren't the same.");
            }
        }

        this.initialGraph = initialGraph;
    }

    /**
     * Sets whether verbose output should be produced.
     */
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    /**
     * Sets the output stream that output (except for log output) should be sent to.
     * By detault System.out.
     */
    public void setOut(PrintStream out) {
        this.out = out;
    }

    /**
     * Returns the output stream that output (except for log output) should be sent to.
     */
    public PrintStream getOut() {
        return out;
    }

    /**
     * Returns the set of preset adjacenies for the algorithm; edges not in this adjacencies graph
     * will not be added.
     */
    public Graph getAdjacencies() {
        return adjacencies;
    }

    /**
     * Sets the set of preset adjacenies for the algorithm; edges not in this adjacencies graph
     * will not be added.
     */
    public void setAdjacencies(Graph adjacencies) {
        this.adjacencies = adjacencies;
    }

    /**
     * Returns the depth for the forward reevaluation step.
     */
    public int getDepth() {
        return depth;
    }

    /**
     * -1 for unlimited depth, otherwise a number >= 0. In the forward reevaluation step, subsets of neighbors up to
     * depth in size are considered. Limiting depth can speed up the algorithm.
     */
    public void setDepth(int depth) {
        this.depth = depth;
    }


    /**
     * A bound on cycle length.
     */
    public int getCycleBound() {
        return cycleBound;
    }

    /**
     * A bound on cycle length.
     *
     * @param cycleBound The bound, >= 1, or -1 for unlimited.
     */
    public void setCycleBound(int cycleBound) {
        if (!(cycleBound == -1 || cycleBound >= 1))
            throw new IllegalArgumentException("Cycle bound needs to be -1 or >= 1: " + cycleBound);
        this.cycleBound = cycleBound;
    }

    /**
     * Creates a new processors pool with the specified number of threads
     */
    public void setNumProcessors(int numProcessors) {
        this.pool = new ForkJoinPool(numProcessors);
    }

    /**
     * True iff edges that cause linear dependence are ignored.
     */
    public boolean isIgnoreLinearDependent() {
        if (gesScore instanceof SemBicScore) {
            return ((SemBicScore) gesScore).isIgnoreLinearDependent();
        }

        throw new UnsupportedOperationException("Operation supported only for SemBicScore.");
    }

    public void setIgnoreLinearDependent(boolean ignoreLinearDependent) {
        if (gesScore instanceof SemBicScore) {
            ((SemBicScore) gesScore).setIgnoreLinearDependent(ignoreLinearDependent);
        } else {
            throw new UnsupportedOperationException("Operation supported only for SemBicScore.");
        }
    }

    /**
     * If non-null, edges not adjacent in this graph will not be added.
     */
    public void setBoundGraph(Graph boundGraph) {
        this.boundGraph = GraphUtils.replaceNodes(boundGraph, getVariables());
    }

    //===========================PRIVATE METHODS========================//

    //Sets the discrete scoring function to use.
    private void setGesScore(GesScore gesScore) {
        this.gesScore = gesScore;
        this.variables = gesScore.getVariables();
        this.sampleSize = gesScore.getSampleSize();
        this.discrete = gesScore.isDiscrete();
    }


    // Simultaneously finds the first edge to add to an empty graph and finds all length 1 paths that are
    // not canceled by other paths (the "effect edges")
    private Graph getEffectEdges(final List<Node> nodes) {
        long start = System.currentTimeMillis();
        final Graph effectEdgesGraph = new EdgeListGraphSingleConnections(nodes);
        final Set<Node> emptySet = new HashSet<>(0);

        final int[] count = new int[1];

        class EffectTask extends RecursiveTask<Boolean> {
            private int chunk;
            private int from;
            private int to;

            public EffectTask(int chunk, int from, int to) {
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            @Override
            protected Boolean compute() {
                if (to - from <= chunk) {
                    for (int i = from; i < to; i++) {
                        synchronized (count) {
                            if (((count[0]++) + 1) % 1000 == 0) out.println("Initializing effect edges: " + count[0]);
                        }

                        Node y = nodes.get(i);

                        for (int j = i + 1; j < nodes.size(); j++) {
                            Node x = nodes.get(j);
//
                            if (existsKnowledge()) {
                                if (getKnowledge().isForbidden(x.getName(), y.getName()) && getKnowledge().isForbidden(y.getName(), x.getName())) {
                                    continue;
                                }

                                if (!validSetByKnowledge(y, emptySet)) {
                                    continue;
                                }
                            }

                            if (adjacencies != null && !adjacencies.isAdjacentTo(x, y)) {
                                continue;
                            }

//                            double bump;

                            int child = hashIndices.get(y);
                            int parent = hashIndices.get(x);
                            double bump = gesScore.localScoreDiff(child, new int[]{}, parent);

                            if (gesScore.isEffectEdge(bump)) {
                                final Edge edge = Edges.undirectedEdge(x, y);
                                if (boundGraph != null && !boundGraph.isAdjacentTo(edge.getNode1(), edge.getNode2()))
                                    continue;
                                effectEdgesGraph.addEdge(edge);
                            }

                            if (bump > 0.0) {
                                Arrow arrow1 = new Arrow(bump, x, y, emptySet, emptySet);
                                Arrow arrow2 = new Arrow(bump, y, x, emptySet, emptySet);

                                sortedArrows.add(arrow1);
                                addLookupArrow(x, y, arrow1);

                                sortedArrows.add(arrow2);
                                addLookupArrow(y, x, arrow2);
                            }
                        }
                    }

                    return true;
                } else {
                    int mid = (to - from) / 2;

                    List<EffectTask> tasks = new ArrayList<>();

                    tasks.add(new EffectTask(chunk, from, from + mid));
                    tasks.add(new EffectTask(chunk, from + mid, to));

                    invokeAll(tasks);

                    return true;
                }
            }

        }

        buildIndexing(nodes);
        pool.invoke(new EffectTask(minChunk, 0, nodes.size()));

        long stop = System.currentTimeMillis();

        out.println("Elapsed getEffectEdges = " + (stop - start) + " ms");

        return effectEdgesGraph;
    }

    /**
     * Forward equivalence search.
     *
     * @param graph The graph in the state prior to the forward equivalence search.
     */
    private void fes(Graph graph) {
        TetradLogger.getInstance().log("info", "** FORWARD EQUIVALENCE SEARCH");

        while (!sortedArrows.isEmpty()) {
            Arrow arrow = sortedArrows.first();
            sortedArrows.remove(arrow);

            Node x = arrow.getA();
            Node y = arrow.getB();

            if (graph.isAdjacentTo(x, y)) {
                continue;
            }

            if (!arrow.getNaYX().equals(getNaYX(x, y, graph))) {
                continue;
            }

            if (!getTNeighbors(x, y, graph).containsAll(arrow.getHOrT())) {
                continue;
            }

            if (!validInsert(x, y, arrow.getHOrT(), arrow.getNaYX(), graph)) {
                continue;
            }

            Set<Node> t = arrow.getHOrT();
            double bump = arrow.getBump();

            if (!insert(x, y, t, graph, bump)) {
                continue;
            }
            score += bump;

            clearArrow(x, y);

            Set<Node> visited = rebuildPatternRestricted(graph, x, y);
            Set<Node> toProcess = new HashSet<>();

            for (Node node : visited) {
                final Set<Node> neighbors = getNeighbors(node, graph);
                final Set<Node> storedNeighbors = this.neighbors.get(node);

                if (!(neighbors.equals(storedNeighbors))) {
                    toProcess.add(node);
                    this.neighbors.put(node, neighbors);
                }
            }

            Edge xy = graph.getEdge(x, y);

            if (xy.pointsTowards(x)) {
                toProcess.add(x);
            } else if (xy.pointsTowards(y)) {
                toProcess.add(y);
            }

            reevaluateForward(graph, toProcess);

            storeGraph(graph);
        }
    }

    // Returns the set of nodes {x} U {y} U adj(x) U adj(y).
    private Set<Node> adjNodes(Graph graph, Node x, Node y) {
        Set<Node> adj = new HashSet<>();
        adj.addAll(graph.getAdjacentNodes(x));
        adj.addAll(graph.getAdjacentNodes(y));
        adj.add(x);
        adj.add(y);
        return adj;
    }

    /**
     * Backward equivalence search.
     *
     * @param graph The graph in the state after the forward equivalence search.
     */
    private void bes(Graph graph) {
        TetradLogger.getInstance().log("info", "** BACKWARD EQUIVALENCE SEARCH");

        neighbors.clear();

        while (!sortedArrows.isEmpty()) {
            Arrow arrow = sortedArrows.first();
            sortedArrows.remove(arrow);

            Node x = arrow.getA();
            Node y = arrow.getB();

            if (!graph.isAdjacentTo(x, y)) continue;

            if (!arrow.getNaYX().equals(getNaYX(x, y, graph))) {
                continue;
            }

            if (!validDelete(y, arrow.getHOrT(), arrow.getNaYX(), graph)) {
                continue;
            }

            Set<Node> h = arrow.getHOrT();
            double bump = arrow.getBump();

            delete(x, y, h, graph, bump, arrow.getNaYX());
            score += bump;

            Set<Node> visited = rebuildPatternRestricted(graph, x, y);
            Set<Node> toProcess = new HashSet<>();

            for (Node node : visited) {
                final Set<Node> nayx = getNeighbors(node, graph);
                final Set<Node> storedNayx = this.neighbors.get(node);

                if (!nayx.equals(storedNayx)) {
                    toProcess.add(node);
                    this.neighbors.put(node, nayx);
                }
            }

            toProcess.add(x);
            toProcess.add(y);

            storeGraph(graph);

            for (Node node : toProcess) {
                for (Node p : graph.getParents(node)) {
                    reevaluateBackward(graph, node, p);
                }
            }
        }
    }

    // Returns true if knowledge is not empty.
    private boolean existsKnowledge() {
        return !knowledge.isEmpty();
    }

    // Initiaizes the sorted arrows and lookup arrows lists for the backward search.
    private void initializeArrowsBackward(Graph graph) {
        sortedArrows.clear();
        lookupArrows.clear();

        for (Edge edge : graph.getEdges()) {
            Node x = edge.getNode1();
            Node y = edge.getNode2();

            if (existsKnowledge()) {
                if (!getKnowledge().noEdgeRequired(x.getName(), y.getName())) {
                    continue;
                }
            }

            if (Edges.isDirectedEdge(edge)) {
                calculateArrowsBackward(x, y, graph);
            } else {
                calculateArrowsBackward(x, y, graph);
                calculateArrowsBackward(y, x, graph);
            }
        }
    }

    // Calcuates new arrows based on changes in the graph for the forward search.
    private void reevaluateForward(final Graph graph, final Set<Node> nodes) {
        List<Node> _nodes = new ArrayList<>(nodes);

        List<OrderedPair<Node>> pairs = new ArrayList<>();

        for (final Node x : _nodes) {
            List<Node> adj;

            if (isFaithfulnessAssumed()) {
                adj = effectEdgesGraph.getAdjacentNodes(x);
            } else {
                adj = variables;
            }

            for (Node w : adj) {
                pairs.add(new OrderedPair<>(w, x));
            }
        }

        class AdjTask extends RecursiveTask<Boolean> {
            private final List<OrderedPair<Node>> pairs;
            private int from;
            private int to;

            public AdjTask(List<OrderedPair<Node>> pairs, int from, int to) {
                this.pairs = pairs;
                this.from = from;
                this.to = to;
            }

            @Override
            protected Boolean compute() {
                if (to - from <= 25) {
                    for (int _w = from; _w < to; _w++) {
                        final OrderedPair<Node> p = pairs.get(_w);
                        Node w = p.getFirst();
                        Node x = p.getSecond();

                        if (w == x) continue;

                        if (adjacencies != null && !(adjacencies.isAdjacentTo(w, x))) {
                            continue;
                        }

                        if (!graph.isAdjacentTo(w, x)) {
                            calculateArrowsForward(w, x, graph);
                        }
                    }

                    return true;
                } else {
                    int mid = (to - from) / 2;

                    List<AdjTask> tasks = new ArrayList<>();

                    tasks.add(new AdjTask(pairs, from, from + mid));
                    tasks.add(new AdjTask(pairs, from + mid, to));

                    invokeAll(tasks);

                    return true;
                }
            }
        }

        final AdjTask task = new AdjTask(pairs, 0, pairs.size());

        pool.invoke(task);

    }

    // Calculates the new arrows for an a->b edge.
    private void calculateArrowsForward(final Node a, final Node b, Graph graph) {
        if (isFaithfulnessAssumed() && !effectEdgesGraph.isAdjacentTo(a, b)) return;
        if (adjacencies != null && !adjacencies.isAdjacentTo(a, b)) return;
//        this.neighbors.put(b, getNeighbors(b, graph));

        if (existsKnowledge()) {
            if (getKnowledge().isForbidden(a.getName(), b.getName())) {
                return;
            }
        }

        clearArrow(a, b);

        final Set<Node> naYX = getNaYX(a, b, graph);

        if (!isClique(naYX, graph)) return;

        List<Node> TNeighbors = getTNeighbors(a, b, graph);

        final int _depth = Math.min(TNeighbors.size(), depth == -1 ? 1000 : depth);

        for (int i = 0; i <= _depth; i++) {
            final ChoiceGenerator gen = new ChoiceGenerator(TNeighbors.size(), i);
            int[] choice;
            List<Set<Node>> subsets = new ArrayList<>();

            while ((choice = gen.next()) != null) {
                Set<Node> T = GraphUtils.asSet(choice, TNeighbors);

                Set<Node> union = new HashSet<>(naYX);
                union.addAll(T);

                if (!isClique(union, graph)) continue;
                subsets.add(T);

                if (existsKnowledge()) {
                    if (!validSetByKnowledge(b, T)) {
                        continue;
                    }
                }

                double bump = insertEval(a, b, T, naYX, graph, hashIndices);

                if (bump > 0.0) {
                    Arrow arrow = new Arrow(bump, a, b, T, naYX);
                    sortedArrows.add(arrow);
                    addLookupArrow(a, b, arrow);
                }
            }
        }
    }

//    private void calculateArrowsForward2(final Node a, final Node b, final Graph graph) {
//        if (isFaithfulnessAssumed() && !effectEdgesGraph.isAdjacentTo(a, b)) return;
//        if (adjacencies != null && !adjacencies.isAdjacentTo(a, b)) return;
//
//        if (existsKnowledge()) {
//            if (getKnowledge().isForbidden(a.getName(), b.getName())) {
//                return;
//            }
//        }
//
//        final Set<Node> naYX = getNaYX(a, b, graph);
//        final List<Node> TNeighbors = getTNeighbors(a, b, graph);
//
////        if (!isClique(naYX, graph)) return;
//
//        final int _depth = Math.min(TNeighbors.size(), depth == -1 ? 1000 : depth);
//
//        clearArrow(a, b);
//
//        final DepthChoiceGenerator gen = new DepthChoiceGenerator(TNeighbors.size(), _depth);
//        int[] choice;
//
//        while ((choice = gen.next()) != null) {
//            Set<Node> T = GraphUtils.asSet(choice, TNeighbors);
//
//            Set<Node> union = new HashSet<>(T);
//            union.addAll(naYX);
//
//            if (!isClique(union, graph)) continue;
//
//            if (existsKnowledge()) {
//                if (!validSetByKnowledge(b, T)) {
//                    continue;
//                }
//            }
//
//            double bump = insertEval(a, b, T, naYX, graph, hashIndices);
//
//            if (bump > 0.0) {
//                Arrow arrow = new Arrow(bump, a, b, T, naYX);
//                sortedArrows.add(arrow);
//                addLookupArrow(a, b, arrow);
//            }
//        }
//    }

    // Reevaluates arrows after removing an edge from the graph.
    private void reevaluateBackward(final Graph graph, final Node x, final Node y) {
        class BackwardTask extends RecursiveTask<Boolean> {
            private List<Node> nodes;
            private Map<Node, Integer> hashIndices;
            private int chunk;
            private int from;
            private int to;

            public BackwardTask(List<Node> nodes, int chunk, int from, int to,
                                Map<Node, Integer> hashIndices) {
                this.nodes = new ArrayList<>(nodes);
                this.hashIndices = new HashMap<>(hashIndices);
                this.chunk = chunk;
                this.from = from;
                this.to = to;
            }

            @Override
            protected Boolean compute() {
                if (to - from <= chunk) {
                    for (int _w = from; _w < to; _w++) {
                        final Node w = nodes.get(_w);

                        if (w == x) continue;
                        if (w == y) continue;

                        calculateArrowsBackward(w, x, graph);
                        calculateArrowsBackward(w, y, graph);
                    }

                    return true;
                } else {
                    int mid = (to - from) / 2;

                    List<BackwardTask> tasks = new ArrayList<>();

                    tasks.add(new BackwardTask(nodes, chunk, from, from + mid, hashIndices));
                    tasks.add(new BackwardTask(nodes, chunk, from + mid, to, hashIndices));

                    invokeAll(tasks);

                    return true;
                }
            }
        }

        Set<Node> _adj = adjNodes(graph, x, y);
        final List<Node> adj = new ArrayList<>(_adj);

        pool.invoke(new BackwardTask(adj, minChunk, 0, adj.size(), hashIndices));
    }

    // Calculates the arrows for the removal in the backward direction.
    private void calculateArrowsBackward(Node a, Node b, Graph graph) {
        if (existsKnowledge()) {
            if (!getKnowledge().noEdgeRequired(a.getName(), b.getName())) {
                return;
            }
        }

        clearArrow(a, b);

        Set<Node> naYX = getNaYX(a, b, graph);
        List<Node> _naYX = new ArrayList<>(naYX);

        final int _depth = Math.min(_naYX.size(), depth == -1 ? 1000 : depth);

        for (int i = 0; i <= _depth; i++) {
            final ChoiceGenerator gen = new ChoiceGenerator(_naYX.size(), i);
            int[] choice;

            while ((choice = gen.next()) != null) {
                Set<Node> diff = GraphUtils.asSet(choice, _naYX);

                if (!isClique(diff, graph)) continue;

                Set<Node> h = new HashSet<>(_naYX);
                h.removeAll(diff);

                if (existsKnowledge()) {
                    if (!validSetByKnowledge(b, h)) {
                        continue;
                    }
                }

                double bump = deleteEval(a, b, diff, naYX, graph, hashIndices);

                if (bump >= 0.0) {
                    Arrow arrow = new Arrow(bump, a, b, h, naYX);
                    sortedArrows.add(arrow);
                    addLookupArrow(a, b, arrow);
                }
            }
        }
    }

//    private void calculateArrowsBackward(Node a, Node b, Graph graph) {
//        if (a == b) {
//            return;
//        }
//
//        if (!graph.isAdjacentTo(a, b)) {
//            return;
//        }
//
//        Edge e = graph.getEdge(a, b);
//
//        if (existsKnowledge()) {
//            if (!getKnowledge().noEdgeRequired(a.getName(), b.getName())) {
//                return;
//            }
//        }
//
//        Set<Node> naYX = getNaYX(a, b, graph);
//
//        clearArrow(a, b);
//
//        List<Node> _naYX = new ArrayList<>(naYX);
//
//        DepthChoiceGenerator gen = new DepthChoiceGenerator(_naYX.size(), _naYX.size());
//        int[] choice;
//
//        while ((choice = gen.next()) != null) {
//            Set<Node> H = GraphUtils.asSet(choice, _naYX);
//
//            Set<Node> diff = new HashSet<>(naYX);
//            diff.removeAll(H);
//
//            if (!isClique(diff, graph)) continue;
//
//            if (existsKnowledge()) {
//                if (!validSetByKnowledge(b, H)) {
//                    continue;
//                }
//            }
//
//            double bump = deleteEval(a, b, H, naYX, graph, hashIndices);
//
//            if (bump >= 0.0) {
//                Arrow arrow = new Arrow(bump, a, b, H, naYX);
//                sortedArrows.add(arrow);
//                addLookupArrow(a, b, arrow);
//            }
//        }
//    }

    public void setSamplePrior(double samplePrior) {
        if (gesScore instanceof LocalDiscreteScore) {
            ((LocalDiscreteScore) gesScore).setSamplePrior(samplePrior);
        }
    }

    public void setStructurePrior(double structurePrior) {
        if (gesScore instanceof LocalDiscreteScore) {
            ((LocalDiscreteScore) gesScore).setStructurePrior(structurePrior);
        }
    }

    /**
     * The top n graphs found by the algorithm, where n is numPatternsToStore.
     */
    public LinkedList<ScoredGraph> getTopGraphs() {
        return topGraphs;
    }

    // Basic data structure for an arrow a->b considered for additiom or removal from the graph, together with
    // associated sets needed to make this determination. For both forward and backward direction, NaYX is needed.
    // For the forward direction, T neighbors are needed; for the backward direction, H neighbors are needed.
    // See Chickering (2002). The score difference resulting from added in the edge (hypothetically) is recorded
    // as the "bump".
    private static class Arrow implements Comparable<Arrow> {
        private double bump;
        private Node a;
        private Node b;
        private Set<Node> hOrT;
        private Set<Node> naYX;

        public Arrow(double bump, Node a, Node b, Set<Node> hOrT, Set<Node> naYX) {
            this.bump = bump;
            this.a = a;
            this.b = b;
            this.hOrT = hOrT;
            this.naYX = naYX;
        }

        public double getBump() {
            return bump;
        }

        public Node getA() {
            return a;
        }

        public Node getB() {
            return b;
        }

        public Set<Node> getHOrT() {
            return hOrT;
        }

        public Set<Node> getNaYX() {
            return naYX;
        }

        // Sorting by bump, high to low. The problem is the SortedSet contains won't add a new element if it compares
        // to zero with an existing element, so for the cases where the comparison is to zero (i.e. have the same
        // bump, we need to determine as quickly as possible a determinate ordering (fixed) ordering for two variables.
        // The fastest way to do this is using a hash code, though it's still possible for two Arrows to have the
        // same hash code but not be equal. If we're paranoid, in this case we calculate a determinate comparison
        // not equal to zero by keeping a list. This last part is commened out by default.
        public int compareTo(Arrow arrow) {
            final int compare = Double.compare(arrow.getBump(), getBump());

            if (compare == 0) {
                int hashcode1 = hashCode();
                int hashcode2 = arrow.hashCode();
                return Integer.compare(hashcode1, hashcode2);
            }

            return compare;
        }

        public boolean equals(Object o) {
            if (!(o instanceof Arrow)) {
                return false;
            }

            Arrow a = (Arrow) o;

            return a.a.equals(this.a) && a.b.equals(this.b) && a.hOrT.equals(this.hOrT) && a.naYX.equals(this.naYX);
        }

        public int hashCode() {
            return 11 * a.hashCode() + 13 * b.hashCode() + 17 * hOrT.hashCode() + 19 * naYX.hashCode();
        }

        public String toString() {
            return "Arrow<" + a + "->" + b + " bump = " + bump + " t/h = " + hOrT + " naYX = " + naYX + ">";
        }
    }

    // Get all nodes that are connected to Y by an undirected edge and not adjacent to X.
    private static List<Node> getTNeighbors(Node x, Node y, Graph graph) {
        List<Edge> yEdges = graph.getEdges(y);
        List<Node> tNeighbors = new ArrayList<>();

        for (Edge edge : yEdges) {
            if (!Edges.isUndirectedEdge(edge)) {
                continue;
            }

            Node z = edge.getDistalNode(y);

            if (graph.isAdjacentTo(z, x)) {
                continue;
            }

            tNeighbors.add(z);
        }

        return tNeighbors;
    }

    // Get all nodes that are connected to Y by an undirected edge.
    private static Set<Node> getNeighbors(Node y, Graph graph) {
        List<Edge> yEdges = graph.getEdges(y);
        Set<Node> neighbors = new HashSet<>();

        for (Edge edge : yEdges) {
            if (!Edges.isUndirectedEdge(edge)) {
                continue;
            }

            Node z = edge.getDistalNode(y);

            neighbors.add(z);
        }

        return neighbors;
    }

    // Evaluate the Insert(X, Y, T) operator (Definition 12 from Chickering, 2002).
    private double insertEval(Node x, Node y, Set<Node> t, Set<Node> naYX, Graph graph,
                              Map<Node, Integer> hashIndices) {
        Set<Node> parents = new HashSet<>(graph.getParents(y));

        Set<Node> set = new HashSet<>(naYX);
        set.addAll(t);
        set.addAll(parents);
        set.remove(x);

        return scoreGraphChange(y, set, x, hashIndices);
    }

    // Evaluate the Delete(X, Y, T) operator (Definition 12 from Chickering, 2002).
    private double deleteEval(Node x, Node y, Set<Node> diff, Set<Node> naYX, Graph graph,
                              Map<Node, Integer> hashIndices) {
        Set<Node> parents = new HashSet<>(graph.getParents(y));
//        parents.remove(x);

        Set<Node> set = new HashSet<>(diff);
        set.addAll(parents);
        set.remove(x);

        return -scoreGraphChange(y, set, x, hashIndices);
    }


    // Do an actual insertion. (Definition 12 from Chickering, 2002).
    private boolean insert(Node x, Node y, Set<Node> t, Graph graph, double bump) {
        if (graph.isAdjacentTo(x, y)) {
            return false; // The initial graph may already have put this edge in the graph.
//            throw new IllegalArgumentException(x + " and " + y + " are already adjacent in the graph.");
        }

        Edge trueEdge = null;

        if (trueGraph != null) {
            Node _x = trueGraph.getNode(x.getName());
            Node _y = trueGraph.getNode(y.getName());
            trueEdge = trueGraph.getEdge(_x, _y);
        }

        if (boundGraph != null && !boundGraph.isAdjacentTo(x, y)) return false;
        graph.addDirectedEdge(x, y);

        if (log) {
            String label = trueGraph != null && trueEdge != null ? "*" : "";
            TetradLogger.getInstance().log("insertedEdges", graph.getNumEdges() + ". INSERT " + graph.getEdge(x, y) +
                    " " + t + " " + bump + " " + label);
        }

        int numEdges = graph.getNumEdges();
        if (numEdges % 1000 == 0) out.println("Num edges added: " + numEdges);

        if (verbose) {
            String label = trueGraph != null && trueEdge != null ? "*" : "";
            out.println(graph.getNumEdges() + ". INSERT " + graph.getEdge(x, y) +
                    " " + t + " " + bump + " " + label);
        }

        for (Node _t : t) {
            Edge oldEdge = graph.getEdge(_t, y);

            if (oldEdge == null) throw new IllegalArgumentException("Not adjacent: " + _t + ", " + y);

            graph.removeEdge(_t, y);
            if (boundGraph != null && !boundGraph.isAdjacentTo(_t, y)) continue;
            graph.addDirectedEdge(_t, y);

            if (log && verbose) {
                TetradLogger.getInstance().log("directedEdges", "--- Directing " + oldEdge + " to " +
                        graph.getEdge(_t, y));
                out.println("--- Directing " + oldEdge + " to " +
                        graph.getEdge(_t, y));
            }
        }

        return true;
    }

    // Do an actual deletion (Definition 13 from Chickering, 2002).
    private void delete(Node x, Node y, Set<Node> H, Graph graph, double bump, Set<Node> naYX) {

        Edge trueEdge = null;

        if (trueGraph != null) {
            Node _x = trueGraph.getNode(x.getName());
            Node _y = trueGraph.getNode(y.getName());
            trueEdge = trueGraph.getEdge(_x, _y);
        }

        Edge oldxy = graph.getEdge(x, y);

        graph.removeEdge(x, y);
//        rebuildPatternRestricted(graph, x, y);


        if (verbose) {
            int numEdges = graph.getNumEdges();
            if (numEdges % 1000 == 0) out.println("Num edges (backwards) = " + numEdges);
        }

        Set<Node> diff = new HashSet<>(naYX);
        diff.removeAll(H);

        if (verbose) {
            int numEdges = graph.getNumEdges();
            if (numEdges % 1000 == 0) out.println("Num edges (backwards) = " + numEdges);

            String label = trueGraph != null && trueEdge != null ? "*" : "";
            String message = (graph.getNumEdges()) + ". DELETE " + oldxy +
                    " H = " + H + " NaYX = " + naYX + " diff = " + diff + " (" + bump + ") " + label;
            TetradLogger.getInstance().log("deletedEdges", message);
            out.println(message);
        }

        for (Node h : H) {
            Edge oldyh = graph.getEdge(y, h);
//
            if (oldyh != null && Edges.isUndirectedEdge(oldyh)) {
                if (!graph.isAdjacentTo(y, h)) throw new IllegalArgumentException("Not adjacent: " + y + ", " + h);

                graph.removeEdge(y, h);
                graph.addDirectedEdge(y, h);

                Edge edge = graph.getEdge(y, h);

                if (log && verbose) {
                    TetradLogger.getInstance().log("directedEdges", "--- Directing " + oldyh + " to " +
                            edge);
                    out.println("--- Directing " + oldyh + " to " + edge);
                }
            }

            Edge oldxh = graph.getEdge(x, h);

            if (oldxh != null && Edges.isUndirectedEdge(oldxh)) {
                if (!graph.isAdjacentTo(x, h)) throw new IllegalArgumentException("Not adjacent: " + x + ", " + h);

                graph.removeEdge(x, h);
                graph.addDirectedEdge(x, h);

                Edge edge = graph.getEdge(x, h);

                if (log && verbose) {
                    TetradLogger.getInstance().log("directedEdges", "--- Directing " + oldxh + " to " +
                            edge);
                    out.println("--- Directing " + oldxh + " to " + edge);
                }
            }
        }
    }

    // Test if the candidate insertion is a valid operation
    // (Theorem 15 from Chickering, 2002).
    private boolean validInsert(Node x, Node y, Set<Node> s, Set<Node> naYX, Graph graph) {
        Set<Node> union = new HashSet<>(s);
        union.addAll(naYX);

        // Note s U NaYX must be a clique, but this has already been checked. Nevertheless, at this
        // point it must be verified that all nodes in s U NaYX are neighbors of Y, since some of
        // the edges g---Y may have been oriented in the interim.
        return allNeighbors(y, union, graph) && !existsUnblockedSemiDirectedPath(y, x, union, graph, cycleBound);

    }

    // Returns true if all of the members of 'union' are neighbors of y.
    private boolean allNeighbors(Node y, Set<Node> union, Graph graph) {
        for (Node n : union) {
            Edge e = graph.getEdge(y, n);
            if (e == null) {
                return false;
            }
            if (!Edges.isUndirectedEdge(e)) {
                return false;
            }
        }

        return true;
    }

    // Test if the candidate deletion is a valid operation (Theorem 17 from Chickering, 2002).
    private boolean validDelete(Node y, Set<Node> h, Set<Node> naXY, Graph graph) {
        Set<Node> set = new HashSet<>(naXY);
        set.removeAll(h);
        return isClique(set, graph) && allNeighbors(y, set, graph);
    }

    // Adds edges required by knowledge.
    private void addRequiredEdges(Graph graph) {
        if (!existsKnowledge()) return;

        for (Iterator<KnowledgeEdge> it = getKnowledge().requiredEdgesIterator(); it.hasNext(); ) {
            KnowledgeEdge next = it.next();

            Node nodeA = graph.getNode(next.getFrom());
            Node nodeB = graph.getNode(next.getTo());

            if (!graph.isAncestorOf(nodeB, nodeA)) {
                graph.removeEdges(nodeA, nodeB);
                graph.addDirectedEdge(nodeA, nodeB);
                TetradLogger.getInstance().log("insertedEdges", "Adding edge by knowledge: " + graph.getEdge(nodeA, nodeB));
            }
        }
        for (Edge edge : graph.getEdges()) {
            final String A = edge.getNode1().getName();
            final String B = edge.getNode2().getName();

            if (knowledge.isForbidden(A, B)) {
                Node nodeA = edge.getNode1();
                Node nodeB = edge.getNode2();
                if (nodeA == null || nodeB == null) throw new NullPointerException();

                if (graph.isAdjacentTo(nodeA, nodeB) && !graph.isChildOf(nodeA, nodeB)) {
                    if (!graph.isAncestorOf(nodeA, nodeB)) {
                        graph.removeEdges(nodeA, nodeB);
                        graph.addDirectedEdge(nodeB, nodeA);
                        TetradLogger.getInstance().log("insertedEdges", "Adding edge by knowledge: " + graph.getEdge(nodeB, nodeA));
                    }
                }

                if (!graph.isChildOf(nodeA, nodeB) && getKnowledge().isForbidden(nodeA.getName(), nodeB.getName())) {
                    if (!graph.isAncestorOf(nodeA, nodeB)) {
                        graph.removeEdges(nodeA, nodeB);
                        graph.addDirectedEdge(nodeB, nodeA);
                        TetradLogger.getInstance().log("insertedEdges", "Adding edge by knowledge: " + graph.getEdge(nodeB, nodeA));
                    }
                }
            } else if (knowledge.isForbidden(B, A)) {
                Node nodeA = edge.getNode2();
                Node nodeB = edge.getNode1();
                if (nodeA == null || nodeB == null) throw new NullPointerException();

                if (graph.isAdjacentTo(nodeA, nodeB) && !graph.isChildOf(nodeA, nodeB)) {
                    if (!graph.isAncestorOf(nodeA, nodeB)) {
                        graph.removeEdges(nodeA, nodeB);
                        graph.addDirectedEdge(nodeB, nodeA);
                        TetradLogger.getInstance().log("insertedEdges", "Adding edge by knowledge: " + graph.getEdge(nodeB, nodeA));
                    }
                }
                if (!graph.isChildOf(nodeA, nodeB) && getKnowledge().isForbidden(nodeA.getName(), nodeB.getName())) {
                    if (!graph.isAncestorOf(nodeA, nodeB)) {
                        graph.removeEdges(nodeA, nodeB);
                        graph.addDirectedEdge(nodeB, nodeA);
                        TetradLogger.getInstance().log("insertedEdges", "Adding edge by knowledge: " + graph.getEdge(nodeB, nodeA));
                    }
                }
            }
        }
    }

    // Use background knowledge to decide if an insert or delete operation does not orient edges in a forbidden
    // direction according to prior knowledge. If some orientation is forbidden in the subset, the whole subset is
    // forbidden.
    private boolean validSetByKnowledge(Node y, Set<Node> subset) {
        for (Node node : subset) {
            if (getKnowledge().isForbidden(node.getName(), y.getName())) {
                return false;
            }
        }
        return true;
    }

    // Find all nodes that are connected to Y by an undirected edge that are adjacent to X (that is, by undirected or
    // directed edge).
    private static Set<Node> getNaYX(Node x, Node y, Graph graph) {
        List<Edge> yEdges = graph.getEdges(y);
        Set<Node> nayx = new HashSet<>();

        for (Edge edge : yEdges) {
            if (!Edges.isUndirectedEdge(edge)) {
                continue;
            }

            Node z = edge.getDistalNode(y);

            if (!graph.isAdjacentTo(z, x)) {
                continue;
            }

            nayx.add(z);
        }

        return nayx;
    }

    // Returns true iif the given set forms a clique in the given graph.
    private static boolean isClique(Set<Node> nodes, Graph graph) {
        List<Node> _nodes = new ArrayList<>(nodes);
        for (int i = 0; i < _nodes.size() - 1; i++) {
            for (int j = i + 1; j < _nodes.size(); j++) {
                if (!graph.isAdjacentTo(_nodes.get(i), _nodes.get(j))) {
                    return false;
                }
            }
        }

        return true;
    }

    private static boolean isClique(Set<Node> clique, Node newNode, Graph graph) {
        for (Node z : clique) {
            if (!graph.isAdjacentTo(z, newNode)) {
                return false;
            }
        }

        return true;
    }

    // Returns true if a path consisting of undirected and directed edges toward 'to' exists of
    // length at most 'bound'. Cycle checker in other words.
    private boolean existsUnblockedSemiDirectedPath(Node from, Node to, Set<Node> cond, Graph G, int bound) {
        Queue<Node> Q = new LinkedList<>();
        Set<Node> V = new HashSet<>();
        Q.offer(from);
        V.add(from);
        Node e = null;
        int distance = 0;

        while (!Q.isEmpty()) {
            Node t = Q.remove();
            if (t == to) return true;

            if (e == t) {
                e = null;
                distance++;
                if (distance > (bound == -1 ? 1000 : bound)) return true;
            }

            for (Node u : G.getAdjacentNodes(t)) {
                Edge edge = G.getEdge(t, u);
                Node c = traverseSemiDirected(t, edge);
                if (c == null) continue;
                if (cond.contains(c)) continue;
                if (c == to) return true;

                if (!V.contains(c)) {
                    V.add(c);
                    Q.offer(c);

                    if (e == null) {
                        e = u;
                    }
                }
            }

        }

        return false;
    }

    // Used to find semidirected paths for cycle checking.
    private static Node traverseSemiDirected(Node node, Edge edge) {
        if (node == edge.getNode1()) {
            if (edge.getEndpoint1() == Endpoint.TAIL) {
                return edge.getNode2();
            }
        } else if (node == edge.getNode2()) {
            if (edge.getEndpoint2() == Endpoint.TAIL) {
                return edge.getNode1();
            }
        }
        return null;
    }

    // Runs the Meek rules on just the changed nodes.
    private Set<Node> rebuildPatternRestricted(Graph graph, Node x, Node y) {
        Set<Node> visited = new HashSet<>();

        visited.addAll(reorientNode(graph, x));
        visited.addAll(reorientNode(graph, y));

        if (TetradLogger.getInstance().isEventActive("rebuiltPatterns")) {
            TetradLogger.getInstance().log("rebuiltPatterns", "Rebuilt pattern = " + graph);
        }

        return visited;
    }

    // Runs Meek rules on just the changed nodes.
    private Set<Node> reorientNode(Graph graph, Node a) {
        List<Node> nodes = graph.getAdjacentNodes(a);
        nodes.add(a);

        List<Edge> edges = graph.getEdges(a);
        SearchGraphUtils.basicPatternRestricted2(a, graph);
        addRequiredEdges(graph);
        Set<Node> visited = meekOrientRestricted(graph, nodes, getKnowledge());

        List<Edge> newEdges = graph.getEdges(a);
        newEdges.removeAll(edges); // The newly oriented edges.

        for (Edge edge : newEdges) {
            if (Edges.isUndirectedEdge(edge)) {
                Node _node = edge.getDistalNode(a);
                visited.addAll(reorientNode(graph, _node));
            }
        }

        return visited;
    }

    // Runs Meek rules on just the changed nodes.
    private Set<Node> meekOrientRestricted(Graph graph, List<Node> nodes, IKnowledge knowledge) {
        MeekRules rules = new MeekRules();
        rules.setKnowledge(knowledge);
        rules.orientImplied(graph, nodes);
        return rules.getVisited();
    }

    // Maps nodes to their indices for quick lookup.
    private void buildIndexing(List<Node> nodes) {
        this.hashIndices = new ConcurrentHashMap<>();
        for (Node node : nodes) {
            this.hashIndices.put(node, variables.indexOf(node));
        }
    }

    // Removes information associated with an edge x->y.
    private void clearArrow(Node x, Node y) {
        final OrderedPair<Node> pair = new OrderedPair<>(x, y);
        final Set<Arrow> lookupArrows = this.lookupArrows.get(pair);

        if (lookupArrows != null) {
            sortedArrows.removeAll(lookupArrows);
        }

        this.lookupArrows.remove(pair);
    }

    // Adds the given arrow for the adjacency i->j. These all are for i->j but may have
    // different T or H or NaYX sets, and so different bumps.
    private void addLookupArrow(Node i, Node j, Arrow arrow) {
        OrderedPair<Node> pair = new OrderedPair<>(i, j);
        Set<Arrow> arrows = lookupArrows.get(pair);

        if (arrows == null) {
            arrows = new ConcurrentSkipListSet<>();
            lookupArrows.put(pair, arrows);
        }

        arrows.add(arrow);
    }

    //===========================SCORING METHODS===================//

    /**
     * Scores the given DAG, up to a constant.
     */
    public double scoreDag(Graph dag) {
        buildIndexing(dag.getNodes());

        double score = 0.0;

        for (Node y : dag.getNodes()) {
            Set<Node> parents = new HashSet<>(dag.getParents(y));
            int index = hashIndices.get(y);
            int parentIndices[] = new int[parents.size()];
            Iterator<Node> pi = parents.iterator();
            int count = 0;
            while (pi.hasNext()) {
                Node nextParent = pi.next();
                parentIndices[count++] = hashIndices.get(nextParent);
            }

            score += gesScore.localScore(index, parentIndices);
        }
        return score;
    }

    private double scoreGraphChange(Node y, Set<Node> parents,
                                    Node x, Map<Node, Integer> hashIndices) {
        int yIndex = hashIndices.get(y);

        if (parents.contains(x)) throw new IllegalArgumentException();

        int[] parentIndices = new int[parents.size()];

        int count = 0;
        for (Node parent : parents) {
            parentIndices[count++] = hashIndices.get(parent);
        }

        return gesScore.localScoreDiff(yIndex, parentIndices, hashIndices.get(x));
    }

    private int sampleSize() {
        return this.sampleSize;
    }

    private List<Node> getVariables() {
        return variables;
    }

    private boolean isDiscrete() {
        return discrete;
    }

    // Stores the graph, if its score knocks out one of the top ones.
    private void storeGraph(Graph graph) {
        if (getNumPatternsToStore() > 0) {
            Graph graphCopy = new EdgeListGraphSingleConnections(graph);
            topGraphs.addLast(new ScoredGraph(graphCopy, score));
        }

        if (topGraphs.size() == getNumPatternsToStore() + 1) {
            topGraphs.removeFirst();
        }
    }
}






