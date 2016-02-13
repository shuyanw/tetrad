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

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeType;

import java.util.*;

/**
 * Implements the continuous BIC score for FGS.
 *
 * @author Joseph Ramsey
 */
public class GraphScore implements GesScore {

    private final Graph dag;

    // The variables of the covariance matrix.
    private List<Node> variables;

    // True if verbose output should be sent to out.
    private boolean verbose = false;

    /**
     * Constructs the score using a covariance matrix.
     */
    public GraphScore(Graph dag) {
        this.dag = dag;

        this.variables = new ArrayList<>();

        for (Node node : dag.getNodes()) {
            if (node.getNodeType() == NodeType.MEASURED) {
                this.variables.add(node);
            }
        }

        Collections.sort(variables);
    }

    /**
     * Calculates the sample likelihood and BIC score for i given its parents in a simple SEM model
     */
    public double localScore(int i, int[] parents) {
        Node child = variables.get(i);
        Set<Node> scoreParents = getVariableSet(parents);
        Set<Node> graphParents = new HashSet<>(dag.getParents(child));
        graphParents.retainAll(scoreParents);

        Set<Node> missing = new HashSet<>(scoreParents);
        missing.removeAll(graphParents);

        Set<Node> extra = new HashSet<>(dag.getParents(child));
        extra.removeAll(scoreParents);

        int error = -(missing.size() + extra.size());


        return error;
    }


    private Set<Node> getVariableSet(int[] indices) {
        Set<Node> variables = new HashSet<>();
        for (int i : indices) {
            variables.add(this.variables.get(i));
        }
        return variables;
    }

    private List<Node> getVariableList(int[] indices) {
        List<Node> variables = new ArrayList<>();
        for (int i : indices) {
            variables.add(this.variables.get(i));
        }
        return variables;
    }


    @Override
    public synchronized double localScoreDiff(int i, int[] parents, int extra) {


        Node y = variables.get(i);
        Node x = variables.get(extra);
        List<Node> scoreParents = getVariableList(parents);
        List<Node> allvars = new ArrayList<>(scoreParents);
        allvars.add(x);

//
//        double diff = score1(y, x, scoreParents);
//        double diff = score2(y, x, scoreParents);
        double diff = score3(x, y, scoreParents);
//        double diff = score4  (y, x, scoreParents);

//        System.out.println("Score diff for " + x + "-->" + y + " given " + scoreParents + " = " + diff);

        return diff;
    }

    private double score1(Node y, Node x, List<Node> scoreParents) {
        double diff = 0;

        if (dag.isDConnectedTo(x, y, scoreParents)) {
            diff += 1;
        }
        return diff;
    }

    private double score2(Node y, Node x, List<Node> scoreParents) {
        double diff = 0;

        if (dag.isDConnectedTo(x, y, scoreParents)) {
            diff += 1;
        }

        List<Node> yUnionScoreParents = new ArrayList<>();
        yUnionScoreParents.add(y);
        yUnionScoreParents.addAll(scoreParents);

        for (Node z : scoreParents) {
            if (
                    dag.isDConnectedTo(x, y, scoreParents) &&
                            dag.isDConnectedTo(z, y, scoreParents) &&
//                            !dag.isDConnectedTo(x, z, scoreParents) &&
                            dag.isDConnectedTo(z, y, scoreParents) &&
                            dag.isDConnectedTo(x, z, yUnionScoreParents)
                    ) {
                diff += 1;
                break;
            }
        }


        return diff;
    }

    private double scoreb(Node y, Node x, List<Node> scoreParents) {
        double diff = 0;

        if (dag.isDConnectedTo(x, y, Collections.EMPTY_LIST)) {
            diff += 1;
        }

        List<Node> yUnionScoreParents = new ArrayList<>();
        yUnionScoreParents.add(y);
        yUnionScoreParents.addAll(scoreParents);

        for (Node z : scoreParents) {
            if (
                    dag.isDConnectedTo(x, y, Collections.EMPTY_LIST) &&
                            dag.isDConnectedTo(z, y, Collections.EMPTY_LIST) &&
//                            !dag.isDConnectedTo(x, z, scoreParents) &&
                            dag.isDConnectedTo(z, y, Collections.EMPTY_LIST) &&
                            dag.isDConnectedTo(x, z, Collections.singletonList(y))
                    ) {
                diff += 1;
                break;
            }
        }


        return diff;
    }


    private double score3(Node x, Node y, List<Node> scoreParents) {
        List<Node> yUnionScoreParents = new ArrayList<>();
        yUnionScoreParents.add(y);
        yUnionScoreParents.addAll(scoreParents);
        int numDependencies = 0;
        int numIndependencies = 0;

        for (Node z : scoreParents) {
            List<Node> _scoreParents = new ArrayList<>(scoreParents);
            _scoreParents.remove(z);

            if (dag.isDConnectedTo(z, y, new ArrayList<>(scoreParents))) {
                numDependencies++;
            } else {
                numIndependencies++;
            }

            if (dag.isDConnectedTo(x, z, yUnionScoreParents)) {
                numDependencies++;
            } else {
                numIndependencies++;
            }
        }

        boolean xyConnected = dag.isDConnectedTo(x, y, scoreParents);

        if (xyConnected) {
            numDependencies++;
        } else {
            numIndependencies++;
        }

        if (isClique(scoreParents, dag)) {
            numIndependencies += 1;
        }

        if (xyConnected) {
            return numDependencies;
        } else {
            return -numIndependencies;
        }
    }

    private boolean isClique(List<Node> nodes, Graph graph) {
        for (int i = 0; i < nodes.size() - 1; i++) {
            for (int j = i + 1; j < nodes.size(); j++) {
                if (!graph.isAdjacentTo(nodes.get(i), nodes.get(j))) {
                    return false;
                }
            }
        }

        return true;
    }

    private double score4(Node x, Node y, List<Node> scoreParents) {
        double score = 0;

        List<Node> yUnionScoreParents = new ArrayList<>();
        yUnionScoreParents.add(y);
        yUnionScoreParents.addAll(scoreParents);

        if (dag.isDConnectedTo(x, y, yUnionScoreParents)) {
            score = 1;
        } else {
            score = -1;
        }

        double up = 2;
        double down = .5;

        for (Node z : scoreParents) {
//            if (dag.isDConnectedTo(x, z, yUnionScoreParents)) {
//                score *= 2;
//            } else {
//                score /= 3;
//            }

//            if (dag.isDConnectedTo(z, y, scoreParents)) {
//                score *= up;
//            } else {
//                score *= down;
//            }

            if (dag.isDConnectedTo(x, z, yUnionScoreParents)) {
                score *= up;
            } else {
                score *= down;
            }
        }

        System.out.println("Score for " + y + " given " + x + " and " + scoreParents + " is " + score);

        return score;
    }

    // Peter's score.
    private double score5(Node y, Node x, List<Node> scoreParents) {
        double diff = 0;

        if (dag.isDSeparatedFrom(x, y, scoreParents)) {
            diff += 1;
        }

        for (Node t : scoreParents) {
            List<Node> sepset = dag.getSepset(x, t);

            if (sepset != null && !sepset.contains(y)) {
                diff += 1;
            } else {
                diff -= 1;
            }
        }

        return diff;
    }


    private List<Node> getParentsInDag(Node y, List<Node> allvars) {
        return allvars;
//        Set<Node> parents = new HashSet<>();
//        parents.addAll(dag.getParents(y));
//        parents.retainAll(allvars);
//        return new ArrayList<>(parents);
    }

    int[] append(int[] parents, int extra) {
        int[] all = new int[parents.length + 1];
        System.arraycopy(parents, 0, all, 0, parents.length);
        all[parents.length] = extra;
        return all;
    }

    /**
     * Specialized scoring method for a single parent. Used to speed up the effect edges search.
     */

    public double localScore(int i, int parent) {
        return localScore(i, new int[]{parent});
    }

    /**
     * Specialized scoring method for no parents. Used to speed up the effect edges search.
     */
    public double localScore(int i) {
        return localScore(i, new int[]{});
    }

    @Override
    public boolean isEffectEdge(double bump) {
        return true;
    }

    public DataSet getDataSet() {
        throw new UnsupportedOperationException();
    }

    public boolean isVerbose() {
        return verbose;
    }

    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }

    @Override
    public List<Node> getVariables() {
        return variables;
    }

    @Override
    public int getSampleSize() {
        return 0;
    }

    @Override
    public boolean isDiscrete() {
        return false;
    }
}



