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
                            !dag.isDConnectedTo(x, z, scoreParents) &&
                            dag.isDConnectedTo(z, y, scoreParents) &&
                            dag.isDConnectedTo(x, z, yUnionScoreParents)
                    )
            {
                diff += 1;
            }
        }

//        System.out.println("Score diff for " + x + "-->" + y + " given " + scoreParents + " = " + diff);

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



