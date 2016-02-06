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
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.IndependenceFact;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.DepthChoiceGenerator;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.cmu.tetrad.util.TetradVector;

import java.io.PrintStream;
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
        this.variables = dag.getNodes();
    }

    @Override
    public double localScore(int i, int[] parents, int[] all) {
        Node child = variables.get(i);
        List<Node> scoreParents = getVariableList(parents);
        int numAgreements = 0;

        List<Node> P = dag.getParents(child);

        for (int j = 0; j < scoreParents.size(); j++) {
            Node x = scoreParents.get(j);

            if (dag.isDConnectedTo(x, child, P)) {
                numAgreements++;
            } else {
                numAgreements--;
            }
        }

        int score = numAgreements;

        System.out.println("Score for " + variables.get(i) + " | " + getVariableSet(parents) + " is " + (score));

        return score;
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
    public double localScoreDiff(int i, int[] parents, int extra) {


        Node child = variables.get(i);
        Node newVar = variables.get(extra);
        List<Node> scoreParents = getVariableList(parents);
        List<Node> allvars = new ArrayList<>(scoreParents);

        double diff = -1;

        List<Node> childParents = dag.getParents(child);
        childParents.retainAll(allvars);

        if (!dag.isDSeparatedFrom(newVar, child, childParents)) {
            diff = 1;
        }

//        for (Node node : allvars) {
//            List<Node> nodeParents = dag.getParents(node);
//            nodeParents.retainAll(allvars);
//
//            if (dag.isDSeparatedFrom(newVar, node, nodeParents)) {
//                diff = -1;
//            }
//        }

        System.out.println("Score diff for " + newVar + "-->" + child + " given " + scoreParents + " = " + diff);

        return diff;
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



