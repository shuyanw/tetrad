package edu.pitt.dbmi.algo.bootstrap;

import java.util.ArrayList;
import java.util.List;

import edu.cmu.tetrad.bayes.BayesIm;
import edu.cmu.tetrad.bayes.BayesPm;
import edu.cmu.tetrad.bayes.MlBayesIm;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DiscreteVariable;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Parameters;

/**
 * 
 * Apr 17, 2017 11:13:21 AM
 * 
 * @author Chirayu (Kong) Wongchokprasitti, PhD
 * 
 */
public class TestBootstrapTest {

    /**
     * @param args
     */
    public static void main(String[] args) {
	double structurePrior = 1, samplePrior = 1;
	boolean faithfulnessAssumed = false;
	
	int numVars = 20;
	int edgesPerNode = 2;
	int numCases = 1000;
	int numBootstrapSamples = 100;
	boolean verbose = true;

	Graph dag = makeDAG(numVars, edgesPerNode);

	BayesPm pm = new BayesPm(dag, 2, 3);
	BayesIm im = new MlBayesIm(pm, MlBayesIm.RANDOM);

	System.out.println("Generating data");

	DataSet data = im.simulateData(numCases, 0, false);

	Parameters parameters = new Parameters();
	parameters.set("structurePrior", structurePrior);
	parameters.set("samplePrior", samplePrior);
	parameters.set("faithfulnessAssumed", faithfulnessAssumed);
	parameters.set("numPatternsToStore", 0);
	parameters.set("verbose", false);
	
	BootstrapTest bootstrapTest = new BootstrapTest(data, BootstrapAlgName.GFCI, numBootstrapSamples);
	bootstrapTest.setVerbose(verbose);
	bootstrapTest.setParameters(parameters);
	Graph resultGraph = bootstrapTest.search();
	System.out.println(resultGraph.toString());
	
    }

    private static Graph makeDAG(int numVars, double edgesPerNode) {
	final int numEdges = (int) (numVars * edgesPerNode);

	//System.out.println("Making list of vars");

	List<Node> vars = new ArrayList<>();

	for (int i = 0; i < numVars; i++) {
	    vars.add(new DiscreteVariable(Integer.toString(i)));
	}

	//System.out.println("Making dag");
	return GraphUtils.randomGraph(vars, 0, numEdges, 30, 15, 15, false);// randomGraphRandomForwardEdges(vars,
									    // 0,numEdges);
    }
    
}
