package edu.pitt.dbmi.algo.bootstrap;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.RecursiveAction;

import edu.cmu.tetrad.data.CovarianceMatrixOnTheFly;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.BDeuScore;
import edu.cmu.tetrad.search.Fges;
import edu.cmu.tetrad.search.GFci;
import edu.cmu.tetrad.search.IndTestChiSquare;
import edu.cmu.tetrad.search.IndTestFisherZ;
import edu.cmu.tetrad.search.IndependenceTest;
import edu.cmu.tetrad.search.Rfci;
import edu.cmu.tetrad.search.Score;
import edu.cmu.tetrad.search.SemBicScore;
import edu.cmu.tetrad.util.Parameters;

/**
 * 
 * Mar 19, 2017 9:45:44 PM
 * 
 * @author Chirayu (Kong) Wongchokprasitti, PhD
 * 
 */
public class BootstrapSearchAction extends RecursiveAction {

    private static final long serialVersionUID = -5781260555185260539L;

    private int dataSetId;

    private int workLoad;

    private BootstrapAlgName algName;

    private Parameters parameters;

    private final BootstrapSearch bootstrapSearch;

    private boolean verbose;

    public BootstrapSearchAction(int dataSetId, int workLoad,
	    BootstrapAlgName algName, Parameters parameters,
	    BootstrapSearch bootstrapSearch, boolean verbose) {
	this.dataSetId = dataSetId;
	this.workLoad = workLoad;
	this.algName = algName;
	this.parameters = parameters;
	this.bootstrapSearch = bootstrapSearch;
	this.verbose = verbose;
    }

    private Graph learnGraph(DataSet dataSet) {
	Score score = null;
	
	if (dataSet.isContinuous()) {
	    ICovarianceMatrix cov = new CovarianceMatrixOnTheFly(dataSet);
	    SemBicScore semBicScore = new SemBicScore(cov);
	    double penaltyDiscount = 4.0;
	    if (parameters.get("penaltyDiscount") != null) {
		penaltyDiscount = parameters.getDouble("penaltyDiscount");
	    }
	    semBicScore.setPenaltyDiscount(penaltyDiscount);
	    score = semBicScore;
	} else if (dataSet.isDiscrete()) {
	    BDeuScore bDeuScore = new BDeuScore(dataSet);
	    double samplePrior = 1.0;
	    if (parameters.get("samplePrior") != null) {
		samplePrior = parameters.getDouble("samplePrior");
	    }
	    double structurePrior = 1.0;
	    if (parameters.get("structurePrior") != null) {
		structurePrior = parameters.getDouble("structurePrior");
	    }
	    bDeuScore.setSamplePrior(samplePrior);
	    bDeuScore.setStructurePrior(structurePrior);
	    score = bDeuScore;
	}

	IndependenceTest independenceTest = null;
	if (dataSet.isContinuous()) {
	    independenceTest = new IndTestFisherZ(dataSet, parameters.getDouble("alpha", 0.5));
	} else if (dataSet.isDiscrete()) {
	    independenceTest = new IndTestChiSquare(dataSet, parameters.getDouble("alpha", 0.5));
	}

	if (algName == BootstrapAlgName.FGES) {
	    Fges fges = new Fges(score);
	    fges.setFaithfulnessAssumed(parameters.getBoolean("faithfulnessAssumed", true));
	    fges.setNumPatternsToStore(parameters.getInt("numPatternsToStore", 0));
	    return fges.search();
	} else if (algName == BootstrapAlgName.GFCI) {
	    GFci gFci = new GFci(independenceTest, score);
	    gFci.setMaxDegree(parameters.getInt("maxDegree", 5));
	    gFci.setMaxPathLength(parameters.getInt("maxPathLength", -1));
	    gFci.setFaithfulnessAssumed(parameters.getBoolean("faithfulnessAssumed", true));
	    gFci.setCompleteRuleSetUsed(parameters.getBoolean("completeRuleSetUsed", false));
	    return gFci.search();
	} else if (algName == BootstrapAlgName.RFCI) {
	    Rfci rfci = new Rfci(independenceTest);
	    rfci.setCompleteRuleSetUsed(parameters.getBoolean("completeRuleSetUsed", true));
	    rfci.setDepth(parameters.getInt("depth", 3));
	    rfci.setMaxPathLength(parameters.getInt("maxPathLength", -1));
	    return rfci.search();
	} else {
	    throw new IllegalArgumentException(
		    "Bootstrap Search does not support the " + algName
			    + " algorithm yet.");
	}
    }

    @Override
    protected void compute() {
	if (workLoad < 2) {
	    long start, stop;
	    start = System.currentTimeMillis();
	    if (verbose) {
		System.out.println("thread started ... ");
	    }

	    DataSet dataSet = bootstrapSearch.getBootstrapDataset(dataSetId);
	    Graph graph = learnGraph(dataSet);

	    stop = System.currentTimeMillis();
	    if (verbose) {
		System.out
			.println("processing time of bootstrap for thread id : "
				+ dataSetId
				+ " was: "
				+ (stop - start)
				/ 1000.0 + " sec");
	    }
	    bootstrapSearch.addPAG(graph);
	} else {
	    BootstrapSearchAction task1 = new BootstrapSearchAction(dataSetId,
		    workLoad / 2, algName, parameters,
		    bootstrapSearch, verbose);
	    BootstrapSearchAction task2 = new BootstrapSearchAction(dataSetId
		    + workLoad / 2, workLoad - workLoad / 2, algName,
		    parameters, bootstrapSearch, verbose);

	    List<BootstrapSearchAction> tasks = new ArrayList<>();
	    tasks.add(task1);
	    tasks.add(task2);

	    invokeAll(tasks);
	}

    }

}
