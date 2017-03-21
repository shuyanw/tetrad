package edu.pitt.dbmi.algo.bootstrap;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.ForkJoinPool;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;

//MP: These libraries are required for multi-threading
import edu.cmu.tetrad.util.ForkJoinPoolInstance;
import edu.cmu.tetrad.util.Parameters;

/**
 * Created by Mahdi Pakdaman Naeini on 1/13/17.
 */
public class BootstrapSearch {

    private BootstrapAlgName algName = BootstrapAlgName.RFCI;
    private int numBootstrap = 5;
    private boolean runParallel = false;
    private boolean verbose = false;
    private List<Graph> PAGs = new ArrayList<>();
    private ForkJoinPool pool = null;
    private DataSet data = null;

    private final List<DataSet> bootstrapDatasets = Collections
	    .synchronizedList(new ArrayList<DataSet>());

    private Parameters parameters;
    
    public BootstrapSearch(DataSet data) {
	this.data = data;
	pool = ForkJoinPoolInstance.getInstance().getPool();
    }

    public void addPAG(Graph pag) {
	PAGs.add(pag);
    }

    public DataSet getBootstrapDataset(int id) {
	return bootstrapDatasets.get(id);
    }

    public void addBootstrapDataset(DataSet dataSet) {
	bootstrapDatasets.add(dataSet);
    }

    public void setAlgorithm(BootstrapAlgName algName) {
	this.algName = algName;
    }

    public void setVerbose(boolean verbose) {
	this.verbose = verbose;
    }

    public void setRunningMode(boolean runParallel) {
	this.runParallel = runParallel;
    }

    public void setNumOfBootstrap(int numBootstrap) {
	this.numBootstrap = numBootstrap;
	bootstrapDatasets.clear();
	BootstrapDatasetAction task = new BootstrapDatasetAction(
		numBootstrap, data, this);
	pool.invoke(task);
	while (!pool.isQuiescent()) {
	    try {
		Thread.sleep(100);
	    } catch (InterruptedException e) {
		e.printStackTrace();
	    }
	}
    }

    public Parameters getParameters() {
        return parameters;
    }

    public void setParameters(Parameters parameters) {
        this.parameters = parameters;
    }

    public List<Graph> search() {

	long start, stop;
	if (!this.runParallel) {
	    // Running in the sequential form
	    if (verbose) {
		System.out
			.println("Running Bootstraps in Sequential Mode, numBoostrap = "
				+ numBootstrap);
	    }
	    for (int i1 = 0; i1 < this.numBootstrap; i1++) {
		start = System.currentTimeMillis();

		BootstrapSearchAction task = new BootstrapSearchAction(i1, 1,
			algName, parameters,
			this, true);
		task.compute();

		stop = System.currentTimeMillis();
		if (verbose) {
		    System.out.println("processing time of bootstrap : "
			    + (stop - start) / 1000.0 + " sec");
		}
	    }
	} else {
	    // Running in the parallel multiThread form
	    if (verbose) {
		System.out
			.println("Running Bootstraps in Parallel Mode, numBoostrap = "
				+ numBootstrap);
	    }

	    Parameters parameters = getParameters();

	    BootstrapSearchAction task = new BootstrapSearchAction(0,
		    numBootstrap, algName,
		    parameters, this, true);

	    pool.invoke(task);
	}
	return PAGs;
    }

}
