package edu.pitt.dbmi.algo.bootstrap;

import edu.cmu.tetrad.data.CovarianceMatrixOnTheFly;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.ICovarianceMatrix;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.util.Parameters;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

/**
 * Created by mahdi on 1/16/17.
 */
public class BootstrapExperiment {
    
    private DataSet data;
    private String outputDir = null;
    private String outputBaseFileName = null;
    private String probFileName = null;
    private String probFileNameDirectedPath = null;
    private String truthFileName = null;
    private String truthFileNameDirectedPath = null;
    private Graph trueGraph = null;
    private double penaltyDiscount = 4.0;
    private int maxPathLength = -1;

    private Parameters parameters;
    
    private boolean runParallel = true;
    private int numBootstrapSamples;
    private int depth = 5;
    private double alpha = 0.001;
    private BootstrapAlgName algName = BootstrapAlgName.RFCI;
    private List<Graph> PAGs = new ArrayList<Graph>();
    private Graph estGraph = null;
    private boolean verbose = false;
    private boolean overrideResults = true;

    // Output File handles:
    private PrintWriter outProb = null; // Causality based on direct
					// connectivity
    private PrintWriter outProbDirectedPath = null; // causality based on the
						    // existence of directed
						    // path
    private PrintWriter outTruth = null; // causality based on the existence of
					 // directed path
    private PrintWriter outTruthDirectedPath = null; // causality based on the
						     // existence of directed
						     // path

    public void setOverrideResults(boolean overrideResults) {
	this.overrideResults = overrideResults;
    }

    public void setMaxPathLength(int maxPathLength) {
	this.maxPathLength = maxPathLength;
    }

    public void setPenaltyDiscount(double penaltyDiscount) {
	this.penaltyDiscount = penaltyDiscount;
    }

    public void setParallelMode(boolean runParallel) {
	this.runParallel = runParallel;
    }

    public void setTrueGraph(Graph trueGraph) {
	this.trueGraph = trueGraph;
    }

    public void setVerbose(boolean verbose) {
	this.verbose = verbose;
    }

    public Parameters getParameters() {
        return parameters;
    }

    public void setParameters(Parameters parameters) {
        this.parameters = parameters;
    }

    public BootstrapExperiment(DataSet data, BootstrapAlgName algName,
	    int numBootstrapSamples) {
	this.data = data;
	this.algName = algName;
	this.numBootstrapSamples = numBootstrapSamples;
    }

    public void runExperiment() {
	long start, stop;
	if (checkProbFileExists() == -1) {
	    return;
	}
	BootstrapSearch bootstrapSearchObj = new BootstrapSearch(data);
	start = System.currentTimeMillis();
	bootstrapSearchObj.setAlgorithm(algName);
	bootstrapSearchObj.setRunningMode(runParallel);
	bootstrapSearchObj.setVerbose(verbose);
	bootstrapSearchObj.setNumOfBootstrap(numBootstrapSamples);
	bootstrapSearchObj.setParameters(parameters);
	estGraph = learnGraph(data);
	if (trueGraph != null) {
	    estGraph = GraphUtils.replaceNodes(estGraph, trueGraph.getNodes());
	}

	PAGs = bootstrapSearchObj.search();

	if (trueGraph != null) {
	    List<Graph> PAGs2 = new ArrayList<Graph>();
	    for (Graph pag : PAGs) {
		pag = GraphUtils.replaceNodes(pag, trueGraph.getNodes());
		PAGs2.add(pag);
	    }
	    PAGs = PAGs2;
	}

	if (verbose) {
	    System.out.println("Bootstrap size is : " + PAGs.size());
	}
	stop = System.currentTimeMillis();
	if (verbose) {
	    System.out.println("processing time of total bootstrapping : "
		    + (stop - start) / 1000.0 + " sec");
	}

	start = System.currentTimeMillis();
	writeProbFile();
	stop = System.currentTimeMillis();
	if (verbose) {
	    System.out.println("probDistribution finished in " + (stop - start)
		    + " ms");
	}

	start = System.currentTimeMillis();
	writeProbFileBasedOnDirectedPath();
	stop = System.currentTimeMillis();
	if (verbose) {
	    System.out
		    .println("probDistribution based on Directed Path finished in "
			    + (stop - start) + " ms");
	}

	if (trueGraph != null) {
	    start = System.currentTimeMillis();
	    writeTruth(trueGraph);
	    stop = System.currentTimeMillis();
	    if (verbose) {
		System.out.println("Writing true graph finished in "
			+ (stop - start) + " ms");
	    }

	    start = System.currentTimeMillis();
	    writeTruthBasedOnDirectedPath(trueGraph);
	    stop = System.currentTimeMillis();
	    if (verbose) {
		System.out
			.println("Writing true graph based on Directed Path finished in "
				+ (stop - start) + " ms");
	    }

	}

	closeOutputFiles();
    }

    private int writeTruth(Graph trueGraph) {
	File truthFile = new File(outputDir, truthFileName);
	try {
	    outTruth = new PrintWriter(truthFile);
	} catch (FileNotFoundException e) {
	    e.printStackTrace();
	    return -1;
	}
	// Write Header:
	print("A, B, truth={0-7}" + " \n", outTruth);

	Graph complete = new EdgeListGraph(trueGraph.getNodes());
	complete.fullyConnect(Endpoint.TAIL);

	for (Edge e : complete.getEdges()) {
	    int trueType = 0;

	    Node n1 = e.getNode1();
	    Node n2 = e.getNode2();

	    // compute true edge type for each pair of nodes

	    if (trueGraph.getEdge(n1, n2) == null) {
		trueType = 0;
	    } else {
		Endpoint p1 = trueGraph.getEdge(n1, n2).getEndpoint1();
		Endpoint p2 = trueGraph.getEdge(n1, n2).getEndpoint2();

		if (p1 == Endpoint.TAIL && p2 == Endpoint.ARROW) // A -> B
		    trueType = 1;

		else if (p1 == Endpoint.ARROW && p2 == Endpoint.TAIL) // A <- B
		    trueType = 2;

		else if (p1 == Endpoint.CIRCLE && p2 == Endpoint.ARROW) // A o->
									// B
		    trueType = 3;

		else if (p1 == Endpoint.ARROW && p2 == Endpoint.CIRCLE) // A <-o
									// B
		    trueType = 4;

		else if (p1 == Endpoint.CIRCLE && p2 == Endpoint.CIRCLE) // A
									 // o-o
									 // B
		    trueType = 5;

		else if (p1 == Endpoint.ARROW && p2 == Endpoint.ARROW) // A <->
								       // B
		    trueType = 6;

		else if (p1 == Endpoint.TAIL && p2 == Endpoint.TAIL) // A -- B
		    trueType = 7;

	    }

	    print(n1 + ", " + n2 + ", " + trueType + "\n", outTruth);
	}
	return 0;
    }

    private int writeTruthBasedOnDirectedPath(Graph trueGraph) {
	File truthFileDirectedPath = new File(outputDir,
		truthFileNameDirectedPath);
	try {
	    outTruthDirectedPath = new PrintWriter(truthFileDirectedPath);
	} catch (FileNotFoundException e) {
	    e.printStackTrace();
	    return -1;
	}
	// Write Header:
	print("A, B, truth={0-3}" + " \n", outTruthDirectedPath);
	Graph complete = new EdgeListGraph(trueGraph.getNodes());
	complete.fullyConnect(Endpoint.TAIL);

	for (Edge e : complete.getEdges()) {
	    int trueType = 0;

	    Node n1 = e.getNode1();
	    Node n2 = e.getNode2();

	    // compute true edge type for each pair of nodes
	    if (trueGraph.existsDirectedPathFromTo(n1, n2)) {
		if (trueGraph.existsDirectedPathFromTo(n2, n1)) {
		    // There is a bidirectional causal relation between two
		    // nodes according to the true PAG
		    trueType = 3;
		} else {
		    trueType = 1;
		}
	    } else if (trueGraph.existsDirectedPathFromTo(n2, n1)) {
		trueType = 2;
	    } else {
		// There is no causal relation between the two nodes based on
		// the true PAG
		trueType = 0;
	    }
	    print(n1 + ", " + n2 + ", " + trueType + "\n", outTruthDirectedPath);
	}
	return 0;
    }

    private int writeTruthBasedOnDirectedPath(HashMap<String, HashMap> trueGraph) {
	// An axuiliary function to write the output file of true causal
	// relationships for the yeast dataset
	// trueGraph is a hashMap that records the adjancyList of true causal
	// relationships

	File truthFileDirectedPath = new File(outputDir,
		truthFileNameDirectedPath);
	try {
	    outTruthDirectedPath = new PrintWriter(truthFileDirectedPath);
	} catch (FileNotFoundException e) {
	    e.printStackTrace();
	    return -1;
	}
	// Write Header:
	print("A, B, truth={1-3}" + " \n", outTruthDirectedPath);
	System.out.println("MP: on top of complete Graph!");
	Graph complete = new EdgeListGraph(estGraph.getNodes());
	complete.fullyConnect(Endpoint.TAIL);
	System.out.println("MP: Start Edging!");
	for (Edge e : complete.getEdges()) {
	    int trueType = 0;
	    Node n1 = e.getNode1();
	    Node n2 = e.getNode2();

	    // compute true edge type for each pair of nodes
	    String n1Name = n1.getName();
	    String n2Name = n2.getName();
	    // Check if the n1->n2 is not among true causal relationships
	    // if n1->n2 is not among true causal relationships then continue

	    if (trueGraph.containsKey(n1Name)
		    && trueGraph.get(n1Name).containsKey(n2Name)) {

		if (trueGraph.containsKey(n2Name)
			&& trueGraph.get(n2Name).containsKey(n1Name)) {
		    // System.out.println(n1Name+" , "+ n2Name
		    // +" bidirected causal relation!");
		    trueType = 3; // This means there is a bidirectional true
				  // causal relationship between two nodes
		    // according to the interventional data in yeast!
		    trueGraph.get(n2Name).put(n1Name,
			    (int) trueGraph.get(n2Name).get(n1Name) + 1);
		} else {
		    trueType = 1;
		}
		trueGraph.get(n1Name).put(n2Name,
			(int) trueGraph.get(n1Name).get(n2Name) + 1);
	    } else if (trueGraph.containsKey(n2Name)
		    && trueGraph.get(n2Name).containsKey(n1Name)) {
		trueType = 2;
		trueGraph.get(n2Name).put(n1Name,
			(int) trueGraph.get(n2Name).get(n1Name) + 1);
	    } else {
		continue;
	    }
	    print(n1 + ", " + n2 + ", " + trueType + "\n", outTruthDirectedPath);
	}
	return 0;
    }

    private int writeProbFile() {
	File probFile = new File(outputDir, probFileName);
	try {
	    outProb = new PrintWriter(probFile);
	} catch (FileNotFoundException e) {
	    e.printStackTrace();
	    return -1;
	}
	// Write Header
	print("A, B, freq(A-nil-B), freq(A --> B), freq(B --> A), freq(A o-> B), freq(B o-> A), "
		+ "freq(A o-o B), freq(A <-> B), freq(A --- B), "
		+ algName
		+ "_out={0-7}" + " \n", outProb);

	System.out.println("MP: on top of complete Graph!");
	Graph complete = new EdgeListGraph(estGraph.getNodes());
	complete.fullyConnect(Endpoint.TAIL);
	System.out.println("MP: Start Edging!");

	for (Edge e : complete.getEdges()) {
	    double AnilB = 0.0;
	    double AtoB = 0.0;
	    double BtoA = 0.0;
	    double ACtoB = 0.0;
	    double BCtoA = 0.0;
	    double AccB = 0.0;
	    double AbB = 0.0;
	    double AuB = 0.0;
	    int estType = 0;

	    Node n1 = e.getNode1();
	    Node n2 = e.getNode2();

	    // compute probability for each edge type
	    Edge e1 = new Edge(n1, n2, Endpoint.NULL, Endpoint.NULL);
	    AnilB = getProbability(e1);

	    e1 = new Edge(n1, n2, Endpoint.TAIL, Endpoint.ARROW);
	    AtoB = getProbability(e1);

	    e1 = new Edge(n1, n2, Endpoint.ARROW, Endpoint.TAIL);
	    BtoA = getProbability(e1);

	    e1 = new Edge(n1, n2, Endpoint.CIRCLE, Endpoint.ARROW);
	    ACtoB = getProbability(e1);

	    e1 = new Edge(n1, n2, Endpoint.ARROW, Endpoint.CIRCLE);
	    BCtoA = getProbability(e1);

	    e1 = new Edge(n1, n2, Endpoint.CIRCLE, Endpoint.CIRCLE);
	    AccB = getProbability(e1);

	    e1 = new Edge(n1, n2, Endpoint.ARROW, Endpoint.ARROW);
	    AbB = getProbability(e1);

	    e1 = new Edge(n1, n2, Endpoint.TAIL, Endpoint.TAIL);
	    AuB = getProbability(e1);

	    if (estGraph.getEdge(n1, n2) == null)
		estType = 0;

	    else {

		Endpoint p1 = estGraph.getEdge(n1, n2).getEndpoint1();
		Endpoint p2 = estGraph.getEdge(n1, n2).getEndpoint2();

		if (p1 == Endpoint.TAIL && p2 == Endpoint.ARROW) // A -> B
		    estType = 1;

		else if (p1 == Endpoint.ARROW && p2 == Endpoint.TAIL) // A <- B
		    estType = 2;

		else if (p1 == Endpoint.CIRCLE && p2 == Endpoint.ARROW) // A o->
									// B
		    estType = 3;

		else if (p1 == Endpoint.ARROW && p2 == Endpoint.CIRCLE) // A <-o
									// B
		    estType = 4;

		else if (p1 == Endpoint.CIRCLE && p2 == Endpoint.CIRCLE) // A
									 // o-o
									 // B
		    estType = 5;

		else if (p1 == Endpoint.ARROW && p2 == Endpoint.ARROW) // A <->
								       // B
		    estType = 6;

		else if (p1 == Endpoint.TAIL && p2 == Endpoint.TAIL) // A -- B
		    estType = 7;
	    }
	    print(n1 + ", " + n2 + ", " + AnilB + ", " + AtoB + ", " + BtoA
		    + ", " + ACtoB + ", " + BCtoA + ", " + AccB + ", " + AbB
		    + ", " + AuB + ", " + estType + "\n", outProb);
	}
	return 0;
    }

    private int writeProbFileBasedOnDirectedPath() {
	File probFileDirectedPath = new File(outputDir,
		probFileNameDirectedPath);
	try {
	    outProbDirectedPath = new PrintWriter(probFileDirectedPath);
	} catch (FileNotFoundException e) {
	    e.printStackTrace();
	    return -1;
	}
	// Write Header
	print("A, B, freq(A-->B), freq(B-->A), freq(A<->B), " + algName
		+ " \n", outProbDirectedPath);
	if (verbose) {
	    System.out.println("MP: on top of complete Graph!");
	}
	Graph complete = new EdgeListGraph(estGraph.getNodes());
	complete.fullyConnect(Endpoint.TAIL);
	if (verbose) {
	    System.out.println("MP: Start Edging!");
	}

	for (Edge e : complete.getEdges()) {
	    double AtoB = 0.0;
	    double BtoA = 0.0;
	    double AbB = 0.0;

	    Node n1 = e.getNode1();
	    Node n2 = e.getNode2();

	    // compute probability for each causal edge type
	    Edge e1;
	    e1 = new Edge(n1, n2, Endpoint.TAIL, Endpoint.ARROW);
	    AtoB = getProbabilityBasedOnDirectedPath(e1);

	    e1 = new Edge(n1, n2, Endpoint.ARROW, Endpoint.TAIL);
	    BtoA = getProbabilityBasedOnDirectedPath(e1);

	    e1 = new Edge(n1, n2, Endpoint.ARROW, Endpoint.ARROW);
	    AbB = getProbabilityBasedOnDirectedPath(e1);

	    // MP: Starting of the part that should be edited
	    int estType = 0;
	    if (estGraph.existsDirectedPathFromTo(n1, n2)) {
		if (estGraph.existsDirectedPathFromTo(n2, n1)) {
		    // There is a bidirectional causal relation between two
		    // nodes according to the estimated Graph
		    estType = 3;
		} else {
		    estType = 1;
		}
	    } else if (estGraph.existsDirectedPathFromTo(n2, n1)) {
		estType = 2;
	    } else {
		// There is no causal relation between two nodes based on
		// estimated Graph (i.e., generated by algorithm e.g., FGES,
		// RFCI, etc.)
		estType = -1;
	    }
	    print(n1 + ", " + n2 + ", " + AtoB + ", " + BtoA + ", " + AbB
		    + ", " + estType + "\n", outProbDirectedPath);
	}
	return 0;
    }

    private void closeOutputFiles() {
	if (this.outProb != null) {
	    this.outProb.close();
	}
	if (this.outProbDirectedPath != null) {
	    this.outProbDirectedPath.close();
	}
	if (this.outTruth != null) {
	    this.outTruth.close();
	}
	if (this.outTruthDirectedPath != null) {
	    this.outTruthDirectedPath.close();
	}
    }

    private void print(String s, PrintWriter out) {
	if (out == null)
	    return;
	out.flush();
	out.print(s);
	out.flush();
    }

    public void setOut(String dirName, String baseFileName) {
	// Here just create file names and store them in the instance variables
	// Ad check availability to check if prob file exist then ignore
	// running!

	File dir = new File(dirName);
	dir.mkdirs();
	outputBaseFileName = baseFileName;
	outputDir = dirName;

	// probFileName: In this file causality is defined based on direct
	// connectivity
	probFileName = "probs_" + baseFileName + ".csv";

	// probFileName_directedPath: In this file causality is defined based on
	// the existence of directed path
	probFileNameDirectedPath = "probs_" + baseFileName + "_directedPath"
		+ ".csv";

	// In this file causality is defined based on the existence of direct
	// connectivity(edge) in the true PAG
	truthFileName = "truth_" + baseFileName + ".csv";

	// In this file causality is defined based on the existence of directed
	// path in the true PAG
	truthFileNameDirectedPath = "truth_" + baseFileName + "_directedPath"
		+ ".csv";
    }

    public int checkProbFileExists() {
	File dir = new File(outputDir);
	dir.mkdirs();
	File probFile = new File(dir, probFileName);
	if (probFile.exists()) {
	    if (overrideResults) {
		if (verbose) {
		    System.out
			    .println("warning: probFile already exists and it will be override!");
		}
	    } else {
		System.out
			.println("Erorr: probFile already exists and it cannot be override until you turn on the overrideResults flag!");
		return -1;
	    }
	}
	return 0;
    }

    private Graph learnGraph(DataSet dataSet) {
	Graph res = null;
	if (algName == BootstrapAlgName.RFCI) {
	    res = learnRFCI(dataSet);
	} else if (algName == BootstrapAlgName.FGES) {
	    res = learnFGES(dataSet);
	} else if (algName == BootstrapAlgName.GFCI) {
	    res = learnGFci(dataSet);
	} else {
	    res = null;
	}

	return res;
    }

    private Graph learnGFci(DataSet dataSet) {
	ICovarianceMatrix cov = new CovarianceMatrixOnTheFly(dataSet);
	IndependenceTest independenceTest = new IndTestFisherZ(cov, alpha);
	// GFci fci = new GFci(independenceTest);

	SemBicScore score = new SemBicScore(cov);
	score.setPenaltyDiscount(penaltyDiscount);
	GFci fci = new GFci(independenceTest, score);

	fci.setVerbose(false);
	fci.setMaxPathLength(maxPathLength);
	fci.setMaxDegree(depth);
	fci.setFaithfulnessAssumed(false);
	fci.setCompleteRuleSetUsed(true);
	Graph outGraph = fci.search();
	return outGraph;
    }

    private Graph learnFGES(DataSet dataSet) {
	double structurePrior = 1, samplePrior = 1;
	BDeuScore score = new BDeuScore(dataSet);
	score.setSamplePrior(samplePrior);
	score.setStructurePrior(structurePrior);
	Fges ges = new Fges(score);
	ges.setVerbose(false);
	ges.setNumPatternsToStore(0);
	ges.setFaithfulnessAssumed(true);
	// ges.setDepth(depth);
	ges.setCycleBound(-1);
	Graph estPag = ges.search();
	return estPag;
    }

    private Graph learnRFCI(DataSet dataSet) {
	ICovarianceMatrix cov = new CovarianceMatrixOnTheFly(dataSet);
	final IndTestFisherZ test = new IndTestFisherZ(cov, alpha);
	// For discrete variables
	// final IndTestChiSquare test = new IndTestChiSquare(dataSet, alpha);
	Rfci fci = new Rfci(test);
	fci.setVerbose(false);
	fci.setCompleteRuleSetUsed(true);
	fci.setDepth(depth);
	Graph estPag = fci.search();
	return estPag;
    }

    private double getProbability(Edge e) {
	int count = 0;

	if (!PAGs.get(0).containsNode(e.getNode1()))
	    throw new IllegalArgumentException();
	if (!PAGs.get(0).containsNode(e.getNode2()))
	    throw new IllegalArgumentException();

	for (Graph g : PAGs) {
	    if (e.getEndpoint1() == Endpoint.NULL
		    || e.getEndpoint2() == Endpoint.NULL) {
		if (!g.isAdjacentTo(e.getNode1(), e.getNode2()))
		    count++;
	    } else {
		if (g.containsEdge(e))
		    count++;
	    }
	}

	return count / (double) PAGs.size();
    }

    private double getProbabilityBasedOnDirectedPath(Edge e) {
	int count = 0;

	if (!PAGs.get(0).containsNode(e.getNode1()))
	    throw new IllegalArgumentException();
	if (!PAGs.get(0).containsNode(e.getNode2()))
	    throw new IllegalArgumentException();

	for (Graph g : PAGs) {
	    if (e.getEndpoint1() == Endpoint.TAIL
		    && e.getEndpoint2() == Endpoint.ARROW) {
		if (g.existsDirectedPathFromTo(e.getNode1(), e.getNode2()))
		    count++;

	    } else if (e.getEndpoint1() == Endpoint.ARROW
		    && e.getEndpoint2() == Endpoint.TAIL) {
		if (g.existsDirectedPathFromTo(e.getNode2(), e.getNode1()))
		    count++;

	    } else if (e.getEndpoint1() == Endpoint.ARROW
		    && e.getEndpoint2() == Endpoint.ARROW) {
		if (g.existsDirectedPathFromTo(e.getNode1(), e.getNode2())
			&& g.existsDirectedPathFromTo(e.getNode2(),
				e.getNode1()))
		    count++;

	    } else if (e.getEndpoint1() == Endpoint.NULL
		    || e.getEndpoint2() == Endpoint.NULL) {
		if (!g.isAdjacentTo(e.getNode1(), e.getNode2()))
		    count++;

	    } else {
		if (g.containsEdge(e))
		    count++;
	    }
	}
	return count / (double) PAGs.size();
    }

}