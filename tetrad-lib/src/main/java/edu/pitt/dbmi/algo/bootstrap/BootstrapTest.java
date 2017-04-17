package edu.pitt.dbmi.algo.bootstrap;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.IKnowledge;
import edu.cmu.tetrad.data.Knowledge2;
import edu.cmu.tetrad.graph.*;
import edu.cmu.tetrad.util.Parameters;

import java.io.PrintStream;
import java.util.List;

/**
 * Created by mahdi on 1/16/17.
 * 
 * Updated: Chirayu Kong Wongchokprasitti, PhD on 4/5/2017
 * 
 */
public class BootstrapTest {

    private PrintStream out = System.out;

    private final BootstrapSearch bootstrapSearch;

    private Parameters parameters;

    private boolean runParallel = true;

    private BootstrapAlgName algName = BootstrapAlgName.RFCI;
    
    private List<Graph> PAGs;

    private boolean verbose = false;
    
    /**
     * Specification of forbidden and required edges.
     */
    private IKnowledge knowledge = new Knowledge2();

    public void setParallelMode(boolean runParallel) {
	this.runParallel = runParallel;
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
     * @return the output stream that output (except for log output) should be sent to.
     */
    public PrintStream getOut() {
        return out;
    }

    public Parameters getParameters() {
	return parameters;
    }

    public void setParameters(Parameters parameters) {
	this.parameters = parameters;
    }

    /**
     * @return the background knowledge.
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
    
    public BootstrapTest(DataSet data, BootstrapAlgName algName,
	    int numBootstrapSamples) {
	this.algName = algName;
	bootstrapSearch = new BootstrapSearch(data);
	bootstrapSearch.setNumOfBootstrap(numBootstrapSamples);
    }

    public Graph search() {
	long start, stop;

	start = System.currentTimeMillis();

	bootstrapSearch.setAlgorithm(algName);
	bootstrapSearch.setRunningMode(runParallel);
	bootstrapSearch.setVerbose(verbose);
	bootstrapSearch.setParameters(parameters);

	if (verbose) {
	    out.println("Bootstrapping on the " + algName + " algorithm");
	}

	PAGs = bootstrapSearch.search();

	if (verbose) {
	    out.println("Bootstrap size is : " + PAGs.size());
	}
	stop = System.currentTimeMillis();
	if (verbose) {
	    out.println("Processing time of total bootstrapping : "
		    + (stop - start) / 1000.0 + " sec");
	}

	start = System.currentTimeMillis();
	Graph graph = generateBootstrapGraph();
	stop = System.currentTimeMillis();
	if (verbose) {
	    out.println("probDistribution finished in " + (stop - start)
		    + " ms");
	}
	
	return graph;
    }

    private Graph generateBootstrapGraph(){
	Graph complete = new EdgeListGraph(PAGs.get(0).getNodes());
	complete.fullyConnect(Endpoint.TAIL);

	Graph graph = new EdgeListGraph();
	
	for (Edge e : complete.getEdges()) {
		double AnilB = 0.0;
		double AtoB = 0.0;
		double BtoA = 0.0;
		double ACtoB = 0.0;
		double BCtoA = 0.0;
		double AccB = 0.0;
		double AbB = 0.0;
		double AuB = 0.0;

		Node n1 = e.getNode1();
		Node n2 = e.getNode2();
		
		if(!graph.containsNode(n1))graph.addNode(n1);
		if(!graph.containsNode(n2))graph.addNode(n2);
	
		Edge edge = null;
		double maxEdgeProb = 0;
		
		// compute probability for each edge type
		Edge e1 = new Edge(n1, n2, Endpoint.NULL, Endpoint.NULL);
		AnilB = getProbability(e1);

		e1 = new Edge(n1, n2, Endpoint.TAIL, Endpoint.ARROW);
		AtoB = getProbability(e1);
		if(AtoB > maxEdgeProb){
		    edge = e1;
		    maxEdgeProb = AtoB;
		}

		e1 = new Edge(n1, n2, Endpoint.ARROW, Endpoint.TAIL);
		BtoA = getProbability(e1);
		if(BtoA > maxEdgeProb){
		    edge = e1;
		    maxEdgeProb = BtoA;
		}

		e1 = new Edge(n1, n2, Endpoint.CIRCLE, Endpoint.ARROW);
		ACtoB = getProbability(e1);
		if(ACtoB > maxEdgeProb){
		    edge = e1;
		    maxEdgeProb = ACtoB;
		}

		e1 = new Edge(n1, n2, Endpoint.ARROW, Endpoint.CIRCLE);
		BCtoA = getProbability(e1);
		if(BCtoA > maxEdgeProb){
		    edge = e1;
		    maxEdgeProb = BCtoA;
		}

		e1 = new Edge(n1, n2, Endpoint.CIRCLE, Endpoint.CIRCLE);
		AccB = getProbability(e1);
		if(AccB > maxEdgeProb){
		    edge = e1;
		    maxEdgeProb = AccB;
		}

		e1 = new Edge(n1, n2, Endpoint.ARROW, Endpoint.ARROW);
		AbB = getProbability(e1);
		if(AbB > maxEdgeProb){
		    edge = e1;
		    maxEdgeProb = AbB;
		}

		e1 = new Edge(n1, n2, Endpoint.TAIL, Endpoint.TAIL);
		AuB = getProbability(e1);
		if(AuB > maxEdgeProb){
		    edge = e1;
		    maxEdgeProb = AuB;
		}

		if(edge != null){
		    edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeTypeProbability.EdgeType.nil, AnilB));
		    edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeTypeProbability.EdgeType.ta, AtoB));
		    edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeTypeProbability.EdgeType.at, BtoA));
		    edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeTypeProbability.EdgeType.ca, ACtoB));
		    edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeTypeProbability.EdgeType.ac, BCtoA));
		    edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeTypeProbability.EdgeType.cc, AccB));
		    edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeTypeProbability.EdgeType.aa, AbB));
		    edge.addEdgeTypeProbability(new EdgeTypeProbability(EdgeTypeProbability.EdgeType.tt, AuB));
		
		    graph.addEdge(edge);
		}
		
	}
	
	return graph;
    }

    private double getProbability(Edge e) {
	int count = 0;

	for (Graph g : PAGs) {
	    if(!g.containsNode(e.getNode1())){
		throw new IllegalArgumentException();
	    }
	    if(!g.containsNode(e.getNode2())){
		throw new IllegalArgumentException();
	    }
	    
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

}