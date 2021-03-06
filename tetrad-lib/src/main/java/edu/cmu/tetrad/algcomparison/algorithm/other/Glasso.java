package edu.cmu.tetrad.algcomparison.algorithm.other;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.annotation.AlgType;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.TetradMatrix;
import edu.pitt.dbmi.algo.bootstrap.BootstrapEdgeEnsemble;
import edu.pitt.dbmi.algo.bootstrap.GeneralBootstrapTest;
import java.util.ArrayList;
import java.util.List;

/**
 * @author jdramsey
 */
@edu.cmu.tetrad.annotation.Algorithm(
        name = "GLASSO",
        command = "glasso",
        algoType = AlgType.produce_undirected_graphs
)
public class Glasso implements Algorithm {

    static final long serialVersionUID = 23L;

    public Graph search(DataModel ds, Parameters parameters) {
    	if (parameters.getInt("bootstrapSampleSize") < 1) {
            DoubleMatrix2D cov = new DenseDoubleMatrix2D(DataUtils.getContinuousDataSet(ds)
                    .getCovarianceMatrix().toArray());

            edu.cmu.tetrad.search.Glasso glasso = new edu.cmu.tetrad.search.Glasso(cov);
            glasso.setMaxit((int) parameters.getInt("maxit"));
            glasso.setIa(parameters.getBoolean("ia"));
            glasso.setIs(parameters.getBoolean("is"));
            glasso.setItr(parameters.getBoolean("itr"));
            glasso.setIpen(parameters.getBoolean("ipen"));
            glasso.setThr(parameters.getDouble("thr"));
            glasso.setRhoAllEqual(1.0);

            edu.cmu.tetrad.search.Glasso.Result result = glasso.search();
            TetradMatrix wwi = new TetradMatrix(result.getWwi().toArray());

            List<Node> variables = ds.getVariables();
            Graph resultGraph = new EdgeListGraph(variables);

            for (int i = 0; i < variables.size(); i++) {
                for (int j = i + 1; j < variables.size(); j++) {
                    if (wwi.get(i, j) != 0.0 && wwi.get(i, j) != 0.0) {
                        resultGraph.addUndirectedEdge(variables.get(i), variables.get(j));
                    }
                }
            }

            return resultGraph;
        } else {
            Glasso algorithm = new Glasso();

            DataSet data = (DataSet) ds;
            GeneralBootstrapTest search = new GeneralBootstrapTest(data, algorithm,
                    parameters.getInt("bootstrapSampleSize"));

            BootstrapEdgeEnsemble edgeEnsemble = BootstrapEdgeEnsemble.Highest;
            switch (parameters.getInt("bootstrapEnsemble", 1)) {
                case 0:
                    edgeEnsemble = BootstrapEdgeEnsemble.Preserved;
                    break;
                case 1:
                    edgeEnsemble = BootstrapEdgeEnsemble.Highest;
                    break;
                case 2:
                    edgeEnsemble = BootstrapEdgeEnsemble.Majority;
            }
            search.setEdgeEnsemble(edgeEnsemble);
            search.setParameters(parameters);
            search.setVerbose(parameters.getBoolean("verbose"));
            return search.search();
        }
    }

    public Graph getComparisonGraph(Graph graph) {
        return GraphUtils.undirectedGraph(graph);
    }

    public String getDescription() {
        return "GLASSO (Graphical LASSO)";
    }

    @Override
    public DataType getDataType() {
        return DataType.Mixed;
    }

    @Override
    public List<String> getParameters() {
        List<String> params = new ArrayList<>();
        params.add("maxit");
        params.add("ia");
        params.add("is");
        params.add("itr");
        params.add("ipen");
        params.add("thr");
        // Bootstrapping
        params.add("bootstrapSampleSize");
        params.add("bootstrapEnsemble");
        params.add("verbose");
        return params;
    }
}
