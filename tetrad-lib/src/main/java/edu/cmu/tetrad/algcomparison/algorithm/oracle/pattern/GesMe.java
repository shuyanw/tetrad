package edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern;

import com.itemanalysis.psychometrics.factoranalysis.EstimationMethod;
import com.itemanalysis.psychometrics.factoranalysis.ExploratoryFactorAnalysis;
import com.itemanalysis.psychometrics.factoranalysis.RotationMethod;
import edu.cmu.tetrad.algcomparison.algorithm.Algorithm;
import edu.cmu.tetrad.algcomparison.score.ScoreWrapper;
import edu.cmu.tetrad.algcomparison.utils.HasKnowledge;
import edu.cmu.tetrad.algcomparison.utils.TakesInitialGraph;
import edu.cmu.tetrad.data.*;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.search.Fges;
import edu.cmu.tetrad.util.Parameters;
import edu.cmu.tetrad.util.TetradMatrix;
import org.apache.commons.math3.linear.RealMatrix;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * FGES (the heuristic version).
 *
 * @author jdramsey
 */
public class GesMe implements Algorithm, TakesInitialGraph, HasKnowledge {

    static final long serialVersionUID = 23L;
    private boolean compareToTrue = false;
    private ScoreWrapper score;
    private Algorithm initialGraph = null;
    private IKnowledge knowledge = new Knowledge2();

    public GesMe() {
        this.score = score;
        this.compareToTrue = false;
    }

    @Override
    public Graph search(DataModel dataSet, Parameters parameters) {
        Graph initial = null;

        if (initialGraph != null) {
            initial = initialGraph.search(dataSet, parameters);
        }

        int numFactors = parameters.getInt("numFactors");

        int methodIndex = parameters.getInt("estimationMathod");
        EstimationMethod method;

        switch (methodIndex) {
            case 1:
                method = EstimationMethod.GLS;
                break;
            case 2:
                method = EstimationMethod.MINRES;
                break;
            case 3:
                method = EstimationMethod.ML;
                break;
            case 4:
                method = EstimationMethod.PRINCOMP;
                break;
            case 5:
                method = EstimationMethod.WLS;
                break;
            default:
                throw new IllegalArgumentException("Unexpected method index: " + methodIndex);

        }

        int rotationIndex = 1;//parameters.getInt("rotationMethod");
        RotationMethod rotationMethod;

        switch (rotationIndex) {
            case 1:
                rotationMethod = RotationMethod.BENTLER_Q;
                break;
            default:
                throw new IllegalArgumentException("Unrecognized rotation method.");
        }


        RealMatrix cor_r = new CorrelationMatrix((DataSet) dataSet).getMatrix().getRealMatrix();

        ExploratoryFactorAnalysis analysis = new ExploratoryFactorAnalysis(cor_r, numFactors);
        analysis.estimateParameters(method, rotationMethod);

        double[][] _L = new double[((DataSet) dataSet).getNumColumns()][numFactors];

        for (int i = 0; i < ((DataSet) dataSet).getNumColumns(); i++) {
            for (int j = 0; j < numFactors; j++) {
                _L[i][j] = analysis.getFactorMethod().getFactorLoadingAt(i, j);
            }
        }

        TetradMatrix L = new TetradMatrix(_L);

        CovarianceMatrix covariances = new CovarianceMatrix(dataSet.getVariables(), L.times(L.transpose()),
                ((DataSet) dataSet).getNumRows());
        SemBicScoreDeterministic score = new SemBicScoreDeterministic(covariances);
        score.setDeterminismThreshold(parameters.getDouble("determimismThreshold"));
        score.setPenaltyDiscount(parameters.getDouble("penaltyDiscount"));
        edu.cmu.tetrad.search.Fges search = new Fges(score);
        search.setFaithfulnessAssumed(parameters.getBoolean("faithfulnessAssumed"));
        search.setKnowledge(knowledge);
        search.setVerbose(parameters.getBoolean("verbose"));
        search.setMaxDegree(parameters.getInt("maxDegree"));
//        search.setSymmetricFirstStep(parameters.getBoolean("symmetricFirstStep"));

        Object obj = parameters.get("printStream");
        if (obj instanceof PrintStream) {
            search.setOut((PrintStream) obj);
        }

        if (initial != null) {
            search.setInitialGraph(initial);
        }

        return search.search();
    }

    @Override
    public Graph getComparisonGraph(Graph graph) {
        if (compareToTrue) {
            return new EdgeListGraph(graph);
        } else {
            return SearchGraphUtils.patternForDag(new EdgeListGraph(graph));
        }
    }

    @Override
    public String getDescription() {
        return "GES-ME";
    }

    @Override
    public DataType getDataType() {
        return DataType.Continuous;
    }

    @Override
    public List<String> getParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("penaltyDiscount");
        parameters.add("numFactors");
        parameters.add("estimationMathod");
        parameters.add("determimismThreshold");
        return parameters;
    }

    @Override
    public IKnowledge getKnowledge() {
        return knowledge;
    }

    @Override
    public void setKnowledge(IKnowledge knowledge) {
        this.knowledge = knowledge;
    }

    public void setCompareToTrue(boolean compareToTrue) {
        this.compareToTrue = compareToTrue;
    }
}
