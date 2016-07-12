package edu.cmu.tetrad.algcomparison.mixed.pag;

import edu.cmu.tetrad.algcomparison.interfaces.Algorithm;
import edu.cmu.tetrad.algcomparison.interfaces.DataType;
import edu.cmu.tetrad.algcomparison.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jdramsey on 6/4/16.
 */
public class MixedGfciMixedScore implements Algorithm {
    public Graph search(DataSet dataSet, Parameters parameters) {
        MixedBicScore score = new MixedBicScore(dataSet);
        GFci pc = new GFci(score);
        return pc.search();
    }

    @Override
    public Graph getComparisonGraph(Graph dag) {
        return new DagToPag(dag).convert();
    }

    public String getDescription() {
        return "GFCI using the Mixed BIC score";
    }

    @Override
    public DataType getDataType() {
        return DataType.Mixed;
    }

    @Override
    public List<String> usesParameters() {
        return new ArrayList<>();
    }
}
