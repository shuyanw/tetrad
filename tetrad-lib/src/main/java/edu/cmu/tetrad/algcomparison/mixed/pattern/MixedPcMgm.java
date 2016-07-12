package edu.cmu.tetrad.algcomparison.mixed.pattern;

import edu.cmu.tetrad.algcomparison.interfaces.Algorithm;
import edu.cmu.tetrad.algcomparison.interfaces.DataType;
import edu.cmu.tetrad.algcomparison.Parameters;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.search.*;
import edu.pitt.csb.mgm.MGM;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by jdramsey on 6/4/16.
 */
public class MixedPcMgm implements Algorithm {
    public Graph search(DataSet ds, Parameters parameters) {
        MGM m = new MGM(ds, new double[]{
                parameters.getDouble("mgmParam1"),
                parameters.getDouble("mgmParam2"),
                parameters.getDouble("mgmParam3")
        });
        Graph gm = m.search();
        IndependenceTest indTest = new IndTestMixedLrt(ds, parameters.getDouble("alpha"));
        Pc pcs = new Pc(indTest);
        pcs.setDepth(-1);
        pcs.setInitialGraph(gm);
        pcs.setVerbose(false);
        return pcs.search();
    }

    public Graph getComparisonGraph(Graph dag) {
        return SearchGraphUtils.patternForDag(dag);
    }

    public String getDescription() {
        return "PC-Stable, using the Mixed LRT test, with the output of MGM as an intial graph";
    }

    @Override
    public DataType getDataType() {
        return DataType.Mixed;
    }

    @Override
    public List<String> usesParameters() {
        List<String> parameters = new ArrayList<>();
        parameters.add("mgmParam1");
        parameters.add("mgmParam2");
        parameters.add("mgmParam3");
        parameters.add("alpha");
        return parameters;
    }
}
