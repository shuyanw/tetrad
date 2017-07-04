package edu.cmu.tetrad.test.joes.examples;

import edu.cmu.tetrad.algcomparison.algorithm.oracle.pattern.GesMe;
import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.GraphUtils;
import edu.cmu.tetrad.sem.SemIm;
import edu.cmu.tetrad.sem.SemPm;
import edu.cmu.tetrad.util.Parameters;
import org.junit.Test;

/**
 * Created by user on 7/3/17.
 */
public class TestGesMe {

    @Test
    public void test1() {

        Graph graph = GraphUtils.randomGraph(10, 0, 10,
                100, 100, 100, false);
        SemPm pm = new SemPm(graph);
        SemIm im = new SemIm(pm);
        DataSet dataSet = im.simulateData(1000, false);

        GesMe alg = new GesMe();
        Parameters parameters = new Parameters();
        Graph out = alg.search(dataSet, parameters);


    }
}
