package edu.pitt.dbmi.algo.bootstrap;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;

import java.io.BufferedReader;
import java.io.FileReader;

import edu.cmu.tetrad.data.DelimiterType;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import edu.cmu.tetrad.data.DataSet;
import edu.cmu.tetrad.data.DataUtils;
import edu.cmu.tetrad.graph.Edge;
import edu.cmu.tetrad.graph.EdgeListGraph;
import edu.cmu.tetrad.graph.Endpoint;
import edu.cmu.tetrad.graph.Graph;
import edu.cmu.tetrad.graph.Node;
import edu.cmu.tetrad.graph.NodeEqualityMode;


//MP: These libraries are required for multi-threading
import java.lang.Thread;
import java.lang.Runtime;
import edu.cmu.tetrad.data.DataReader;
import edu.cmu.tetrad.data.CovarianceMatrixOnTheFly;
import edu.cmu.tetrad.data.ICovarianceMatrix;

import java.util.concurrent.*;

import edu.cmu.tetrad.search.*;
import edu.cmu.tetrad.util.ForkJoinPoolInstance;

/**
 * Created by Mahdi Pakdaman Naeini on 1/13/17.
 */
public class BootstrapSearch {
    private  String algorithm = "rfci";
    private  int numBootstrap = 5;
    private  boolean runParallel = false;
    private  double alpha = 0.001;
    private  int depth = 5;
    private  boolean verbose = false;
    private  List<Graph> PAGs = new ArrayList<Graph>();
    private  ForkJoinPool pool = null;
    private  DataSet data = null;
    private double penaltyDiscount = 4.0;
    private int maxPathLength = -1;

    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }
    public void setAlgorithm(String algorithm) {
        this.algorithm = algorithm;
    }
    public void setDepth(int depth) {
        this.depth = depth;
    }
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }
    public void setRunningMode(boolean runParallel) {
        this.runParallel = runParallel;
    }
    public void setNumOfBootstrap(int numBootstrap) {
        this.numBootstrap = numBootstrap;
    }
    public void setMaxPathLength(int maxPathLength){this.maxPathLength = maxPathLength;}
    public void setPenaltyDiscount(double penaltyDiscount){this.penaltyDiscount = penaltyDiscount;}


    public BootstrapSearch(DataSet data){
        this.data = data;
        pool = ForkJoinPoolInstance.getInstance().getPool();
    }

    public List<Graph> search(){
        class ParallelBootstrapSearch extends RecursiveAction{
            private int count = -1;
//            private DataSet bootstrapSample = null;

            public ParallelBootstrapSearch(int count){
                this.count = count;
            }

            private Graph learnGraph(DataSet dataSet){
                if (algorithm.equals("rfci")){
                    return learnRFCI(dataSet);
                }else if (algorithm.equals("fges")){
                    return learnFGES(dataSet);
                }else if (algorithm.equals("gfci")){
                     return learnGFci(dataSet);
                }else{
                    return null;
                }
            }

            private Graph learnGFci(DataSet dataSet){
                ICovarianceMatrix cov = new CovarianceMatrixOnTheFly(dataSet);
                IndependenceTest independenceTest = new IndTestFisherZ(cov, alpha);
//        GFci fci = new GFci(independenceTest);

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

            private Graph learnFGES(DataSet dataSet){
                double structurePrior = 1, samplePrior = 1;
                BDeuScore score = new BDeuScore(dataSet);
                score.setSamplePrior(samplePrior);
                score.setStructurePrior(structurePrior);
                Fges ges = new Fges(score);
                ges.setVerbose(false);
                ges.setNumPatternsToStore(0);
                ges.setFaithfulnessAssumed(true);
//			ges.setDepth(depth);
                ges.setCycleBound(-1);
                Graph estPag = ges.search();
                return estPag;

            }

            private Graph learnRFCI(DataSet dataSet) {
                ICovarianceMatrix cov = new CovarianceMatrixOnTheFly(dataSet);
                final IndTestFisherZ test = new IndTestFisherZ(cov, alpha);
                //		For discrete variables
                //		final IndTestChiSquare test = new IndTestChiSquare(dataSet, alpha);

                Rfci fci = new Rfci(test);
                fci.setVerbose(false);
                fci.setCompleteRuleSetUsed(true);
                fci.setDepth(depth);
                Graph estPag = fci.search();
                return estPag;
            }

            @Override
            protected void compute() {
                ParallelBootstrapSearch task1 = null, task2 = null;
                if (this.count<2){
                    DataSet bootstrapData_ = DataUtils.getBootstrapSample(data, data.getNumRows());
                    long start, stop;
                    start =  System.currentTimeMillis();
                    if (verbose) {
                        System.out.println("thread started ... ");
                    }
                    Graph outGraph = this.learnGraph(bootstrapData_);
                    stop =  System.currentTimeMillis();
                    if (verbose){
                        System.out.println("processing time of bootstrap for thread id : " +this.count+" was: "+ (stop-start)/1000.0 +" sec");
                    }
                    PAGs.add(outGraph);
                }else{
                    task1 = new ParallelBootstrapSearch(this.count/2);
                    task2 = new ParallelBootstrapSearch(this.count - this.count/2);

                    List<ParallelBootstrapSearch> tasks = new ArrayList<>();
                    tasks.add(task1);
                    tasks.add(task2);
                    invokeAll(tasks);
                }

            }
        }

        long start, stop;
        if (!this.runParallel){
            //Running in the sequential form
            if (verbose){
                System.out.println("Running Bootstraps in Sequential Mode, numBoostrap = "+numBootstrap);
            }
            for (int i1 = 0; i1 < this.numBootstrap; i1++) {
//                DataSet bootstrapData = DataUtils.getBootstrapSample(data, data.getNumRows());
                start =  System.currentTimeMillis();
                ParallelBootstrapSearch aSingleTask = new ParallelBootstrapSearch(1);
                aSingleTask.compute();
                stop =  System.currentTimeMillis();
                if (verbose){
                    System.out.println("processing time of bootstrap : " + (stop-start)/1000.0 +" sec");
                }
            }
        }else{
            //Running in the parallel multiThread form
            if (verbose){
                System.out.println("Running Bootstraps in Parallel Mode, numBoostrap = "+numBootstrap);
            }
//            DataSet bootstrapData = DataUtils.getBootstrapSample(data, data.getNumRows());
            ParallelBootstrapSearch task = new ParallelBootstrapSearch(numBootstrap);
            pool.invoke(task);
        }
        return PAGs;
    }
}

// Trash Codes:

//    protected void compute1() {
//        ParallelBootstrapSearch task = null;
//        if (this.count>1) {
////                    DataSet bootstrapData_ = DataUtils.getBootstrapSample(data, data.getNumRows());
//            task = new ParallelBootstrapSearch( this.count - 1);
//            task.fork();
//        }
//        long start, stop;
//        start =  System.currentTimeMillis();
//        System.out.println("thread started booboo ... ");
//        Graph outGraph = this.learnGraph(this.bootstrapSample);
//        stop =  System.currentTimeMillis();
//        if (verbose){
//            System.out.println("processing time of bootstrap for thread id : " +this.count+" was: "+ (stop-start)/1000.0 +" sec");
//        }
//        if (this.count > 1) {
//            task.join();
//        }
//        PAGs.add(outGraph);
//    }
