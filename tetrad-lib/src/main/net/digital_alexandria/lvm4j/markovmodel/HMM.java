/**
 * lvm4j: a Java implementation of various latent variable models.
 * <p>
 * Copyright (C) 2015 - 2016 Simon Dirmeier
 * <p>
 * This file is part of lvm4j.
 * <p>
 * lvm4j is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * <p>
 * lvm4j is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * <p>
 * You should have received a copy of the GNU General Public License
 * along with lvm4j.  If not, see <http://www.gnu.org/licenses/>.
 */


package net.digital_alexandria.lvm4j.markovmodel;

import net.digital_alexandria.lvm4j.DiscreteStateMarkovModel;
import net.digital_alexandria.lvm4j.edges.WeightedArc;
import net.digital_alexandria.lvm4j.util.File;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;


/**
 * Central HMM class, that contains states, transitions etc.
 *
 * @author Simon Dirmeier {@literal s@simon-dirmeier.net}
 */
public final class HMM implements DiscreteStateMarkovModel
{
    // the order of the HMM -> number of previous states that are
    // considered for prediction
    int order;
    // latent variables
    final List<LatentHMMNode<Character, String>> STATES;
    // observed variables
    final List<HMMNode<Character, String>> OBSERVATIONS;
    // arcs between states
    final List<WeightedArc> TRANSITIONS;
    // arcs between states and observations
    final List<WeightedArc> EMISSIONS;
    // flag whether HMM is trained or raw
    boolean isTrained;

    HMM()
    {
        this.STATES = new ArrayList<>();
        this.OBSERVATIONS = new ArrayList<>();
        this.TRANSITIONS = new ArrayList<>();
        this.EMISSIONS = new ArrayList<>();
        this.isTrained = false;
    }

    /**
     * Train the HMM using two files: a file of observations and a file of
     * latent states that emit these observations.
     *
     * @param states a mapping from the id of a state to the
     * real state sequence
     * @param observations a mapping from the id of an observation to the
     * real observations sequence
     */
    @Override
    public void train(Map<String, String> states,
                      Map<String, String> observations)
    {
        HMMTrainer.instance().train(this, states, observations);
        this.isTrained = true;
    }

    /**
     * Predict the most probable latent state sequence using a sequence of
     * observations.  Prediciton is done using the viterbi algorithm.
     *
     * @param y a mapping from the id of an observation to the real
     * observations sequence
     *
     * @return returns a map the predicted states for given sequences
     */
    @Override
    public Map<String, String> predict(Map<String, String> y)
    {
        return HMMPredictor.instance().predict(this, y);
    }

    /**
     * Returns all the log-transformed startProbabilities for the states.
     *
     * @return returns an array
     */
    public final double[] logStartProbabilities()
    {
        double[] probs = new double[this.STATES.size()];
        final double pseudo = 0.000001;
        for (HMMNode<Character, String> s : STATES)
            probs[s.idx()] = Math.log(s.startingProbability() + pseudo);
        return probs;
    }

    /**
     * Returns all the startProbabilities for the states.
     *
     * @return returns an array
     */
    public final double[] startProbabilities()
    {
        double[] probs = new double[this.STATES.size()];
        for (LatentHMMNode<Character, String> s : STATES)
            probs[s.idx()] = s.startingProbability();
        return probs;
    }

    /**
     * Returns the log-transformed stochastic matrix of emissions, i.e. what
     * probabilities do single emissions have
     * for a nodes.
     *
     * @return returns a matrix
     */
    public final double[][] logEmissionMatrix()
    {
        double[][] emissionMatrix =
          new double[this.STATES.size()][this.OBSERVATIONS.size()];
        final double pseudo = 0.000001;
        for (WeightedArc e : EMISSIONS)
            emissionMatrix[e.source().idx()][e.sink().idx()] =
              Math.log(e.weight() + pseudo);
        return emissionMatrix;
    }

    /**
     * Returns the stochastic matrix of emissions, i.e. what probabilities do
     * single emissions have for a nodes.
     *
     * @return returns a matrix
     */
    public final double[][] emissionMatrix()
    {
        double[][] emissionMatrix =
          new double[this.STATES.size()][this.OBSERVATIONS.size()];
        for (WeightedArc e : EMISSIONS)
            emissionMatrix[e.source().idx()][e.sink().idx()] = e.weight();
        return emissionMatrix;
    }

    /**
     * Returns the log-transformed stochastic matrix of transitions,
     * i.e. what probabilties does a transition have given a hidden state.
     *
     * @return returns a matrix
     */
    public final double[][] logTransitionMatrix()
    {
        double[][] transitionsMatrix =
          new double[this.STATES.size()][this.STATES.size()];
        final double pseudo = 0.000001;
        for (WeightedArc t : TRANSITIONS)
            transitionsMatrix[t.source().idx()][t.sink().idx()] =
              Math.log(t.weight() + pseudo);
        return transitionsMatrix;
    }

    /**
     * Returns the stochastic matrix of transitions,
     * i.e. what probabilties does a transition have given a hidden state.
     *
     * @return returns a matrix
     */
    public double[][] transitionMatrix()
    {
        double[][] transitionsMatrix =
          new double[this.STATES.size()][this.STATES.size()];
        for (WeightedArc t : TRANSITIONS)
            transitionsMatrix[t.source().idx()][t.sink().idx()] = t.weight();
        return transitionsMatrix;
    }

    /**
     * Getter for the transitions between the states.
     *
     * @return returns a list of weighted arcs
     */
    public List<WeightedArc> transitions()
    {
        return TRANSITIONS;
    }

    /**
     * Getter for the emissions between hidden states and observations.
     *
     * @return returns a list of weighted arcs
     */
    public List<WeightedArc> emissions()
    {
        return EMISSIONS;
    }

    /**
     * Getter for the hidden states.
     *
     * @return returns a list of nodes
     */
    public List<LatentHMMNode<Character, String>> states()
    {
        return STATES;
    }

    /**
     * Getter of possible obervations.
     *
     * @return returns a list of nodes
     */
    public List<HMMNode<Character, String>> observations()
    {
        return OBSERVATIONS;
    }

    /**
     * Getter for the order of the markov chain.
     *
     * @return returns the order of the markov chain
     */
    public int order() { return order; }

    /**
     * Write the trained HMM parameters to a xml file.
     *
     * @param file the output file
     */
    public void writeHMM(String file)
    {
        File.writeXML(this, file);
    }

    /**
     * Getter for isTrained. True if the HMM has been trained. False if it is
     * the raw HMM.
     *
     * @return returns true if HMM has been trained
     */
    public boolean isTrained() {return isTrained; }

}
