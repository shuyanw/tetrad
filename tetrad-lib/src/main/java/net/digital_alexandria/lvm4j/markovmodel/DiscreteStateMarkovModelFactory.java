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

import net.digital_alexandria.lvm4j.datastructures.Pair;
import net.digital_alexandria.lvm4j.datastructures.Triple;
import net.digital_alexandria.lvm4j.edges.ArcFactory;
import net.digital_alexandria.lvm4j.edges.WeightedArc;
import net.digital_alexandria.lvm4j.util.File;

import java.util.Collections;
import java.util.List;

import static net.digital_alexandria.lvm4j.util.Combinatorial.combinatorial;


/**
 * HMMFactory class: builds and initializes an HMM
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public final class DiscreteStateMarkovModelFactory
{

    // creator for arcs
    private final static ArcFactory _ARC_FACTORY = ArcFactory.instance();

    /**
     * Create a HMM using the provided file. The HMM can be used for training
     * and prediction. If the edges weights are binary training has to be done
     * at first.
     *
     * @param hmmFile the file containing edges/nodes information
     *
     * @return an HMM
     */
    public static HMM hmm(String hmmFile)
    {
        return HMMbuilder(hmmFile);
    }

    /**
     * Create a HMM using the provided parameters. No training is done.
     *
     * @param states list of chars that represent the states
     * @param observations ist of chars that represent the observations
     * @param order the order of the markov chain
     *
     * @return returns the raw HMM (untrained)
     */
    public static HMM hmm(char states[], char observations[], int order)
    {
        return HMMbuilder(states, observations, order);
    }

    private static HMM HMMbuilder(char[] states, char[] observations, int order)
    {
        HMM hmm = new HMM();
        init(hmm, states, observations, order);
        return hmm;
    }

    private static HMM HMMbuilder(String hmmFile)
    {
        HMM hmm = new HMM();
        init(hmm, hmmFile);
        return hmm;
    }

    private static void init(HMM hmm, String hmmFile)
    {
        // get relevant options of the HMM
        HMMParams params = File.parseXML(hmmFile);
        init(hmm, params);
    }

    private static void init(HMM hmm, HMMParams params)
    {
        // set up nodes
        init(hmm, params.states(), params.observations(), params.order());
        // if the XML provided has trained parameter, initialize a trained HMM
        if (params.hasTrainingParams())
            initTrainingParams(hmm,
                               params.emissionProbabilities(),
                               params.transitionProbabilities(),
                               params.startProbabilities());
    }

    private static void init(HMM hmm, char states[], char observations[], int
      order)
    {
        hmm.order = order;
        // recursively get all combinates of strings of over an alphabet
        // state of size order
        List<String> stateList = combinatorial(states, hmm.order);
        // set up nodes
        init(hmm, stateList, observations);
        // if the XML provided has trained parameter, initialize a trained HMM
    }

    private static void init(HMM hmm, List<String> states, char[] observations)
    {
        addStates(hmm, states);
        addObservations(hmm, observations);
        addTransitions(hmm);
        addEmissions(hmm);
    }

    private static void initTrainingParams(
      HMM hmm,
      List<Triple<String, String, Double>> emissions,
      List<Triple<String, String, Double>> transitions,
      List<Pair<String, Double>> startProbabilities)
    {
        // set up the starting probabilities for every state
        for (Pair<String, Double> p : startProbabilities)
        {
            String state = p.getFirst();
            double prob = p.getSecond();
            hmm.STATES.stream()
                      .filter(s -> s.state().equals(state))
                      .forEach(s -> s.startingProbability(prob));
        }
        // set up the transition probabilities from a state to another state
        transitions.stream().forEach(t -> setUpWeights(t, hmm.TRANSITIONS));
        // set up the emission probabilities from a state to an observation
        emissions.stream().forEach(e -> setUpWeights(e, hmm.EMISSIONS));
    }

    @SuppressWarnings("unchecked")
    private static void setUpWeights(Triple<String, String, Double> t,
                                     List<WeightedArc> it)
    {
        String source = t.getFirst();
        String sink = t.getSecond();
        double prob = t.getThird();
        for (WeightedArc a : it)
        {
            HMMNode<String, String> aso = (HMMNode<String, String>) a.source();
            HMMNode<String, String> asi = (HMMNode<String, String>) a.sink();
            if (aso.state().equals(source) && asi.state().equals(sink))
                a.weight(prob);
        }
    }

    @SuppressWarnings("unchecked")
    private static void addStates(HMM hmm, List<String> states)
    {
        Collections.sort(states, (o1, o2) ->
        {
            if (o1.length() != o2.length())
                return o1.length() < o2.length() ? -1 : 1;
            else
                return o1.compareTo(o2);
        });
        for (int i = 0; i < states.size(); i++)
        {
            String s = states.get(i);
            int length = s.length();
            hmm.STATES.add(
              new LatentHMMNode(s.charAt(length - 1), i, s));
        }
    }

    @SuppressWarnings("unchecked")
    private static void addObservations(HMM hmm, char[] observations)
    {
        for (int i = 0; i < observations.length; i++)
            hmm.OBSERVATIONS.add(
              new HMMNode(observations[i], i, String.valueOf(observations[i])));
    }

    private static void addTransitions(HMM hmm)
    {
        for (int i = 0; i < hmm.STATES.size(); i++)
        {
            LatentHMMNode<Character, String> source = hmm.STATES.get(i);
            hmm.STATES.stream().forEach(sink -> addTransition(hmm, source,
                                                              sink));
        }
    }

    private static void addTransition(HMM hmm,
                                      LatentHMMNode<Character, String> source,
                                      HMMNode<Character, String> sink)
    {
        String sourceSeq = source.state();
        String sinkSeq = sink.state();
        int sourceLength = sourceSeq.length();
        int sinkLength = sinkSeq.length();
        if (sourceLength > sinkLength) return;
        if (sourceLength == sinkLength && sourceLength < hmm.order)
            return;
        if (sourceLength != sinkLength && sourceLength + 1 != sinkLength)
            return;
        String sourceSuffix;
        String sinkPrefix;
        if (sourceLength < hmm.order)
            sourceSuffix = sourceSeq;
        else
            sourceSuffix = sourceSeq.substring(1, sourceLength);
        if (sinkLength < hmm.order)
            sinkPrefix = sinkSeq.substring(0, sourceLength);
        else
            sinkPrefix = sinkSeq.substring(0, sinkLength - 1);
        if (!sourceSuffix.equals(sinkPrefix))
            return;
        WeightedArc t = _ARC_FACTORY.weightedArc(source, sink, .0);
        hmm.TRANSITIONS.add(t);
        source.addTransition(t);
    }

    private static void addEmissions(HMM hmm)
    {
        for (LatentHMMNode<Character, String> state : hmm.STATES)
            hmm.OBSERVATIONS
              .stream()
              .forEach(obs -> addEmission(hmm, state, obs));
    }

    private static void addEmission(HMM hmm, LatentHMMNode source, HMMNode sink)
    {
        WeightedArc t = _ARC_FACTORY.weightedArc(source, sink, .0);
        hmm.EMISSIONS.add(t);
        source.addEmission(t);
    }

}
