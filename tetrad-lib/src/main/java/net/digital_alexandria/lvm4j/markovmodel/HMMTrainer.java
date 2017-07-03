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

import net.digital_alexandria.lvm4j.edges.WeightedArc;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.HashMap;
import java.util.List;
import java.util.Map;


/**
 * Training methods for the HMM.
 *
 * @author Simon Dirmeier {@literal s@simon-dirmeier.net}
 */
final class HMMTrainer
{
    private final static Logger _LOGGER =
      LoggerFactory.getLogger(HMMTrainer.class);
    // singleton
    private static HMMTrainer _trainer;

    static HMMTrainer instance()
    {
        if (_trainer == null)
        {
            _LOGGER.info("Instantiating HMMTrainer");
            _trainer = new HMMTrainer();
        }
        return _trainer;
    }

    private HMMTrainer() {}

    /**
     * Train the HMM using two files: a file of observations and a file of
     * latent states that emit these observations.
     *
     * @param statesMap a mapping from the id of a state to the real state
     * sequence
     * @param observationsMap a mapping from the id of an observation to the
     * real observations sequence
     */
    final void train(HMM hmm,
                     Map<String, String> statesMap,
                     Map<String, String> observationsMap)
    {
        _LOGGER.info("Training HMM!");
        initializeEdges(hmm);
        // a mapping of state -> state object
        final Map<String, LatentHMMNode<Character, String>> labelStatesMap =
          nodeMap(hmm.states());
        // a mapping of state  -> observation label -> emission object
        final Map<String, Map<String, WeightedArc>> labelEmissionsMap =
          edgeMap(hmm.emissions());
        // a mapping of state label -> state label -> transition object
        final Map<String, Map<String, WeightedArc>> labelTransitionsMap =
          edgeMap(hmm.transitions());
        /* count observations, states, emissions and transitions.
         * this is replaced with Baum-Welch algorithm when state sequence is
         * not known
         */
        final int order = hmm.order();
        for (Map.Entry<String, String> statesMapEntry : statesMap.entrySet())
        {
            String id = statesMapEntry.getKey();
            // convert ith state sequence to char array for easy access
            String stateSeq = statesMapEntry.getValue();
            // convert ith observation sequence to char array
            if (!observationsMap.containsKey(id))
                throw new IllegalArgumentException("Observation map does not " +
                                                   "contain:" + id);
            String[] obsArr = observationsMap.get(id).split("");
            if (stateSeq.length() != obsArr.length)
                throw new IllegalArgumentException("Observation sequence and " +
                                                   "state sequence not " +
                                                   "equally long");
            // increase the counter of the state the state sequence begins with.
            for (int i = 0; i < order; i++)
            {
                String statePrefix = stateSeq.substring(0, i + 1);
                if (i == 0)
                    incStartStateCnt(statePrefix, labelStatesMap);
                // increase the counter of the emission of state -> observation
                incEdgeCount(statePrefix, obsArr[i], labelEmissionsMap);
                if (i > 0)
                {
                    String lastStatePrefix = stateSeq.substring(0, i);
                    // increase the counter of the transition of lastState
                    // ->state
                    incEdgeCount(lastStatePrefix, statePrefix,
                                 labelTransitionsMap);
                }
            }
            for (int i = order; i < stateSeq.length(); i++)
            {
                String lastState = stateSeq.substring(i - order, i);
                String curState = stateSeq.substring(i - order + 1, i + 1);
                String curObs = obsArr[i];
                // increase the counter of the transition state -> state
                incEdgeCount(lastState, curState, labelTransitionsMap);
                // increase the counter of the emission state -> observation
                incEdgeCount(curState, curObs, labelEmissionsMap);
            }
        }
        normalize(hmm);
    }

    private void initializeEdges(HMM hmm)
    {
        hmm.transitions().forEach(t -> t.weight(0.0));
        hmm.emissions().forEach(e -> e.weight(0.0));
    }

    private <T extends HMMNode<Character, String>> Map<String, T> nodeMap
      (List<T> l)
    {
        Map<String, T> map = new HashMap<>();
        for (T t : l) map.put(t.state(), t);
        return map;
    }

    @SuppressWarnings("unchecked")
    private <T extends WeightedArc> Map<String, Map<String, T>> edgeMap(
      List<T> l)
    {
        Map<String, Map<String, T>> map = new HashMap<>();
        for (T t : l)
        {
            String source = ((HMMNode<Character, String>) t.source()).state();
            String sink = ((HMMNode<Character, String>) t.sink()).state();
            if (!map.containsKey(source)) map.put(source, new HashMap<>());
            map.get(source).put(sink, t);
        }
        return map;
    }

    private void incStartStateCnt(String state, Map<String,
      LatentHMMNode<Character, String>> map)
    {
        map.get(state).increment();
    }

    private <T extends WeightedArc> void incEdgeCount(
      String source, String sink, Map<String, Map<String, T>> map)
    {
        map.get(source).get(sink).increment();
    }

    private void normalize(HMM hmm)
    {
        double initStateCount = 0;
        for (LatentHMMNode<Character, String> s : hmm.states())
        {
            initStateCount += s.startingProbability();
            double cnt = 0.0;
            for (WeightedArc e : s.emissions())
                cnt += e.weight();
            for (WeightedArc e : s.emissions())
            {
                double p = e.weight() / cnt;
                if (!Double.isFinite(p)) p = 0.0;
                e.weight(p);
            }
            cnt = 0.0;
            for (WeightedArc t : s.transitions())
                cnt += t.weight();
            for (WeightedArc t : s.transitions())
            {
                double p = t.weight() / cnt;
                if (!Double.isFinite(p)) p = 0.0;
                t.weight(p);
            }
        }
        for (HMMNode s : hmm.states())
        {
            double p = s.startingProbability() / initStateCount;
            if (!Double.isFinite(p)) p = 0.0;
            s.startingProbability(p);
        }
    }
}
