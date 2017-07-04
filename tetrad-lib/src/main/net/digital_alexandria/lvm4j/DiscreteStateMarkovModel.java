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

package net.digital_alexandria.lvm4j;

import java.util.Map;

/**
 * Interface for all Markov models
 *
 * @author Simon Dirmeier {@literal mail@simon-dirmeier.net}
 */
public interface DiscreteStateMarkovModel
{
    /**
     * Train the HMM using two files: a file of observations and a file of
     * latent states that emit these observations.
     *
     * @param states a mapping from the id of a state to the
     * real state sequence
     * @param observations a mapping from the id of an observation to the
     * real observations sequence
     */
    public void train(Map<String, String> states,
                      Map<String, String> observations);

    /**
     * Predict the most probable latent state sequence using a sequence of
     * observations.  Prediciton is done using the viterbi algorithm.
     *
     * @param y a mapping from the id of an observation to the real
     * observations sequence
     *
     * @return returns a map the predicted states for given sequences
     */
    public Map<String, String> predict(Map<String, String> y);
}
