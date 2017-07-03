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

import java.util.ArrayList;
import java.util.List;


/**
 * HMMParam class stores several attributes a HMM must have.
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public final class HMMParams
{
    // the possible latent variables of an HMM as characters
    // (i.e. discrete space)
    private char[] _states;
    // the possible observed variables of an HMM as characters (i.e. discrete space)
    private char[] _observations;
    // n-th order of the markov chain
    private int _order;
    // is the provided HMM file trained or not
    private boolean _isTrainingParam;
    // weightings in case the HMM is trained
    private List<Pair<String, Double>> _startWeights;
    // transition probabilities in case the HMM is trained
    private List<Triple<String, String, Double>> _transitionWeights;
    // emission probabilities in case the HMM is trained
    private List<Triple<String, String, Double>> _emissionWeights;

    /**
     * Get a new instance of an HMMParams object.
     *
     * @return returns an HMMParams object
     */
    public static HMMParams newInstance() { return new HMMParams(); }

    HMMParams()
    {
        this._isTrainingParam = false;
        this._startWeights = new ArrayList<>();
        this._transitionWeights = new ArrayList<>();
        this._emissionWeights = new ArrayList<>();
    }

    public final int order() { return _order; }

    public final char[] observations() { return _observations; }

    public final char[] states() { return _states; }

    public final void observations(char[] observations)
    {
        this._observations = observations;
    }

    public final void order(int order) { this._order = order; }

    public final void states(char[] states) { this._states = states; }

    public final List<Pair<String, Double>> startProbabilities()
    {
        return _startWeights;
    }

    public final List<Triple<String, String, Double>> transitionProbabilities()
    {
        return _transitionWeights;
    }

    public final List<Triple<String, String, Double>> emissionProbabilities()
    {
        return _emissionWeights;
    }

    public final void setTrainingParam(boolean b)
    {
        this._isTrainingParam = b;
    }

    public final boolean hasTrainingParams()
    {
        return this._isTrainingParam;
    }
}
