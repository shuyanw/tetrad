/**
 * lvm4j: a Java implementation of various latent variable models.
 *
 * Copyright (C) 2015 - 2016 Simon Dirmeier
 *
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

import java.util.ArrayList;
import java.util.List;

/**
 * Node class that has a label and a state, a couple of transitions and a couple of emissions.
 * Most importantly it is latend in an HMM.
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 * @param <T> generic for label of node
 * @param <U> generic for state of node
 */
public final class LatentHMMNode<T,U> extends HMMNode <T,U>
{
    // transitions between latent states
    private final List<WeightedArc> _TRANSITIONS;
    // emissions of a state
    private final List<WeightedArc> _EMISSIONS;

    LatentHMMNode(T label, int idx, U state)
    {
        super(label, idx, state);
        this._TRANSITIONS = new ArrayList<>();
        this._EMISSIONS = new ArrayList<>();
    }

    /**
     * Add a transition to the nodes.
     *
     * @param t the transition to be added
     */
    public final void addTransition(WeightedArc t)
    {
        this._TRANSITIONS.add(t);
    }

    /**
     * Add a emission to the nodes.
     *
     * @param e the transition to be added
     */
    public final void addEmission(WeightedArc e)
    {
        this._EMISSIONS.add(e);
    }

    /**
     * Getter for all the emissions the nodes can produce.
     *
     * @return returns a list of arcs to observed variables
     */
    public final List<WeightedArc> emissions()
    {
        return _EMISSIONS;
    }

    /**
     * Getter for all the transitions the nodes can produce.
     *
     * @return returns a list of arcs to hidden variables
     */
    public final List<WeightedArc> transitions()
    {
        return _TRANSITIONS;
    }
}
