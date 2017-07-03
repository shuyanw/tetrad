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

import net.digital_alexandria.lvm4j.nodes.LabelledNode;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

/**
 * Node class that has a label and a state, a couple of transitions and a
 * couple of emissions.
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 *
 * @param <T> generic for label of node
 * @param <U> generic for state of node
 */
public class HMMNode<T, U> extends LabelledNode<T>
{
    // the state of the HMMNode
    private final U _STATE;
    // the startng probability
    private double _probStart;

    HMMNode(T label, int idx, U state)
    {
        super(label, idx);
        this._STATE = state;
        _probStart = 0.0;
    }

    @Override
    public String toString()
    {
        return new StringBuilder().append(label()).append("-").append(_STATE)
                                  .toString();
    }

    @Override
    public boolean equals(Object o)
    {
        throw new NotImplementedException();
    }

    /**
     * Increments starting probability by one!
     */
    public final void increment()
    {
        this._probStart++;
    }

    /**
     * Getter for starting probability!
     *
     * @return returns the starting probability
     */
    public final double startingProbability()
    {
        return _probStart;
    }

    /**
     * Setter for starting probability!
     *
     * @param d the probability do set
     */
    public final void startingProbability(double d)
    {
        this._probStart = d;
    }

    /**
     * Getter for the state of the nodes.
     *
     * @return returns the state.
     */
    public final U state() { return _STATE; }
}
