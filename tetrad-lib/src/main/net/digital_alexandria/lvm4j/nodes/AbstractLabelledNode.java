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

package net.digital_alexandria.lvm4j.nodes;

/**
 *Abstract class that stores a nodes that as an arbitrary
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 *
 * @param <T> generic for label of node
 */
abstract class AbstractLabelledNode<T> extends AbstractNode
{
    // the label of the nodes
    private final T _LABEL;

    AbstractLabelledNode(T label, int idx)
    {
        super(idx);
        this._LABEL = label;
    }

    /**
     * Getter for the nodes label.
     *
     * @return returns the label
     */
    public T label()
    {
        return this._LABEL;
    }
}
