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

package net.digital_alexandria.lvm4j.edges;

import net.digital_alexandria.lvm4j.nodes.Node;

/**
 * Abstract class that stores two nodes as members
 *
 * @author Simon Dirmeier {@literal s@simon-dirmeier.net}
 */
abstract class AbstractArc implements Arc
{
    // starting point of arc
    private final Node _SOURCE;
    // end point of arc
    private final Node _SINK;

    AbstractArc(Node source, Node sink)
    {
        this._SOURCE = source;
        this._SINK = sink;
    }

    @Override
    public String toString()
    {
        return new StringBuilder()
            .append(_SOURCE).append("->")
            .append(_SINK).toString();
    }

    @Override
    public boolean equals(Object o)
    {
        if (o instanceof AbstractArc)
        {
            AbstractArc t = (AbstractArc) o;
            return t._SINK.equals(this._SINK)
                   && t._SOURCE.equals(this._SOURCE);
        }
        return false;
    }

    /**
     * Getter for the end nodes of the arc!
     *
     * @return returns the sink nodes
     */
    public final Node sink()
    {
        return _SINK;
    }

    /**
     * Getter for the start nodes of the arc!
     *
     * @return returns the source nodes
     */
    public final Node source()
    {
        return _SOURCE;
    }
}
