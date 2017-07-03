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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

/**
 * Factory class that produces arcs!
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public class ArcFactory
{
    private final static Logger _LOGGER = LoggerFactory.getLogger(ArcFactory.class);
    // creator of arcs
    private static ArcFactory _factory;

    private ArcFactory(){}

    /**
     * Get an instance of an ArcFactory
     *
     * @return returns an ArcFactory object
     */
    public static ArcFactory instance()
    {
        if (_factory == null)
        {
            _LOGGER.info("Instantiating ArcFactory");
            _factory = new ArcFactory();
        }
        return _factory;
    }

    /**
     * Instantiates an arc that has two nodes as starting and end points and a weight.
     *
     * @param source source nodes of arc
     * @param sink nodes that arc points so
     * @param weight weight of the arc
     * @return returns a new instance of a WeightedArc
     */
    public WeightedArc weightedArc(Node source, Node sink, double weight)
    {
        return new WeightedArc(source, sink, weight);
    }
}
