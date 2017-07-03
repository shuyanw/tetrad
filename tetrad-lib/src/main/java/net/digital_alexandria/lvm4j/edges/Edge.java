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
 * Interface for all kinds of edges
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public interface Edge
{
    /**
     * Getter for the two terminal nodes of the edges.
     *
     * @return returns the terminal nodes
     */
    Node[] incidentNodes();
}
