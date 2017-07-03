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


/**
 * Interface for all classes that produce clusterings.
 *
 * @author Simon Dirmeier {@literal mail@simon-dirmeier.net}
 */
public interface Cluster
{

    /**
     * Create a clustering of k components from a data-set.
     * <p>
     * The number of cluster-center/components is determined by <code>k</code>,
     * i.e. for <code>k=3</code> 3 clusters will be created.
     *
     * @param k the number of components (cluster-centers)
     *
     * @return returns the clustering
     */
    public Clustering cluster(final int k);
}
