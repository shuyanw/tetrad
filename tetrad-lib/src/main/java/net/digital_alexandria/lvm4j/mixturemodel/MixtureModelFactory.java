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


package net.digital_alexandria.lvm4j.mixturemodel;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * @author Simon Dirmeier {@literal mail@simon-dirmeier.net}
 */
public final class MixtureModelFactory
{
    private MixtureModelFactory() {}

    /**
     * Create a GaussianMixture object with a given matrix that is used
     * for clustering.
     * <p>
     * <code>X</code> is an <code>(n x m)</code> matrix where rows
     * represent samples and columns represent covariates.
     *
     * @param X the matrix for which the mixture is calculated
     *
     * @return returns a GaussianMixture object
     */
    public static GaussianMixture gaussianMixture(double[][] X)
    {
        return new GaussianMixture(X);
    }


    /**
     * Create a GaussianMixture object with a given matrix that is used
     * for clustering.
     * <p>
     * <code>X</code> is an <code>(n x m)</code> matrix where rows
     * represent samples and columns represent covariates.
     *
     * @param X the matrix for which the mixture is calculated
     *
     * @return returns a GaussianMixture object
     */
    public static GaussianMixture gaussianMixture(INDArray X)
    {
        return new GaussianMixture(X);
    }
}
