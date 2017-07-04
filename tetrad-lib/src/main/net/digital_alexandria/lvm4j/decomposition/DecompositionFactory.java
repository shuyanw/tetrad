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

package net.digital_alexandria.lvm4j.decomposition;

import org.nd4j.linalg.api.ndarray.INDArray;

/**
 * @author Simon Dirmeier {@literal mail@simon-dirmeier.net}
 */
public final class DecompositionFactory
{

    private DecompositionFactory() {}

    /**
     * Create a PCA object with a given matrix that is used for the dimension
     * reduction.
     *
     * @param X the matrix for which the PCA is calculated
     *
     * @return returns an PCA object
     */
    public static PCA pca(double[][] X)
    {
        return new PCA(X);
    }

    /**
     * Create a FactorAnalysis object with a given matrix that is used for
     * creation of a latent space.
     *
     * @param X the matrix for which the FA is calculated
     *
     * @return returns an FA object
     */
    public static FactorAnalysis factorAnalysis(double[][] X)
    {
        return new FactorAnalysis(X);
    }

    /**
     * Create a FactorAnalysis object with a given matrix that is used for
     * creation of a latent space.
     *
     * @param X the matrix for which the FA is calculated
     *
     * @return returns an FA object
     */
    public static FactorAnalysis factorAnalysis(INDArray X)
    {
        return new FactorAnalysis(X);
    }
}
