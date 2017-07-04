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


package net.digital_alexandria.lvm4j.util;

import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;
import org.nd4j.linalg.api.ndarray.INDArray;

import static java.lang.Math.sqrt;


/**
 *
 * @author Simon Dirmeier {@literal mail@simon-dirmeier.net}
 */
public final class Matrix
{
    /**
     * Scale a matrix column-wise.
     * <p>
     * Scaling can be either centering the data, dividing by the column's
     * standard deviation
     * or both.
     *
     * @param X the input matrix to scale
     * @param center if true substracts column means for every column
     * @param scale if true divids by column standard deviation for every
     * column
     *
     * @return returns a copy of the scaled matrix
     */
    public static INDArray scale(INDArray X, boolean center, boolean scale)
    {
        INDArray Xn = X.dup();
        INDArray means = Xn.mean(0);
        INDArray vars = null;
        if (scale) vars = Xn.var(0);
        for (int i = 0; i < Xn.columns(); i++)
        {
            INDArray col = Xn.getColumn(i);
            if (center)
                col.assign(col.sub(means.getDouble(i)));
            if (scale)
                col.assign(col.div(sqrt(vars.getDouble(i))));
        }
        return Xn;
    }

    /**
     * Compute a singular value decompositon.
     *
     * @param X matrix that is going to be decomposed
     * @return returns the SVD matrices
     */
    public static SimpleSVD svd(INDArray X)
    {
        return new SimpleMatrix(
          X.rows(), X.columns(), true, X.data().asDouble())
          .svd(true);
    }

}
