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

import net.digital_alexandria.lvm4j.Decomposition;
import org.ejml.simple.SimpleSVD;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

import java.util.ArrayList;
import java.util.List;

import static java.lang.Math.sqrt;
import static net.digital_alexandria.lvm4j.util.Matrix.scale;
import static net.digital_alexandria.lvm4j.util.Matrix.svd;


/**
 * Class that calculates a PCA.
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public final class PCA implements Decomposition
{
    // the input data matrix
    private final INDArray _X;
    // number of rows
    private final int _N;
    // number of columns
    private final int _P;
    // the computed transformation matrix
    private final INDArray _LOADINGS;
    // the standard deviations of each component
    private final List<Double> _SD;
    // the transformed matrix
    private final INDArray _SCORES;

    PCA(final double X[][])
    {
        this(Nd4j.create(X));
    }

    PCA(final INDArray X)
    {
        this._X = X;
        this._N = _X.rows();
        this._P = _X.columns();

        SimpleSVD svd = svd(_X);

        this._LOADINGS = Nd4j.create(
          svd.getV()
             .transpose()
             .getMatrix().getData(),
          new int[]{_P, _P}, 'r');

        // add standard deviations
        this._SD = new ArrayList<>();
        final double sq = sqrt(_X.rows() - 1);
        for (int i = 0; i < _X.columns(); i++)
            _SD.add(svd.getW().get(i, i) / sq);

        this._SCORES = _X.mmul(_LOADINGS);
    }

    /**
     * Computes the rotation matrix of the original dataset using the first
     * k principal components.
     *
     * @param K the number of principal components
     *
     * @return returns the rotation matrix.
     */
    @Override
    public final INDArray run(final int K)
    {
        final int[] cols = new int[K];
        for (int i = 0; i < cols.length; i++)
        {
            cols[i] = i;
        }
        return this._SCORES.getColumns(cols);
    }

    /**
     * Getter for the loadings matrix.
     *
     * @return returns the loadings matrix
     */
    public final INDArray loadings()
    {
        return this._LOADINGS;
    }

    /**
     * Getter for the standard deviations of the singular values.
     *
     * @return returns the standard deviations of the singular values
     */
    public final List<Double> standardDeviations()
    {
        return this._SD;
    }
}
