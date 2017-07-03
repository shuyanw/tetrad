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

import net.digital_alexandria.lvm4j.Cluster;
import net.digital_alexandria.lvm4j.Clustering;
import net.digital_alexandria.lvm4j.MixtureComponents;
import net.digital_alexandria.lvm4j.MixtureModel;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;

/**
 * @author Simon Dirmeier {@literal mail@simon-dirmeier.net}
 */
public class GaussianMixture implements Cluster, MixtureModel
{

    private final INDArray _X;
    private final int _N;
    private final int _P;

    GaussianMixture(double[][] X)
    {
        this(Nd4j.create(X));
    }

    GaussianMixture(INDArray X)
    {
        this._X = X;
        this._N = X.rows();
        this._P = X.columns();
    }


    @Override
    public Clustering cluster(int k)
    {
        return null;
    }

    @Override
    public MixtureComponents fit(int k)
    {
        return null;
    }
}
