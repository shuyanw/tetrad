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
 * along with lvm4j. If not, see <http://www.gnu.org/licenses/>.
 */


package net.digital_alexandria.lvm4j.decomposition;

import net.digital_alexandria.lvm4j.Decomposition;
import org.ejml.simple.SimpleMatrix;
import org.ejml.simple.SimpleSVD;
import org.nd4j.linalg.api.ndarray.INDArray;
import org.nd4j.linalg.factory.Nd4j;
import org.nd4j.linalg.inverse.InvertMatrix;
import org.nd4j.linalg.ops.transforms.Transforms;

import static java.lang.Double.max;
import static java.lang.Math.sqrt;
import static net.digital_alexandria.lvm4j.util.Math.log;
import static net.digital_alexandria.lvm4j.util.Math.sum;
import static net.digital_alexandria.lvm4j.util.Matrix.scale;
import static org.nd4j.linalg.factory.Nd4j.diag;

/**
 * Class that calculates a factor analysis.
 *
 * The algorithm was largely taken from:
 *  Barber - Bayesian Reasoning and Machine Learning,
 *  Chapter 21. 1, Algorithm 21. 1
 *
 * @author Simon Dirmeier {@literal simon.dirmeier@gmx.de}
 */
public final class FactorAnalysis implements Decomposition
{
    // stop criterion: threshold if likelihood does not change any more
    private static final double _THRESHOLD = 0.0001;
    // stop criterion: maximal number of iterations
    private static final int _MAXIT = 10000;
    // a small pseudocount for numerical stability
    private static final double PSEUDO_COUNT = 1e-12;

    // the input data matrix
    private final INDArray _X;
    // number of rows
    private final int _N;
    // number of columns
    private final int _P;

    // factor loading matrix
    private INDArray _f;
    // vcov matrix of Gaussian errors
    private INDArray _psi;

    FactorAnalysis(final double X[][])
    {
        this(Nd4j.create(X));
    }

    FactorAnalysis(final INDArray X)
    {
        this._X = scale(X, true, false);
        this._N = _X.rows();
        this._P = _X.columns();
    }

    /**
     * Runs the factor analysis and computes the score matrix with K factors.
     *
     * @param K the number of factors that are computed for the score matrix.
     *
     * @return returns the transformed matrix
     */
    @Override
    public final INDArray run(final int K)
    {
        fit(K);
        return decomp(K);
    }

    private void fit(final int K)
    {
        final INDArray vars = _X.var(0);
        INDArray F;
        INDArray psis = Nd4j.eye(_P);

        double loglik = Double.MIN_VALUE;
        double oldLoglik;

        int niter = 0;
        do
        {
            oldLoglik = loglik;
            final INDArray sqrt_psis = sqrtPsis(psis);

            INDArray X = _X.dup();
            SimpleSVD svd;
            double unexp = .0;
            INDArray s, V;

            // TODO: CHANGE TO OWN METHOD USING
            {
                svd = svd
                  (X, sqrt_psis.data().asDouble());
                s = getSingularValues(svd.getW(), K);
                V = getRightSingularVectors(svd.getV(), K);
                unexp = unexplainedVariance(svd.getW(), K);
            }

            // update the likelihood
            loglik = proploglik(s, unexp, psis);
            // update the factor matrix and variances AFTERWARDS
            F = factorUpdate(s, V, sqrt_psis);
            psis = vcovUpdate(vars, F);
        }
        while (niter++ < _MAXIT && Math.abs(loglik - oldLoglik) > _THRESHOLD);

        // set the member variables to the computed values
        this._f = F;
        this._psi = psis;
    }

    private INDArray sqrtPsis(INDArray psis)
    {
        return Transforms.sqrt(diag(psis)).add(PSEUDO_COUNT);
    }

    private SimpleSVD svd(final INDArray X, final double[] sdevs)
    {
        final double nsqrt = sqrt(_N);
        for (int i = 0; i < sdevs.length; i++)
        {
            X.getColumn(i).assign(X.getColumn(i).div(sdevs[i] * nsqrt));
        }
        return new SimpleMatrix(_N, _P, true, X.data().asDouble()).svd(true);
    }

    private INDArray getSingularValues(final SimpleMatrix S,
                                       final int k)
    {
        final double[] diag = S.extractDiag().getMatrix().data;
        INDArray s = Nd4j.create(k);
        for (int i = 0; i < k; i++) s.getColumn(i).assign(Math.pow(diag[i], 2));
        return s;
    }

    private INDArray getRightSingularVectors(final SimpleMatrix V,
                                             final int K)
    {
        return Nd4j.create
          (V.transpose().extractMatrix(0, K, 0, V.numCols()).getMatrix().data,
           new int[]{K, V.numCols()}, 'r');
    }

    private double unexplainedVariance(final SimpleMatrix S, final int k)
    {
        double var = 0.0;
        for (int i = k; i < S.numCols(); i++)
        {
            final double sv = S.get(i, i);
            var += sv * sv;
        }
        return var;
    }

    private INDArray factorUpdate(final INDArray singVals,
                                  final INDArray V,
                                  final INDArray sqrt_psis)
    {
        final int sNcol = singVals.columns();
        INDArray m = Nd4j.create(sNcol);
        for (int i = 0; i < sNcol; i++)
        {
            m.getColumn(i)
             .assign(sqrt(Math.max(singVals.getDouble(i) - 1, 0.)));
        }
        m = Nd4j.diag(m);
        INDArray W = V.transpose().mmul(m).transpose();
        for (int i = 0; i < W.columns(); i++)
        {
            W.getColumn(i).assign(W.getColumn(i).mul(sqrt_psis.getDouble(i)));
        }
        return W;
    }

    private INDArray vcovUpdate(final INDArray var, final INDArray w)
    {
        INDArray fft = Transforms.pow(w, 2);
        INDArray psi = fft.sum(0);
        for (int i = 0; i < psi.columns(); i++)
        {
            psi
              .getColumn(i)
              .assign(max(var.getDouble(i) - psi.getDouble(i), PSEUDO_COUNT));
        }
        return Nd4j.diag(psi);
    }

    private double proploglik(final INDArray s,
                              final double unexp,
                              final INDArray psis)
    {
        return sum(log(s.data().asDouble())) +
               unexp +
               sum(log(diag(psis).data().asDouble()));
    }

    private INDArray decomp(final int K)
    {
        INDArray m = Nd4j.eye(K);
        INDArray psiw = wpsi();
        INDArray b = psiw.mmul(this._f.transpose());
        INDArray coz = coz(K, m, b);
        return _X.mmul(psiw.transpose()).mmul(coz);
    }

    private INDArray wpsi()
    {
        INDArray psiw = this._f.dup();
        INDArray ps = diag(this._psi);
        for (int i = 0; i < psiw.rows(); i++)
            psiw.getRow(i).assign(psiw.getRow(i).div(ps));
        return psiw;
    }

    private INDArray coz(final int K, final INDArray m, final INDArray b)
    {
        return (K == 1)
          ? Nd4j.ones(1).div(m.add(b))
          : InvertMatrix.invert(m.add(b), false);
    }
}
