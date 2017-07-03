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

package net.digital_alexandria.lvm4j.markovmodel;


import net.digital_alexandria.lvm4j.enums.ExitCode;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.Map;
import java.util.TreeMap;

import static net.digital_alexandria.lvm4j.util.Math.sumEquals;
import static net.digital_alexandria.lvm4j.util.System.exit;


/**
 * Class that predict the hidden state sequence given a sequence of
 * observations.
 *
 * @author Simon Dirmeier {@literal s@simon-dirmeier.net}
 */
final class HMMPredictor
{
    private final static Logger _LOGGER =
      LoggerFactory.getLogger(HMMPredictor.class);
    // singleton
    private static HMMPredictor _predictor;

    private HMMPredictor() {}

    static HMMPredictor instance()
    {
        if (_predictor == null)
        {
            _LOGGER.info("Instantiating HMMPredictor");
            _predictor = new HMMPredictor();
        }
        return _predictor;
    }

    /**
     * Predict the most probable latent state sequence using a sequence of
     * observations.  Prediciton is done using the viterbi algorithm.
     *
     * @param observations a mapping from the id of an observation to the real
     * observations sequence
     */
    final Map<String, String> predict(HMM hmm, Map<String, String> observations)
    {
        _LOGGER.info("Predicting latent variables using viterbi!");
        checkMatrices(hmm);
        /* calculate log matrices, because then probabilities can be added
         * instead of multiplied: double has maybe too low precision!
         */
        double[] startProbabilities = hmm.logStartProbabilities();
        double[][] transitionsMatrix = hmm.logTransitionMatrix();
        double[][] emissionsMatrix = hmm.logEmissionMatrix();
        Map<String, String> statesMap = new TreeMap<>();
        for (Map.Entry<String, String> entry : observations.entrySet())
        {
            char[] obs = entry.getValue().toUpperCase().toCharArray();
            String stateSequence = viterbi(hmm, obs,
                                           startProbabilities,
                                           transitionsMatrix,
                                           emissionsMatrix);
            statesMap.put(entry.getKey(), stateSequence);
        }
        return statesMap;
    }

    private void checkMatrices(HMM hmm)
    {
        final double delta = 0.01;
        final double probSum = 1.0;
        final double altProbSum = 0.0;
        double[] startProbabilities = hmm.startProbabilities();
        if (!sumEquals(startProbabilities, delta, probSum))
            exit("Sum of starting probabilities does not equal 1.00!",
                 ExitCode.EXIT_ERROR);
        double[][] transitionsMatrix = hmm.transitionMatrix();
        for (double row[] : transitionsMatrix)
        {
            if (!(sumEquals(row, delta, probSum) ||
                  sumEquals(row, delta, altProbSum)))
                exit("Sum of transition probabilities does not equal 1.00!",
                     ExitCode.EXIT_ERROR);
        }
        double[][] emissionsMatrix = hmm.emissionMatrix();
        for (double row[] : emissionsMatrix)
        {
            if (!(sumEquals(row, delta, probSum) ||
                  sumEquals(row, delta, altProbSum)))
                exit("Sum of emission probabilities does not equal 1.00!",
                     ExitCode.EXIT_ERROR);
        }
    }

    /**
     * Calculate the most probable state sequence.
     *
     * @param obs array of observations
     * @param startProbabilities array of starting probabilities for the states
     * @param transitionsMatrix matrix of state transition probabilities
     * @param emissionMatrix matrix of emission probabilities
     *
     * @return returns the sequence of predicted states
     */
    private String viterbi(HMM hmm,
                           char[] obs,
                           double[] startProbabilities,
                           double[][] transitionsMatrix,
                           double[][] emissionMatrix)

    {
        // code array of characters (observations) to integers for faster
        // access
        int encodedObservations[] = charToInt(hmm, obs);
        // matrix of paths of probabilities
        double probabilityPath[][] = new double[hmm.states().size()][obs
          .length];
        // matrix of paths of states
        int statePath[][] = new int[hmm.states().size()][obs.length];
        // set the first element of state/probability paths
        initStarts(hmm, probabilityPath, statePath, startProbabilities,
                   emissionMatrix, encodedObservations);
        // fill the two matrices with state and probability paths
        fillPathMatrices(hmm, statePath, probabilityPath, encodedObservations,
                         transitionsMatrix, emissionMatrix);
        // backtrace the matrices to predict the most probable state sequence
        return backtrace(hmm, probabilityPath, statePath, encodedObservations);
    }

    private void fillPathMatrices(HMM hmm,
                                  int[][] statePath,
                                  double[][] probabilityPath,
                                  int[] encodedObservations,
                                  double[][] transitionsMatrix,
                                  double[][] emissionMatrix)
    {
        for (int i = 1; i < encodedObservations.length; i++)
        {
            for (int j = 0; j < hmm.states().size(); j++)
                setPaths(statePath,
                         probabilityPath,
                         i, j,
                         transitionsMatrix,
                         emissionMatrix,
                         encodedObservations[i]);
        }
    }

    private void initStarts(HMM hmm,
                            double[][] probs,
                            int[][] states,
                            double[] startProbabilities,
                            double[][] emissionMatrix,
                            int[] encodedObservations)
    {
        for (int i = 0; i < hmm.states().size(); i++)
        {
            probs[i][0] = startProbabilities[i] +
                          emissionMatrix[i][encodedObservations[0]];
            states[i][0] = 0;
        }
    }

    /**
     * Backtrace the probability and state path matrices to create the most
     * probable latent state sequence.
     *
     * @param probabilityPath matrix of paths of probabilities
     * @param statePath matrix of paths of states
     *
     * @return return a string of the predicted latent states
     */
    private String backtrace(HMM hmm,
                             final double[][] probabilityPath,
                             final int[][] statePath,
                             int[] encodedObservations)
    {
        // sequence of most probable state indeces
        int statesIdx[] = new int[encodedObservations.length];
        // sequence of according labels
        char statesLabel[] = new char[encodedObservations.length];
        // set idx of last state
        statesIdx[encodedObservations.length - 1] =
          getMostProbableEndingStateIdx(probabilityPath);
        // set label of last state
        statesLabel[encodedObservations.length - 1] = hmm.states().get
          (statesIdx[encodedObservations.length - 1]).label();
        // backtrace state/probability paths to get most probable state
        // sequence
        for (int i = encodedObservations.length - 1; i >= 1; i--)
        {
            statesIdx[i - 1] = statePath[statesIdx[i]][i];
            statesLabel[i - 1] = hmm.states().get(statesIdx[i - 1]).label();
        }
        return String.valueOf(statesLabel);
    }

    // nice function name
    private int getMostProbableEndingStateIdx(double[][] probs)
    {
        double prob = Double.MIN_VALUE;
        int idx = 0;
        for (int i = 0; i < probs.length; i++)
        {
            if (probs[i][probs[i].length - 1] > prob)
            {
                idx = i;
                prob = probs[i][probs[i].length - 1];
            }
        }
        return idx;
    }

    private void setPaths(int[][] statePath,
                          double[][] probabilityPath,
                          final int i, int j,
                          final double[][] transitionsMatrix,
                          final double[][] emissionMatrix, int o)
    {
        double trans = Double.MIN_VALUE;
        int idx = 0;
        // calculate the most probable state following after state s[i-1]
        for (int k = 0; k < transitionsMatrix.length; k++)
        {
            double curr = probabilityPath[k][i - 1] +
                          transitionsMatrix[k][j] +
                          emissionMatrix[j][o];
            if (curr > trans)
            {
                trans = curr;
                idx = k;
            }
        }
        probabilityPath[j][i] = trans;
        statePath[j][i] = idx;
    }

    /**
     * Convert a char array to an array of int indexes.
     *
     * @param hmm the markovmodel
     * @param obs the array of observations to be encoded
     *
     * @return the array of encoded integers
     */
    private int[] charToInt(HMM hmm, char[] obs)
    {
        int idxs[] = new int[obs.length];
        for (int i = 0; i < idxs.length; i++)
        {
            for (int j = 0; j < hmm.observations().size(); j++)
            {
                if (obs[i] == hmm.observations().get(j).label())
                    idxs[i] = hmm.observations().get(j).idx();
            }
        }
        return idxs;
    }
}
