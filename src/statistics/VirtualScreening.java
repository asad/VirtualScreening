/*
 * Copyright (c) 2015, Asad
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright notice, this
 *   list of conditions and the following disclaimer.
 * * Redistributions in binary form must reproduce the above copyright notice,
 *   this list of conditions and the following disclaimer in the documentation
 *   and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
package statistics;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;

/**
 * Java code for enrichvs package in R
 *
 * @author Syed Asad Rahman <asad@ebi.ac.uk>
 */
public class VirtualScreening {

    private static final Logger LOG = Logger.getLogger(VirtualScreening.class.getName());
    private final double alpha = 0.20d;
    private final double top = 0.05d;
    private final double[] scores;
    private final boolean[] lables;
    private final boolean DEBUG = false;

    /**
     *
     * @param scores {-99.9,90.0,-98.0}
     * @param lables {true,false,true....}
     */
    public VirtualScreening(double[] scores, boolean[] lables) {
        this.scores = scores;
        this.lables = lables;
    }

    /*
     * Calculate a RIE for the highly-ranked compounds
     * Ref.: Truchon et al. Evaluating Virtual Screening Methods: 
     * Good and Bad Metrics for the "Early Recognition" Problem.
     * J. Chem. Inf. Model. (2007) 47, 488-508.
     */

    /**
     * Function to culculate the Robust Initial Enhancement (RIE)
     * @param alpha coefficient alpha
     * @param decreasing TRUE if the compounds are ranked by decreasing score
     * @return RIE, in the range from 0 to +Inf.
     */
    
    public double rie(double alpha, boolean decreasing) {
        if (scores.length != lables.length) {
            System.err.println("The number of scores must be equal to the number of labels.");
            return -1.0f;
        }

        double N = lables.length;
        if (DEBUG) {
            System.err.println("N: " + N);
        }
        double n = getPositiveHitsCount(lables);
        if (DEBUG) {
            System.err.println("n: " + n);
        }
        int[] order;
        if (decreasing) {
            order = getDecreasingOdering(scores);
        } else {
            order = getIncreasingOdering(scores);
        }
        if (DEBUG) {
            System.err.println("Rank order: " + order.length);
        }
        int[] m_rank = getPositiveHitsRankedIndex(order, lables);
        if (DEBUG) {
            System.err.println("m_rank: " + m_rank.length);
        }

        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += Math.exp(-alpha * m_rank[i] / N);
        }
        if (DEBUG) {
            System.err.println("Sum: " + s);
        }
        double random_sum = (n / N) * (1 - Math.exp(-alpha)) / (Math.exp(alpha / N) - 1);
        return (s / random_sum);
    }


    /*
     * Calculate a BEDROC for the highly-ranked compounds
     * Ref.: Truchon et al. Evaluating Virtual Screening Methods: 
     * Good and Bad Metrics for the "Early Recognition" Problem.
     * J. Chem. Inf. Model. (2007) 47, 488-508.
     */

    /**
     * Boltzmann-Enhanced Discrimination of ROC (BEDROC)
     * @param alpha coefficient alpha
     * @param decreasing TRUE if the compounds are ranked by decreasing score
     * @return BEDROC, in the range from 0 to 1.
     */
    
    public double bedroc(double alpha, boolean decreasing) {
        if (scores.length != lables.length) {
            System.err.println("The number of scores must be equal to the number of labels.");
            return -1.0f;
        }
        double N = lables.length;
        if (DEBUG) {
            System.err.println("N: " + N);
        }
        double n = getPositiveHitsCount(lables);
        if (DEBUG) {
            System.err.println("n: " + n);
        }
        int[] order;
        if (decreasing) {
            order = getDecreasingOdering(scores);
        } else {
            order = getIncreasingOdering(scores);
        }
        if (DEBUG) {
            System.err.println("Rank order: " + order.length);
        }
        int[] m_rank = getPositiveHitsRankedIndex(order, lables);
        if (DEBUG) {
            System.err.println("m_rank: " + m_rank.length);
        }

        double s = 0.0;
        for (int i = 0; i < n; i++) {
            s += Math.exp(-alpha * m_rank[i] / N);
        }
        if (DEBUG) {
            System.err.println("Sum: " + s);
        }
        double ra = n / N;
        double ri = (N - n) / N;

        double random_sum = ra * Math.exp(-alpha / N) * (1.0 - Math.exp(-alpha)) / (1.0 - Math.exp(-alpha / N));
        double fac = ra * Math.sinh(alpha / 2.0) / (Math.cosh(alpha / 2.0) - Math.cosh(alpha / 2.0 - alpha * ra));
        double cte = 1.0 / (1 - Math.exp(alpha * ri));
        return (s / random_sum * fac + cte);
    }

    /*
     * Calculate a enrichment factor for the highly-ranked compounds
     * Ref.: Truchon et al. Evaluating Virtual Screening Methods: 
     * Good and Bad Metrics for the "Early Recognition" Problem.
     * J. Chem. Inf. Model. (2007) 47, 488-508.
     */

    /**
     * Function to calculate the enrichment factor (EF)
     * @param top threshold ratio of the false positives (when ROC analysis is performed on a top list)
     * @param decreasing TRUE if the compounds are ranked by decreasing score
     * @return EF, in the range from 0 to +Inf.
     */
    
    public double enrichment_factor(double top, boolean decreasing) {
        if (scores.length != lables.length) {
            System.err.println("The number of scores must be equal to the number of labels.");
            return -1.0f;
        }
        double N = lables.length;
        if (DEBUG) {
            System.err.println("N: " + N);
        }
        double n = getPositiveHitsSum(lables);
        if (DEBUG) {
            System.err.println("n: " + n);
        }

        double x_prev = Double.NEGATIVE_INFINITY;
        double fp = 0.0d, tp = 0.0d, fp_prev = 0.0d, tp_prev = 0.0d;

        int[] order;
        if (decreasing) {
            order = getDecreasingOdering(scores);
        } else {
            order = getIncreasingOdering(scores);
        }

        if (DEBUG) {
            System.err.println("Rank order: " + order.length);
        }

        for (int i : order) {
            int j = order[i];
            if (scores[j] != x_prev) {
                if (fp + tp >= N * top) {
                    double n_right = ((fp - fp_prev) + (tp - tp_prev));
                    double rat = (N * top - (fp_prev + tp_prev)) / n_right;
                    double tp_r = tp_prev + rat * (tp - tp_prev);
                    return ((tp_r / (N * top)) / (n / N));
                }
                x_prev = scores[j];
                fp_prev = fp;
                tp_prev = tp;
            }
            if (lables[j]) {
                tp = tp + 1;
            } else {
                fp = fp + 1;
            }
        }
        //  n_right <- (fp - fp_prev) + (tp - tp_prev)
        return (1);
    }

    /*
     * Calculate the AUC of ROC curve for ranked compounds
     * Ref.: Tom Fawcett, An introduction to ROC analysis. 
     * Pattern Recognition Letters 27, 861-874 (2006)
     */

    /**
     * Function to calculate the Area Under the ROC Curve (AUC)
     * @param top threshold ratio of the false positives (when ROC analysis is performed on a top list)
     * @param decreasing TRUE if the compounds are ranked by decreasing score
     * @return AUC, in the range from 0 to 1.
     */
    
    public double auc(double top, boolean decreasing) {
        if (scores.length != lables.length) {
            System.err.println("The number of scores must be equal to the number of labels.");
            return -1.0f;
        }
        double N = lables.length;
        if (DEBUG) {
            System.err.println("N: " + N);
        }
        double n = getPositiveHitsSum(lables);
        if (DEBUG) {
            System.err.println("n: " + n);
        }

        double x_prev = Double.NEGATIVE_INFINITY;
        double fp = 0.0d, tp = 0.0d, fp_prev = 0.0d, tp_prev = 0.0d;
        double area = 0.0d;
        int[] order;
        if (decreasing) {
            order = getDecreasingOdering(scores);
        } else {
            order = getIncreasingOdering(scores);
        }

        if (DEBUG) {
            System.err.println("Rank order: " + order.length);
        }

        for (int i : order) {
            int j = order[i];
            if (scores[j] != x_prev) {
                if (fp >= (N - n) * top) {
                    double rat = ((N - n) * top - fp_prev) / (fp - fp_prev);
                    area = area + rat * (fp - fp_prev) * (tp + tp_prev) / 2.0d;
                    return (area / (n * (N - n) * top));
                }
                area = area + (fp - fp_prev) * (tp + tp_prev) / 2.0d;
                x_prev = scores[j];
                fp_prev = fp;
                tp_prev = tp;
            }
            if (lables[j]) {
                tp = tp + 1;
            } else {
                fp = fp + 1;
            }
        }
        area = area + (fp - fp_prev) * (tp + tp_prev) / 2.0d;
        return (area / (n * (N - n)));
    }

    /*
     * Calculate the Area Under the Accumulation Curve (AUAC)
     * Ref.: Tom Fawcett, An introduction to ROC analysis. 
     * Pattern Recognition Letters 27, 861-874 (2006)
     */

    /**
     * Function to calculate the Area Under the Accumulation Curve (AUAC)
     * @param top threshold ratio of the false positives (when ROC analysis is performed on a top list)
     * @param decreasing TRUE if the compounds are ranked by decreasing score
     * @return AUAC, in the range from 0 to 1.
     */
    
    public double auac(double top, boolean decreasing) {
        if (scores.length != lables.length) {
            System.err.println("The number of scores must be equal to the number of labels.");
            return -1.0f;
        }
        double N = lables.length;
        if (DEBUG) {
            System.err.println("N: " + N);
        }
        double n = getPositiveHitsSum(lables);
        if (DEBUG) {
            System.err.println("n: " + n);
        }

        double x_prev = Double.NEGATIVE_INFINITY;
        double fp = 0.0d, tp = 0.0d, fp_prev = 0.0d, tp_prev = 0.0d;
        double area = 0.0d;
        int[] order;
        if (decreasing) {
            order = getDecreasingOdering(scores);
        } else {
            order = getIncreasingOdering(scores);
        }

        if (DEBUG) {
            System.err.println("Rank order: " + order.length);
        }

        for (int i : order) {
            int j = order[i];
            if (scores[j] != x_prev) {
                if (fp + tp >= N * top) {
                    double n_right = (fp - fp_prev) + (tp - tp_prev);
                    double rat = (N * top - (fp_prev + tp_prev)) / n_right;
                    area = area + rat * n_right * (tp + tp_prev) / 2.0d;
                    return (area / (n * N * top));
                }
                double n_right = (fp - fp_prev) + (tp - tp_prev);
                area = area + n_right * (tp + tp_prev) / 2.0d;
                x_prev = scores[j];
                fp_prev = fp;
                tp_prev = tp;
            }
            if (lables[j]) {
                tp = tp + 1;
            } else {
                fp = fp + 1;
            }
        }
        double n_right = (fp - fp_prev) + (tp - tp_prev);
        area = area + n_right * (tp + tp_prev) / 2.0d;
        return (area / (n * N));
    }

    private double getPositiveHitsSum(boolean[] lables) {
        double count = 0;
        for (boolean b : lables) {
            if (b) {
                count++;
            }
        }
        return count;
    }

    private double getPositiveHitsCount(boolean[] lables) {
        double count = 0;
        for (boolean b : lables) {
            if (b) {
                count++;
            }
        }
        return count;
    }

    private int[] getPositiveHitsRankedIndex(int[] order, boolean[] lables) {
        List<Integer> posIndexList = new LinkedList<>();
        for (int i = 0; i < lables.length; i++) {
            if (lables[i]) {
                posIndexList.add(order[i]);
            }
        }
        int[] posIndex = new int[posIndexList.size()];
        for (int i = 0; i < posIndexList.size(); i++) {
            posIndex[i] = posIndexList.get(i);
        }
        return posIndex;
    }

    /*
     Returns the index of scores in desceding order
     */
    int[] getDecreasingOdering(double[] scores) {
        NaturalRanking naturalRanking = new NaturalRanking(NaNStrategy.FIXED,
                TiesStrategy.AVERAGE);
//        System.out.println("scores " + Arrays.toString(scores));
        final double[] ranks = naturalRanking.rank(scores);
//        System.out.println("ranks " + Arrays.toString(ranks));
        boolean[] visitorFlag = new boolean[scores.length];
        for (int i = 0; i < scores.length; i++) {
            visitorFlag[i] = false;
        }
        double[] sortedArray = ranks.clone();
        Arrays.sort(sortedArray);

        /*
         Reverse the array to make it in descending order
         */
        for (int i = 0; i < sortedArray.length / 2; i++) {
            double temp = sortedArray[i];
            sortedArray[i] = sortedArray[sortedArray.length - i - 1];
            sortedArray[sortedArray.length - i - 1] = temp;
        }

        //System.out.println("sortedArray " + Arrays.toString(sortedArray));
        final int[] rankIndex = new int[scores.length];
        for (int i = sortedArray.length - 1; i >= 0; i--) {
            double value = sortedArray[i];
            int index = getIndex(visitorFlag, value, ranks);
            rankIndex[i] = index;
            visitorFlag[index] = true;
        }
        return rankIndex;
    }

    /*
     Returns the index of scores in desceding order
     */
    int[] getIncreasingOdering(double[] scores) {
        NaturalRanking naturalRanking = new NaturalRanking(NaNStrategy.FIXED,
                TiesStrategy.AVERAGE);
//        System.out.println("scores " + Arrays.toString(scores));
        final double[] ranks = naturalRanking.rank(scores);
//        System.out.println("ranks " + Arrays.toString(ranks));
        boolean[] visitorFlag = new boolean[scores.length];
        for (int i = 0; i < scores.length; i++) {
            visitorFlag[i] = false;
        }
        double[] sortedArray = ranks.clone();
        Arrays.sort(sortedArray);

        //System.out.println("sortedArray " + Arrays.toString(sortedArray));
        final int[] rankIndex = new int[scores.length];
        for (int i = sortedArray.length - 1; i >= 0; i--) {
            double value = sortedArray[i];
            int index = getIndex(visitorFlag, value, ranks);
            rankIndex[i] = index;
            visitorFlag[index] = true;
        }
        return rankIndex;
    }

    int getIndex(boolean[] b, double value, double[] ranks) {
        for (int i = 0; i < ranks.length; i++) {
            if (ranks[i] == value && !b[i]) {
                return i;
            }
        }
        return -1;
    }
}

