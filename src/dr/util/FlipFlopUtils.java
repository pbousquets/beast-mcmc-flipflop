/*
 * FlipFlopUtils.java
 *
 * Author: Pablo Bousquets
 */

package dr.util;

import java.io.*;
import java.util.*;
import dr.math.GammaFunction;

public class FlipFlopUtils {

    public static int[] mapDoubleRangeToInt(double[] y, int start, int stop, int num){
        /*
         * Increase num for less precision loss (more noticeable in the power of ten scale) when remapping back to double!
         */
        return digitize(y, linspace(start, stop, num));
    }

    public static double[] linspace(double start, double stop, int num) {
        /*
         * Creates an array of bins between [start,end], including both start and end.
         * Based on numpy.linspace, but very simplified.
         * Original execution example:
         * >> import numpy as np
         * >> np.linspace(0, 1, 5)
         * array([0. , 0.25, 0.5 , 0.75, 1. ])
         */

        if (num < 2) {
            throw new RuntimeException("Error: num is expected to equal at least 2!");
        }

        double div = num - 1;
        double delta = stop - start;
        double step = delta / div;
        double[] res = new double[(int) num];
        for (int x = 0; x < num; x++) {
            res[x] = x * step + start;
        }

        return res;
    }

    public static int[] digitize(double[] x, double[] bin) {
    /*
        Convert normalized double values to integers. Creates num bins across the [start,stop) range.
        For instance, in the 0-1 range, the behaviour for the limits is:
         - 0 will be assigned the bin 1 (the bin 0 will be used for value < start)
         - 1 will be assigned the last bin (last bin will be used for value >= end)
         This is the default behaviour of numpy.digitize

        Example:
        >>> np.digitize([0, 1, .4], np.linspace(0,1, 100))
        array([1, 100, 40])
     */

        for (int i = 0; i < bin.length - 1; i++) {
            if (bin[i] > bin[i + 1])
                throw new RuntimeException("Error: the bin array should be sorted!");
        }

        int[] digitized_x = new int[x.length];

        for (int x_index = 0; x_index < x.length; x_index++) {
            double val = x[x_index];
            int selected_index = bin.length; // Set the max by default, since the loop won't reach it or go any further.
            // The numbers over it will be covered.

            if (val < bin[0]) {
                selected_index = 0; // The lower limit can't be evaluated in the for loop
            } else {

                for (int index = 1; index < bin.length; index++) {
                    if ((bin[index - 1] <= val) & (val < bin[index])) {
                        selected_index = index;
                        break;
                    }
                }
            }

            digitized_x[x_index] = selected_index;
        }
        return digitized_x;
    }

    public static double[] remapDigitized(int[] x, int num){
        /*
         * Takes the digitized values and converts them back to the original double values,
         * assuming a certain precision loss (proportional to num).
         */

        double[] remappedArray = new double[x.length];
        for(int i=0; i<x.length; i++) {
            remappedArray[i] = (double) x[i]/num;
        }
        return remappedArray;
    }

    ////////////////////////////
    // Log-likelihood related //
    ////////////////////////////

    // beta_lpdf constructor for ints
    public static double[] beta_lpdf(double[] y, double[] alpha, double[] beta) {
        //Loop through each of the possible Z=2S+1 peaks and calculate the log-pdf of the nth beta value in y[N] for the zth beta distribution
        int N = y.length; //number of sites
        int Z = alpha.length; //number of states
        double lgamma_alpha;
        double lgamma_beta;
        double lgamma_alphaplusbeta;
        double[] lpk = new double[N * Z];

        for (int z = 0; z < Z; z++) {  //Precompute the log-gamma constants. Iterate through states
            lgamma_alpha = GammaFunction.lnGamma(alpha[z]);
            lgamma_beta = GammaFunction.lnGamma(beta[z]);
            lgamma_alphaplusbeta = GammaFunction.lnGamma(alpha[z] + beta[z]);

            for (int n = 0; n < N; n++) { //Calculate the log pdf of the beta distribution for the nth dataproint drawn from the zth peak. Iterate through sites
                lpk[n * Z + z] = ((alpha[z] - 1) * Math.log(y[n]) + (beta[z] - 1) * StrictMath.log1p(-y[n]) - lgamma_alpha - lgamma_beta + lgamma_alphaplusbeta);
            }
        }
        return lpk;
    }

    // beta_lpdf method for ints
    public static double[] beta_lpdf(int[] y, double[] alpha, double[] beta) {
        //Loop through each of the possible Z=2S+1 peaks and calculate the log-pdf of the nth beta value in y[N] for the zth beta distribution
        int N = y.length; //number of sites
        int Z = alpha.length; //number of states
        double lgamma_alpha;
        double lgamma_beta;
        double lgamma_alphaplusbeta;
        double[] lpk = new double[N * Z];

        for (int z = 0; z < Z; z++) {  //Precompute the log-gamma constants. Iterate through states
            lgamma_alpha = GammaFunction.lnGamma(alpha[z]);
            lgamma_beta = GammaFunction.lnGamma(beta[z]);
            lgamma_alphaplusbeta = GammaFunction.lnGamma(alpha[z] + beta[z]);

            for (int n = 0; n < N; n++) { //Calculate the log pdf of the beta distribution for the nth dataproint drawn from the zth peak. Iterate through sites
                lpk[n * Z + z] = ((alpha[z] - 1) * Math.log(y[n]) + (beta[z] - 1) * StrictMath.log1p(-y[n]) - lgamma_alpha - lgamma_beta + lgamma_alphaplusbeta);
            }
        }
        return lpk;
    }
}