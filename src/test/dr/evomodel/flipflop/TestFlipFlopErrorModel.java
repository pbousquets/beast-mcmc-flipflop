package test.dr.evomodel.flipflop;
import java.util.Arrays;
import dr.util.FlipFlopUtils;

import junit.framework.TestCase;

/**
 * Test the flipflop error model
 *
 * @author Pablo Bousquets
 */
public class TestFlipFlopErrorModel extends TestCase {

    interface Instance {
        double getStemCellParameter();
        double getDeltaParameter();
        double getEtaParameter();
        double getKappaParameter();

        double[] getSequence();

        double[] getExpectedResult();
    }

    Instance test0 = new Instance() {
        public double getStemCellParameter(){
            return 2;
        };
        public double getDeltaParameter(){
            return 0.05;
        };
        public double getEtaParameter(){
            return 0.95;
        };
        public double getKappaParameter(){
            return 0.9;
        };

        public double[] getSequence(){
            return new double[]{
                    0.1, 0.4, 0.14, 0.12, 0.99, .05, 0.01
            };
        };

        public double[] getExpectedIdealAlpha() {
            return new double[]{
                    0.045 , 0.2475, 0.45  , 0.6525, 0.855
            };
        }
        public double[] getExpectedIdealBeta() {
            return new double[]{
                    0.855 , 0.6525, 0.45  , 0.2475, 0.045
            };
        }
        public double[] getExpectedResult() {
            return new double[]{
                    -0.8984707443537294, 0.21490279494496956, 0.03657211353834405, -0.6749731588762, -2.678222651996067,
                    -2.163589418547541, -0.6873845867301609, -0.5028839756181054, -0.8515979555139683, -2.492016156115154,
                    -1.2132096860858579, -0.02249438812081994, -0.12348331086110648, -0.7576868246093116, -2.6835945590628394,
                    -1.0693292669934076, 0.08551514086610923, -0.051344671979698886, -0.7214190758334251, -2.6831977003924745,
                    -2.4353690352250608, 0.05345425063647635, 1.2505733157926897, 1.9144777899409842, 1.2866780433839557,
                    -0.24435493400317204, 0.7177076889249079, 0.3880660911476625, -0.47479009763750163, -2.629350507127988,
                    1.2866780433839566, 1.9144777899409862, 1.2505733157926906, 0.05345425063647652, -2.4353690352250603
            };
        }

    };

    Instance[] all = {test0}; //Add more instances of tests here

    public void testErrorModel() {
        for (Instance test : all) {
            // Create the parameters
            double delta = test.getDeltaParameter();
            double eta = test.getEtaParameter();
            double kappa = test.getKappaParameter();
            double cells = test.getStemCellParameter();

            // Get the pattern
            double[] seq = test.getSequence();

            // Compute the partials!
            int total_states = 2 * (int) cells + 1;
            double[] ideal_beta = new double [total_states];
            double[] beta = new double [total_states];
            double[] transformed_alpha = new double [total_states];
            double[] transformed_beta = new double [total_states];

            for (int state = 0; state < total_states; state ++){
                ideal_beta[state] = ((double) state) / (total_states - 1); // We'll have 0/N methylated alleles, 1/N allele, 2/N alleles, ... N/N methylated alleles (considering here N = total alleles)
                beta[state] = (eta - delta) * ideal_beta[state] + delta; // Rescale the expected peaks to linear scale

                //Convert mean/dispersion parameterization of a beta distribution
                transformed_alpha[state] = beta[state] * kappa;
                transformed_beta[state] =  (1 - beta[state]) * kappa;
            }

            //Compute log(p(y|z))
            double[] partials = FlipFlopUtils.beta_lpdf(seq, transformed_alpha, transformed_beta); //Format: partials per site per state: [site1state1,site1state2,site1state3,...,site2state1,site2state2,site2state3,...]
            System.out.println(Arrays.toString(partials));

            double[] expected = test.getExpectedResult();
            System.out.println(Arrays.toString(expected));
            assertEquals(expected.length, partials.length, 1e-15);
            
            for (int kk = 0; kk < expected.length; kk+=1) {
                assertEquals(expected[kk], partials[kk], 1e-5);
            }

        }
    };

}
