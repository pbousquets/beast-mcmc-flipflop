package test.dr.evomodel.flipflop;

import dr.evolution.datatype.AFsequence;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 * Test the flipflop model with relative and absolute rates
 *
 * @author Diego Mallo
 */
public class TestFlipFlopModelRates extends TestCase {

    class Test {
        public Test(double age, int nCells, double gamma, double lambda, double mu, boolean useNormalization, boolean useFreqModel, FrequencyModel freqModel){
            this.age = age;
            this.nCells = nCells;
            this.gamma = gamma;
            this.lambda = lambda;
            this.mu = mu;
            this.useNormalization = useNormalization;
            this.nStates = (int) (0.5 * (nCells + 1) * (nCells + 2));

            AFsequence afseq = new AFsequence(this.nStates);

            if (freqModel == null){
                double[] freqs = new double[nStates];
                for (int i = 0; i < freqs.length; i++) {
                    freqs[i] = 1.0 / nStates;
                }
                freqModel = new FrequencyModel(afseq, freqs);
            }

            this.useFreqModel = useFreqModel;

            Parameter stemCellParameter = new Parameter.Default(1, nCells);
            Parameter gammaParam = new Parameter.Default(1, gamma);
            Parameter lambdaParam = new Parameter.Default(1, lambda);
            Parameter muParam = new Parameter.Default(1, mu);

            absoluteModel = new FlipFlopModel("absolute", afseq, stemCellParameter, gammaParam, lambdaParam, muParam, useNormalization, useFreqModel, freqModel);

            Parameter relativeGammaParam = new Parameter.Default(1, gamma/mu);
            Parameter relativeLambdaParam = new Parameter.Default(1,lambda/mu);
            Parameter relativeMuParam = new Parameter.Default(1, mu/mu);

            relativeModel = new FlipFlopModel("relative", afseq, stemCellParameter, relativeGammaParam, relativeLambdaParam,relativeMuParam, useNormalization, useFreqModel, freqModel);
        }

        public void getAbsoluteProbabilityMatrix(double [] matrix) {
            absoluteModel.getTransitionProbabilities(this.age, matrix);
        }

        public void getRelativeProbabilityMatrix(double [] matrix) {
            relativeModel.getTransitionProbabilities(this.age*this.mu, matrix);
        }

        public double [] getRelativeTransitionMatrix() {
            return flattenRateMatrix(relativeModel.getAmat());
        }

        public double [] getAbsoluteTransitionMatrix() {
            return flattenRateMatrix(absoluteModel.getAmat());
        }

        private double [] flattenRateMatrix(double [][]matrix){
            double [] flatmat = new double [nStates*nStates];
            for (int k = 0; k < nStates; ++k) {
                for (int i = 0; i < nStates; ++i) {
                    flatmat[nStates*k+i] = matrix[k][i];
                }
            }
            return flatmat;
        }

        public int getNstates(){return nStates;};

        public FlipFlopModel getRelativeModel(){
            return(relativeModel);
        };

        public FlipFlopModel getAbsoluteModel(){
            return(absoluteModel);
        };

        protected final double age;
        protected final int nCells;
        protected final int nStates;
        protected final double gamma;
        protected final double lambda;
        protected final double mu;
        protected final boolean useNormalization;
        protected final boolean useFreqModel;

        protected final FlipFlopModel absoluteModel;
        protected final FlipFlopModel relativeModel;
    };

    // TEST 1

    Test test1 = new Test(2, 3, 0.05, 1.3, 0.05, false, false, null);
    Test test2 = new Test(0.05, 3, 0.05, 1.3, 0.05, false, false, null);
    Test test3 = new Test(10, 3, 0.05, 1.3, 0.05, false, false, null);



    Test[] all = {test1,test2,test3}; //Add more tests here

    public void testTransitionMatrix(){
        for (Test test : all) {
            double [] relativeMatrix = new double [test.getNstates() * test.getNstates()];
            test.getRelativeProbabilityMatrix(relativeMatrix);

            double [] absoluteMatrix = new double [test.getNstates() * test.getNstates()];
            test.getAbsoluteProbabilityMatrix(absoluteMatrix);

            for (int j = 0; j < test.getNstates()*test.getNstates(); ++j) {
                assertEquals(relativeMatrix[j], absoluteMatrix[j], 1e-7);
            }
        }
    }
}


