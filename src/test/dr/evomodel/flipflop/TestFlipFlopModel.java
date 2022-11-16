package test.dr.evomodel.flipflop;

import dr.evolution.datatype.AFsequence;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.util.Arrays;

/**
 * Test the flipflop error model
 *
 * @author Pablo Bousquets
 */
public class TestFlipFlopModel extends TestCase {

    interface Instance {

        int getNstates();
        double getStemCellParameter();
        double getGammaParameter();
        double getLambdaParameter();
        double getMuParameter();
        boolean getUseFrequencyModel();
        double[] getExpectedResult();

    }

    Instance test0 = new Instance() {
        //MODEL
        protected final int nCells=2;
        protected final int nStates = (int) (0.5 * (nCells + 1) * (nCells + 2));
        protected final double gamma=0.05;
        protected final double lambda=0.95;
        protected final double mu=0.05;
        protected final boolean useFreqModel = true;
        //EXPECTED
        protected final double [] expectedResults = {
                -0.07333333333333335, 0.07333333333333335, 0.0, 0.0, 0.0, 0.0,
                0.16666666666666666, -0.18999999999999997, 0.017499999999999998, 0.005833333333333333, 0.0, 0.0,
                0.0, 0.036666666666666674, -0.053333333333333344, 0.0, 0.016666666666666666, 0.0,
                0.15833333333333333, 0.036666666666666674, 0.0, -0.37, 0.016666666666666666, 0.15833333333333333,
                0.0, 0.0, 0.017499999999999998, 0.005833333333333333, -0.19, 0.16666666666666666,
                0.0, 0.0, 0.0, 0.0, 0.03333333333333333, -0.03333333333333333
        };

        //GETTERS
        public int getNstates(){return this.nStates;};
        public double getStemCellParameter(){
            return this.nCells;
        };
        public double getGammaParameter(){
            return this.gamma;
        };
        public double getLambdaParameter(){
            return this.lambda;
        };
        public double getMuParameter(){
            return this.mu;
        };
        public boolean getUseFrequencyModel(){
            return this.useFreqModel;
        };

        public double[] getExpectedResult() {
            return(this.expectedResults);
        }

        public double[] getExpectedResult(int iTip) {
            return(Arrays.copyOfRange(expectedResults,iTip*nStates,iTip*nStates+nStates));
        }

    };
    Instance test1 = new Instance() {
        //MODEL
        protected final int nCells=2;
        protected final int nStates = (int) (0.5 * (nCells + 1) * (nCells + 2));
        protected final double gamma=0.05;
        protected final double lambda=0.95;
        protected final double mu=0.05;
        protected final boolean useFreqModel = false;
        //EXPECTED
        protected final double [] expectedResults = {
                -0.2, 0.2, 0.0, 0.0, 0.0, 0.0,
                1.0, -2.1, 1.05, 0.05, 0.0, 0.0,
                0.0, 0.1, -0.2, 0.0, 0.1, 0.0,
                0.95, 0.1, 0.0, -2.1, 0.1, 0.95,
                0.0, 0.0, 1.05, 0.05, -2.0999999999999996, 1.0,
                0.0, 0.0, 0.0, 0.0, 0.2, -0.2
        };

        //GETTERS
        public int getNstates(){return this.nStates;};
        public double getStemCellParameter(){
            return this.nCells;
        };
        public double getGammaParameter(){
            return this.gamma;
        };
        public double getLambdaParameter(){
            return this.lambda;
        };
        public double getMuParameter(){
            return this.mu;
        };
        public boolean getUseFrequencyModel(){
            return this.useFreqModel;
        };

        public double[] getExpectedResult() {
            return(this.expectedResults);
        }
    };

    Instance[] all = {test0, test1}; //Add more instances of tests here

    public void testErrorModelClass() {
        for (Instance test : all) {

            Parameter stemCellParameter = new Parameter.Default(1, test.getStemCellParameter());
            Parameter gammaParam = new Parameter.Default(1, test.getGammaParameter());
            Parameter lambdaParam = new Parameter.Default(1, test.getLambdaParameter());
            Parameter muParam = new Parameter.Default(1, test.getMuParameter());
            boolean useFreqModel = test.getUseFrequencyModel();

            int nstates = test.getNstates();
            double[] freqs = new double[nstates];
            Arrays.fill(freqs, (double) 1/nstates);
            freqs[1] = freqs[1]+0.2;
            freqs[2] = freqs[2]-0.15;
            freqs[3] = freqs[3]-0.05;

            AFsequence afseq = new AFsequence(nstates);
            FrequencyModel freqModel = new FrequencyModel(afseq, freqs);

            FlipFlopModel model = new FlipFlopModel("test", stemCellParameter, gammaParam, lambdaParam, muParam, useFreqModel, freqModel);
            model.setupMatrix();

            double[] expectedMat = test.getExpectedResult();
            double[][] mat = model.getAmat();
            double[] flatmat = new double[nstates*nstates];

            for (int k = 0; k < mat.length; ++k) {
                for (int i = 0; i < mat.length; ++i) {
                    flatmat[mat.length*k+i] = mat[k][i];
                }
            }

            for (int j = 0; j < flatmat.length; ++j) {
                assertEquals(flatmat[j], expectedMat[j], 1e-7);
            }
        }
    };
}


