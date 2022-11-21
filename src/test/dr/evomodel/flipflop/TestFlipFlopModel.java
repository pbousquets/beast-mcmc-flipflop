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

    class Test {
        public Test(int age, int nCells, double gamma, double lambda, double mu, boolean useFreqModel, double[] stationaryDistribution, double[] expectedRateMat, double[] expectedTransitionMat){
            this.age = age;
            this.nCells = nCells;
            this.gamma = gamma;
            this.lambda = lambda;
            this.mu = mu;
            this.useFreqModel = useFreqModel;
            this.stationaryDistribution = stationaryDistribution;
            this.nStates = (int) (0.5 * (nCells + 1) * (nCells + 2));
            this.expectedRateMat = expectedRateMat;
            this.expectedTransitionMat = expectedTransitionMat;

            AFsequence afseq = new AFsequence(this.nStates);
            FrequencyModel freqModel = new FrequencyModel(afseq, this.stationaryDistribution);

            Parameter stemCellParameter = new Parameter.Default(1, nCells);
            Parameter gammaParam = new Parameter.Default(1, gamma);
            Parameter lambdaParam = new Parameter.Default(1, lambda);
            Parameter muParam = new Parameter.Default(1, mu);

            model = new FlipFlopModel("test", stemCellParameter, gammaParam, lambdaParam, muParam, useFreqModel, freqModel, null);
        }

        public int getNstates(){return nStates;};

        public double[] getExpectedRateMat(){
            return expectedRateMat;
        };

        public double[] getExpectedTransitionMat(){
            return expectedTransitionMat;
        }
        public FlipFlopModel getModel(){
            return(model);
        };

        protected final int age;
        protected final int nCells;
        protected final int nStates;
        protected final double gamma;
        protected final double lambda;
        protected final double mu;
        protected final boolean useFreqModel;
        protected final double[] stationaryDistribution;
        protected final double[] expectedRateMat;
        protected final double[] expectedTransitionMat;

        protected final FlipFlopModel model;
    };

    // TEST1
    double[] stationaryDistributionTest = new double[]{0.16666666666666666, 0.3666666666666667, 0.016666666666666663, 0.11666666666666665, 0.16666666666666666, 0.16666666666666666};
    double[] expectedRateMatTest1 = {
            -0.07333333333333335, 0.07333333333333335, 0.0, 0.0, 0.0, 0.0,
            0.16666666666666666, -0.18999999999999997, 0.017499999999999998, 0.005833333333333333, 0.0, 0.0,
            0.0, 0.036666666666666674, -0.053333333333333344, 0.0, 0.016666666666666666, 0.0,
            0.15833333333333333, 0.036666666666666674, 0.0, -0.37, 0.016666666666666666, 0.15833333333333333,
            0.0, 0.0, 0.017499999999999998, 0.005833333333333333, -0.19, 0.16666666666666666,
            0.0, 0.0, 0.0, 0.0, 0.03333333333333333, -0.03333333333333333
    };

    double[] expectedTransitionMatTest1 = {
            0.69328497, 0.26284929, 0.03255818, 0.00387858, 0.00228277, 0.00514622
           , 0.60575896, 0.29333402, 0.07214965, 0.00539943, 0.0069886 , 0.01636934
           , 0.15752991, 0.15115086, 0.54539944, 0.00305908, 0.06098829, 0.08187243
           , 0.3619755 , 0.13927455, 0.02489913, 0.01371788, 0.06692113, 0.3932118
           , 0.01886422, 0.01496843, 0.06401687, 0.00424752, 0.19672447, 0.70117849
           , 0.00299875, 0.00257297, 0.01658318, 0.0020206 , 0.13831613, 0.83750837
    };


    Test test1 = new Test(2, 2, 0.05, 0.95, 0.05, true, stationaryDistributionTest, expectedRateMatTest1, expectedTransitionMatTest1);

    // TEST2
    double[] expectedRateMatTest2 = {
            -0.2, 0.2, 0.0, 0.0, 0.0, 0.0,
            1.0, -2.1, 1.05, 0.05, 0.0, 0.0,
            0.0, 0.1, -0.2, 0.0, 0.1, 0.0,
            0.95, 0.1, 0.0, -2.1, 0.1, 0.95,
            0.0, 0.0, 1.05, 0.05, -2.0999999999999996, 1.0,
            0.0, 0.0, 0.0, 0.0, 0.2, -0.2
    };
    double[] expectedTransitionMatTest2 = {
           0.78525757, 0.08144377, 0.12179875, 0.00193187, 0.00463868, 0.00492936
           , 0.41639524, 0.07340708, 0.45193286, 0.00316781, 0.0227272 , 0.03236981
           , 0.05993344, 0.04304122, 0.7921188 , 0.00193187, 0.04304122, 0.05993344
           , 0.38759567, 0.04304122, 0.12179875, 0.01692745, 0.04304122, 0.38759567
           , 0.03236981, 0.0227272 , 0.45193286, 0.00316781, 0.07340708, 0.41639524
           , 0.00492936, 0.00463868, 0.12179875, 0.00193187, 0.08144377, 0.78525757
    };
    Test test2 = new Test(2, 2, 0.05, 0.95, 0.05, false, stationaryDistributionTest, expectedRateMatTest2, expectedTransitionMatTest2);



    Test[] all = {test1, test2}; //Add more tests here

    public void testRateMatrix() {
        for (Test test : all) {
            int nstates = test.getNstates();
            test.model.setupMatrix();
            double[][] rateMat = test.model.getAmat();
            double[] expectedMat = test.getExpectedRateMat();
            double[] flatmat = new double[nstates*nstates];

            for (int k = 0; k < rateMat.length; ++k) {
                for (int i = 0; i < rateMat.length; ++i) {
                    flatmat[rateMat.length*k+i] = rateMat[k][i];
                }
            }

            for (int j = 0; j < flatmat.length; ++j) {
                assertEquals(expectedMat[j], flatmat[j], 1e-7);
            }
        }
    };

    public void testTransitionMatrix(){
        for (Test test : all) {
            int nstates = test.getNstates();
            double[] transitionMat = new double[nstates*nstates];

            test.model.getTransitionProbabilities(test.age, transitionMat);
            double[] expectedMat = test.getExpectedTransitionMat();

            for (int j = 0; j < transitionMat.length; ++j) {
                assertEquals(expectedMat[j], transitionMat[j], 1e-7);
            }
        }
    }
}


