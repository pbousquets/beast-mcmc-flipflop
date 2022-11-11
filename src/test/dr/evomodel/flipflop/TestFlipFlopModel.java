package test.dr.evomodel.flipflop;

import dr.evolution.datatype.AFsequence;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.lang.reflect.Array;
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

        String getSequenceString (int iTip);

        double[] getExpectedResult();

    }

    Instance test0 = new Instance() {
        //DATA
        protected final double[] sequence={0.1, 0.4, 0.14, 0.12, 0.99, .05, 0.01};

        //MODEL
        protected final int nCells=2;
        protected final int nStates = (int) (0.5 * (nCells + 1) * (nCells + 2));
        protected final double gamma=0.1;
        protected final double lambda=1;
        protected final double mu=0.1;

        //EXPECTED
        protected final double [] expectedResults = {
                -0.4, 0.4, 0.0, 0.0, 0.0, 0.0,
                1.1, -2.4000000000000004, 1.2, 0.1, 0.0, 0.0,
                0.0, 0.2, -0.4, 0.0, 0.2, 0.0,
                1.0, 0.2, 0.0, -2.4000000000000004, 0.2, 1.0,
                0.0, 0.0, 1.2, 0.1, -2.4, 1.1,
                0.0, 0.0, 0.0, 0.0, 0.4, -0.4
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
        public double[] getSequence(){
            return this.sequence;
        };

        public String getSequenceString(int iTip){
            return(String.valueOf(this.sequence[iTip])); //Implementation for only one patter/locus
        }

        public double[] getExpectedResult() {
            return(this.expectedResults);
        }

        public double[] getExpectedResult(int iTip) {
            return(Arrays.copyOfRange(expectedResults,iTip*nStates,iTip*nStates+nStates));
        }

    };

    Instance[] all = {test0}; //Add more instances of tests here

    public void testErrorModelClass() {
        for (Instance test : all) {

            Parameter stemCellParameter = new Parameter.Default(1, test.getStemCellParameter());
            Parameter gammaParam = new Parameter.Default(1, test.getGammaParameter());
            Parameter lambdaParam = new Parameter.Default(1, test.getLambdaParameter());
            Parameter muParam = new Parameter.Default(1, test.getMuParameter());
            int nstates = test.getNstates();
            double[] freqs = new double[nstates];
            for(int i = 0; i < freqs.length; i++){
                freqs[i] = 1.0/nstates;
            }
            AFsequence afseq = new AFsequence(nstates);
            FrequencyModel freqModel = new FrequencyModel(afseq, freqs);

            FlipFlopModel model = new FlipFlopModel("test", stemCellParameter, gammaParam, lambdaParam, muParam, freqModel);
            model.setupMatrix();

            double[] expectedMat = test.getExpectedResult();
            double[][] mat = model.getAmat();
            double[] flatmat = new double[nstates*nstates];

            for (int k = 0; k < mat.length; ++k) {
                for (int i = 0; i < mat.length; ++i) {
                    flatmat[mat.length*k+i] = mat[k][i];
                }
            }

            //System.out.println(Arrays.deepToString(mat));
            //System.out.println(Arrays.toString(flatmat));
            //System.out.println(Arrays.toString(expectedMat));

            for (int j = 0; j < flatmat.length; ++j) {
                assertEquals(flatmat[j], expectedMat[j], 1e-10);
            }
        }
    };

}
