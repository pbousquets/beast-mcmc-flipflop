package dr.evomodel.flipflop;


import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GTR;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import junit.framework.TestCase;

public class TestSubstitutionModelEmpiricalFrequencies extends TestCase {
    interface Instance {
        double [] getPi();
        double [] getRates();

        default DataType getDataType(){return Nucleotides.INSTANCE;}

        default double getDistance(){return 1e5;}

        default double getPrecission(){return 1e-5;}
    }
    Instance test0 = new Instance() {
        double [] pi = {0.1,0.4,0.1,0.4};
        double [] rates = {2,3,4,5,66,7};

        public double [] getPi(){return this.pi;}
        public double [] getRates(){return this.rates;}
    };
    Instance test1 = new Instance() {
        double [] pi = {0.2875,0.3373,0.0688,0.3064};
        double [] rates = {0.42,2.11,0.16,0.28,1,0.34};

        public double [] getPi(){return this.pi;}
        public double [] getRates(){return this.rates;}
    };

    Instance[] all = {test0,test1};

    public void testSubstitutionModel() {
        for (Instance test : all){
            FrequencyModel trueFreqs=new FrequencyModel(test.getDataType(),test.getPi());
            double [] rates=test.getRates();
            Parameter rateACValue= new Parameter.Default(rates[0]);
            Parameter rateAGValue=new Parameter.Default(rates[1]);
            Parameter rateATValue=new Parameter.Default(rates[2]);
            Parameter rateCGValue=new Parameter.Default(rates[3]);
            Parameter rateCTValue=new Parameter.Default(rates[4]);
            Parameter rateGTValue=new Parameter.Default(rates[5]);
            testGTR gtrModel=new testGTR(rateACValue,rateAGValue,rateATValue,rateCGValue,rateCTValue,rateGTValue,trueFreqs);

            double[] mat = new double[test.getDataType().getStateCount() * test.getDataType().getStateCount()];
            gtrModel.getTransitionProbabilities(test.getDistance(), mat);
            double [] calculatedFreqs1=gtrModel.calculateEmpiricalFrequencies();

            for (int k = 0; k < mat.length; ++k) {
                assertEquals(mat[k], trueFreqs.getFrequencies()[k%4], test.getPrecission());
            }

            for (int k = 0; k < test.getDataType().getStateCount(); ++k) {
                assertEquals(calculatedFreqs1[k], trueFreqs.getFrequencies()[k], test.getPrecission());
            }


        }
    }
}

class testGTR extends GTR {
    public testGTR(
            Variable rateACValue,
            Variable rateAGValue,
            Variable rateATValue,
            Variable rateCGValue,
            Variable rateCTValue,
            Variable rateGTValue,
            FrequencyModel freqModel) {
        super(rateACValue, rateAGValue, rateATValue, rateCGValue, rateCTValue, rateGTValue, freqModel);
    }

    public double [] calculateEmpiricalFrequencies(){
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }

        int eigenValPos = 0;
        double minEigenVal=Math.abs(Eval[0]);
        double newEigenVal=0;
        for(int i = 1; i < stateCount; i++){
            newEigenVal=Math.abs(Eval[i]);
            if(newEigenVal < minEigenVal){
                eigenValPos = i;
                minEigenVal = newEigenVal;
            }
        }

        double[] empFreq = new double[stateCount];

        for(int i = 0; i < stateCount; i++){
            empFreq[i] = Evec[i][eigenValPos]*Ievc[eigenValPos][i];

        }
        return empFreq;
    }
}
