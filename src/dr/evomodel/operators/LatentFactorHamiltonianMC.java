package dr.evomodel.operators;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.model.*;
import dr.inference.operators.AbstractHamiltonianMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import java.util.Random;

/**
 * Created by max on 12/2/15.
 */
public class LatentFactorHamiltonianMC extends AbstractHamiltonianMCOperator{
    private LatentFactorModel lfm;
    private FullyConjugateMultivariateTraitLikelihood tree;
    private MatrixParameterInterface factors;
    private MatrixParameterInterface loadings;
    private MatrixParameterInterface Precision;
    private int nfac;
    private int ntaxa;
    private int ntraits;
    private double stepSize;
    private int nSteps;
    private boolean diffusionSN=true;
    private Random random;


    public LatentFactorHamiltonianMC(LatentFactorModel lfm, FullyConjugateMultivariateTraitLikelihood tree, double weight, CoercionMode mode, double stepSize, int nSteps, double momentumSd){
        super(mode, momentumSd);
        setWeight(weight);
        this.lfm = lfm;
        this.tree = tree;
        this.factors = lfm.getFactors();
        this.loadings = lfm.getLoadings();
        this.Precision = lfm.getColumnPrecision();
        nfac = lfm.getFactorDimension();
        ntaxa = lfm.getFactors().getColumnDimension();
        ntraits = Precision.getRowDimension();
        this.stepSize = stepSize;
        this.nSteps = nSteps;
        random = new Random(555);
    }



    @Override
    public double getCoercableParameter() {
        return 0;
    }

    @Override
    public void setCoercableParameter(double value) {

    }

    @Override
    public double getRawParameter() {
        return 0;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Latent Factor Hamiltonian Monte Carlo";
    }

    private double[] getMatrix(int element, double[] residual){
        double answer[] = new double[this.nfac];
        for (int i = 0; i < this.nfac ; i++) {
            for (int j = 0; j < ntraits; j++) {
                answer[i] += loadings.getParameterValue(j, i) * Precision.getParameterValue(j, j) *
                        residual[j * ntaxa + element];
            }
        }
        return answer;
    }

    private double[] getGradient(int randel, double[] mean, double[][] prec, double precfactor){
        double[] residual = lfm.getResidual();
        double[] derivative = getMatrix(randel, residual);

        if(diffusionSN){
            for (int i = 0; i <mean.length ; i++) {
                derivative[i] -= (factors.getParameterValue(i, randel) - mean[i])*precfactor;
            }
        }
        else{
            for (int i = 0; i <mean.length ; i++) {
                double sumi = 0;
                for (int j = 0; j <mean.length ; j++) {
                    sumi += prec[i][j]*(factors.getParameterValue(j, randel) - mean[j]);
                }
                derivative[i] -= sumi;
            }
        }
        return derivative;
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        int randel = MathUtils.nextInt(ntaxa);
//        System.out.println(randel);



        double[] mean = tree.getConditionalMean(randel);
        double precfactor = 0;
        double[][] prec = null;
        double rand = random.nextDouble();
//        System.out.println(rand);
        double functionalStepSize = stepSize * rand;

        if(diffusionSN){
            precfactor = tree.getPrecisionFactor(randel);
        }
        else {
            prec = tree.getConditionalPrecision(randel);
        }

        double[] derivative = getGradient(randel, mean, prec, precfactor);
        drawMomentum(lfm.getFactorDimension());

        double prop=0;
        for (int i = 0; i < momentum.length ; i++) {
            prop += momentum[i]*momentum[i] / (2 * getMomentumSd()*getMomentumSd());
        }

        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            momentum[i] = momentum[i] - functionalStepSize / 2 * derivative[i];
        }

        for (int i = 0; i <nSteps ; i++) {
            for (int j = 0; j <lfm.getFactorDimension() ; j++) {
                factors.setParameterValueQuietly(j, randel, factors.getParameterValue(j, randel) + functionalStepSize * momentum[j] / (getMomentumSd() * getMomentumSd()));
            }
//            System.out.println("randel");
//            System.out.println(randel);
            if(factors instanceof FastMatrixParameter) {
                for (int j = 0; j <factors.getParameter(randel).getDimension() ; j++) {
                    factors.fireParameterChangedEvent(randel * factors.getRowDimension() + j, null);
                }
//                factors.fireParameterChangedEvent();
            }
            else{
                for (int j = 0; j <factors.getParameter(randel).getDimension() ; j++) {
                    factors.getParameter(randel).fireParameterChangedEvent(j, null);
                }
//                factors.getParameter(randel).fireParameterChangedEvent();
            }

//            factors.fireParameterChangedEvent();

            if(i!=nSteps){
                derivative=getGradient(randel, mean, prec, precfactor);

                for (int j = 0; j < lfm.getFactorDimension() ; j++) {
                    momentum[j] = momentum[j] - functionalStepSize * derivative[j];
                }
            }
        }

        derivative=getGradient(randel , mean , prec , precfactor);
        for (int i = 0; i <lfm.getFactorDimension() ; i++) {

            momentum[i] = momentum[i] - functionalStepSize / 2 * derivative[i];
        }

        double res=0;
        for (int i = 0; i <momentum.length ; i++) {
            res += momentum[i] * momentum[i] / (2 * getMomentumSd() * getMomentumSd());
        }
        return prop-res;
    }
}
