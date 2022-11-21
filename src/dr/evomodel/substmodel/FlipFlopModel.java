/*
 * FlipFlopModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.substmodel;

import cern.colt.matrix.DoubleMatrix1D;
import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import cern.colt.matrix.linalg.Property;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.RobustEigenDecomposition;
import dr.math.matrixAlgebra.RobustSingularValueDecomposition;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;


/**
 * @author Chieh-Hsi Wu
 *
 * A class for the flipflop model
 */

public class FlipFlopModel extends AbstractSubstitutionModel implements Likelihood, Citable {

    /**
     * Constructor
     * @param name              Model name
     * @param stemCellParam     Number of stem cells
     * @param gammaParam        Gamma parameter
     * @param lambdaParam       Lambda parameter
     * @param muParam           Mu parameter
     * @param useFrequencyModel Whether to use the frequency model or ignore it
     * @param freqModel         Frequency model
     */

    protected boolean updateMatrix = true;
    private Parameter stemCellParam;
    private Variable<Double> gammaParam;
    private Variable<Double> lambdaParam;
    private Variable<Double> muParam;
    private double[] dummyStationaryDistribution;
    private int[][] stateVar;

    public FlipFlopModel(String name, Parameter stemCellParam, Variable gammaParam, Variable lambdaParam, Variable muParam, boolean useFrequencyModel, FrequencyModel freqModel, Parameter parameter) {
        super(name, freqModel.getDataType(), freqModel);
        this.infinitesimalRates = parameter;

        rateCount = stateCount * (stateCount - 1);

        if (parameter != null) {
            if (rateCount != infinitesimalRates.getDimension()) {
                throw new RuntimeException("Dimension of '" + infinitesimalRates.getId() + "' ("
                        + infinitesimalRates.getDimension() + ") must equal " + rateCount);
            }
            addVariable(infinitesimalRates);
        }


        stationaryDistribution = new double[stateCount];
        storedStationaryDistribution = new double[stateCount];


        // None of the variables should be null with the current parser, so no need to evaluate them, right?
        addVariable(stemCellParam);
        this.stemCellParam = stemCellParam;
        addVariable(gammaParam);
        this.gammaParam = gammaParam;
        addVariable(lambdaParam);
        this.lambdaParam = lambdaParam;
        addVariable(muParam);
        this.muParam = muParam;
        generateStateVar();

        updateMatrix = true;

        if (useFrequencyModel){
            setNormalization(true);
        } else {
            setNormalization(false);
        }

        dummyStationaryDistribution = new double[stateCount];
        Arrays.fill(dummyStationaryDistribution, 1);

    }

    /*
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        updateMatrix = true;
    }
    */

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == freqModel)
            return; // freqModel only affects the likelihood calculation at the tree root
        super.handleModelChangedEvent(model, object, index);
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        if (!updateMatrix) {
        updateMatrix = true;
//            fireModelChanged();
//        }
    }

    protected void restoreState() {

        // To restore all this stuff just swap the pointers...

        double[] tmp3 = storedEvalImag;
        storedEvalImag = EvalImag;
        EvalImag = tmp3;

        tmp3 = storedStationaryDistribution;
        storedStationaryDistribution = stationaryDistribution;
        stationaryDistribution = tmp3;

        // Inherited
        updateMatrix = storedUpdateMatrix;
        wellConditioned = storedWellConditioned;

        double[] tmp1 = storedEval;
        storedEval = Eval;
        Eval = tmp1;

        double[][] tmp2 = storedIevc;
        storedIevc = Ievc;
        Ievc = tmp2;

        tmp2 = storedEvec;
        storedEvec = Evec;
        Evec = tmp2;

    }

    protected void storeState() {

        storedUpdateMatrix = updateMatrix;

        storedWellConditioned = wellConditioned;

        System.arraycopy(stationaryDistribution, 0, storedStationaryDistribution, 0, stateCount);
        System.arraycopy(EvalImag, 0, storedEvalImag, 0, stateCount);

        // Inherited
        System.arraycopy(Eval, 0, storedEval, 0, stateCount);
        for (int i = 0; i < stateCount; i++) {
            System.arraycopy(Ievc[i], 0, storedIevc[i], 0, stateCount);
            System.arraycopy(Evec[i], 0, storedEvec[i], 0, stateCount);
        }
    }

    public void getTransitionProbabilities(double distance, double[] matrix) {

        double temp;

        int i, j, k;

        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }

        if (!wellConditioned) {
            Arrays.fill(matrix, 0.0);
            return;
        }


// Eigenvalues and eigenvectors of a real matrix A.
//
// If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is diagonal
// and the eigenvector matrix V is orthogonal. I.e. A = V D V^t and V V^t equals
// the identity matrix.
//
// If A is not symmetric, then the eigenvalue matrix D is block diagonal with
// the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
// lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda]. The columns
// of V represent the eigenvectors in the sense that A*V = V*D. The matrix
// V may be badly conditioned, or even singular, so the validity of the
// equation A = V D V^{-1} depends on the conditioning of V.

        double[][] iexp = popiexp();

        for (i = 0; i < stateCount; i++) {

            if (EvalImag[i] == 0) {
                // 1x1 block
                temp = Math.exp(distance * Eval[i]);
                for (j = 0; j < stateCount; j++) {
                    iexp[i][j] = Ievc[i][j] * temp;
                }
            } else {
                // 2x2 conjugate block
                // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                int i2 = i + 1;
                double b = EvalImag[i];
                double expat = Math.exp(distance * Eval[i]);
                double expatcosbt = expat * Math.cos(distance * b);
                double expatsinbt = expat * Math.sin(distance * b);

                for (j = 0; j < stateCount; j++) {
                    iexp[i][j] = expatcosbt * Ievc[i][j] + expatsinbt * Ievc[i2][j];
                    iexp[i2][j] = expatcosbt * Ievc[i2][j] - expatsinbt * Ievc[i][j];
                }
                i++; // processed two conjugate rows
            }
        }

        int u = 0;
        for (i = 0; i < stateCount; i++) {
            for (j = 0; j < stateCount; j++) {
                temp = 0.0;
                for (k = 0; k < stateCount; k++) {
                    temp += Evec[i][k] * iexp[k][j];
                }
                if (temp < 0.0)
                    matrix[u] = minProb;
                else
                    matrix[u] = temp;
                u++;
            }
        }
        pushiexp(iexp);
    }


    //store rate matrix
    public void storeIntoAmat(){
        int S = (int) stemCellParam.getParameterValue(0);
        double gamma = gammaParam.getValue(0);
        double lam = lambdaParam.getValue(0);
        double mu = muParam.getValue(0);

        for (int down = 0; down < stateCount; down++){ //Represents the new state (out-state)
            int k_down = stateVar[down][0];
            int m_down = stateVar[down][1];
            for (int across = 0; across < stateCount; across++) { //Represents the old state (in-state)
                int k = stateVar[across][0];
                int m = stateVar[across][1];
                if (k == k_down-1 & m == m_down){
                    amat[across][down] = (S-m-k) * (k*lam/(S-1) + 2*mu);
                } else if (k == k_down & m == m_down-1){
                    amat[across][down] = m * (S-m-k) * lam / (S-1);
                } else if (k == k_down+1 & m == m_down-1) {
                    amat[across][down] = k * (m*lam/(S-1)+mu);
                } else if (k == k_down+1 & m == m_down) {
                    amat[across][down] = k * ((S-m-k)*lam/(S-1)+gamma);
                } else if (k == k_down & m == m_down+1) {
                    amat[across][down] = m * (S-m-k) * lam / (S-1);
                } else if (k == k_down-1 & m == m_down+1) {
                    amat[across][down] = m * (k*lam/(S-1)+2*gamma);
                } else if (k == k_down & m == m_down) {
                    amat[across][down] = -(2*((k+m)*(S-m-k)+k*m)*lam/(S-1) + (k+2*m)*gamma + (2*S-(k+2*m))*mu);
                } else {
                    amat[across][down] = 0;
                }
            }
        }

        double[] pi = getPi();

        int i, j;
        for (i = 0; i < stateCount; i++) {
            for (j = i+1; j < stateCount; j++) {
                amat[i][j] = amat[i][j] * pi[j];
            }
        }

        // Copy lower triangle in column-order form (transposed)
        for (j = 0; j< stateCount; j++) {
            for (i = j+1; i < stateCount; i++) {
                amat[i][j] = amat[i][j] * pi[j];
            }
        }
    }



    /**
     * This method will create a matrix with all possible cell states, for parameters k (one met allele) and m (both met alleles)
     * For S=2, possible states are for (k,m): 0,0 0,1 0,2 1,0, 2,0 1,1
     **/
    public void generateStateVar(){
        int S = (int) stemCellParam.getParameterValue(0);
        int num_states = stateCount;

        this.stateVar = new int[num_states][];
        int ii = 0;
        for (int m = 0; m <= S; m++){
            for (int k = 0; k <= S; k++){
                if (k+m <=S) {
                    this.stateVar[ii] = new int[]{k, m};
                    ii++;
                }
            }
        }
    }

    protected void computeStationaryDistribution() {
        if (doNormalization) {
            stationaryDistribution = freqModel.getFrequencies();
        } else {
            stationaryDistribution = dummyStationaryDistribution;
        }
    }

    protected double[] getPi() {
        if (doNormalization) {
            return freqModel.getFrequencies();
        } else {
            return dummyStationaryDistribution;
        }
    }

    public double[][] getAmat(){return amat;}

    public double[] getStationaryDistribution() {
        return stationaryDistribution;
    }


    public void setupMatrix() {

        if (!eigenInitialised) {
            initialiseEigen();
            storedEvalImag = new double[stateCount];
        }

        int i = 0;

        storeIntoAmat();


        makeValid(amat, stateCount);

        // compute eigenvalues and eigenvectors

        RobustEigenDecomposition eigenDecomp;
        try {
            eigenDecomp = new RobustEigenDecomposition(new DenseDoubleMatrix2D(amat), maxIterations);
        } catch (ArithmeticException ae) {
            System.err.println(ae.getMessage());
            wellConditioned = false;
            System.err.println("amat = \n" + new Matrix(amat));
            return;
        }

        DoubleMatrix2D eigenV = eigenDecomp.getV();
        DoubleMatrix1D eigenVReal = eigenDecomp.getRealEigenvalues();
        DoubleMatrix1D eigenVImag = eigenDecomp.getImagEigenvalues();
        DoubleMatrix2D eigenVInv;

        // A better (?) approach to checking diagonalizability comes from:
        //
        // J. Gentle (2007) Matrix Algebra
        //
        // Diagonalizbility Theorem: A matrix A is (complex) diagonalizable iff all distinct eigenvalues \lambda_l
        // with algebraic multiplicity m_l are semi-simple, i.e.
        //
        //          rank(A - \lambda_l I) = n - m_l
        //
        // Equivalently (?), eigenV must be non-singular.
        //
        // SVD is needed to numerically approximate the rank of a matrix, so we can check Algrebra.rank()
        // or Algebra.cond() with almost equal amounts of work.  I don't know which is more reliable. -- MAS

        if (checkConditioning) {
            RobustSingularValueDecomposition svd;
            try {
                svd = new RobustSingularValueDecomposition(eigenV, maxIterations);
            } catch (ArithmeticException ae) {
                System.err.println(ae.getMessage());
                wellConditioned = false;
                return;
            }
            if (svd.cond() > maxConditionNumber) {
                wellConditioned = false;
                return;
            }
        }

        try {
            eigenVInv = alegbra.inverse(eigenV);
        } catch (IllegalArgumentException e) {
            wellConditioned = false;
            return;
        }

        Ievc = eigenVInv.toArray();
        Evec = eigenV.toArray();
        Eval = eigenVReal.toArray();
        EvalImag = eigenVImag.toArray();

        // Check for valid decomposition
        for (i = 0; i < stateCount; i++) {
            if (Double.isNaN(Eval[i]) || Double.isNaN(EvalImag[i]) ||
                    Double.isInfinite(Eval[i]) || Double.isInfinite(EvalImag[i])) {
                wellConditioned = false;
                return;
            } else if (Math.abs(Eval[i]) < 1e-10) {
                Eval[i] = 0.0;
            }
        }

        updateMatrix = false;
        wellConditioned = true;
        // compute normalization and rescale eigenvalues

        computeStationaryDistribution();

        if (doNormalization) {
            double subst = 0.0;

            for (i = 0; i < stateCount; i++)
                subst += -amat[i][i] * stationaryDistribution[i];

            for (i = 0; i < stateCount; i++) {
                Eval[i] /= subst;
                EvalImag[i] /= subst;
            }
        }
    }



    protected void frequenciesChanged() {
    }

    protected void ratesChanged() {
    }

    protected void setupRelativeRates() {
    }

    protected Parameter infinitesimalRates;

    public LogColumn[] getColumns() {
        return new LogColumn[]{
                new FlipFlopModel.LikelihoodColumn(getId())
        };
    }

    protected class LikelihoodColumn extends NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    public void setMaxIterations(int max) {
        maxIterations = max;
    }


    private double[] stationaryDistribution;
    private double[] storedStationaryDistribution;

    protected boolean doNormalization = true;

    protected double[] EvalImag;
    protected double[] storedEvalImag;

    protected boolean wellConditioned = true;
    private boolean storedWellConditioned;

    protected static final double minProb = Property.DEFAULT.tolerance();
    private static final Algebra alegbra = new Algebra(minProb);

    private double maxConditionNumber = 1000;

    private int maxIterations = 1000;

    private boolean checkConditioning = true;

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(probability,this))
            return 0;
        return Double.NEGATIVE_INFINITY;
    }

    /**
     * Needs to be evaluated before the corresponding data likelihood.
     * @return
     */
    public boolean evaluateEarly() {
        return true;
    }

    public String prettyName() {
        return Abstract.getPrettyName(this);
    }

    public void setNormalization(boolean doNormalization) {
        this.doNormalization = doNormalization;
    }

    public void makeDirty() {

    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<Likelihood>(Arrays.asList(this));
    }

    @Override
    public boolean isUsed() {
        return super.isUsed() && isUsed;
    }

    public void setUsed() {
        isUsed = true;
    }

    private boolean isUsed = false;

    private double[] probability = null;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Complex-diagonalizable, irreversible substitution model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.EDWARDS_2011_ANCIENT);
    }
}