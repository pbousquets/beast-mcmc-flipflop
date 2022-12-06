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
import dr.evolution.datatype.DataType;
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

public class FlipFlopModel extends ComplexSubstitutionModel implements Likelihood, Citable {

    /**
     * Constructor
     * @param name              Model name
     * @param stemCellParam     Number of stem cells
     * @param gammaParam        Gamma parameter
     * @param lambdaParam       Lambda parameter
     * @param muParam           Mu parameter
     * @param normalize         Whether to normalize using the out flux of states or not
     * @param useFrequencyModel Whether to use the frequency model or ignore it
     * @param freqModel         Frequency model
     */

    private Parameter stemCellParam;
    private Variable<Double> gammaParam;
    private Variable<Double> lambdaParam;
    private Variable<Double> muParam;
    private double[] dummyStationaryDistribution;
    private int[][] stateVar;
    protected boolean useFrequencyModel;

    public FlipFlopModel(String name, DataType dataType, Parameter stemCellParam, Variable gammaParam, Variable lambdaParam, Variable muParam, boolean normalize, boolean useFrequencyModel, FrequencyModel freqModel) {
        super(name,dataType,freqModel,null);

        addVariable(stemCellParam);
        this.stemCellParam = stemCellParam;
        addVariable(gammaParam);
        this.gammaParam = gammaParam;
        addVariable(lambdaParam);
        this.lambdaParam = lambdaParam;
        addVariable(muParam);
        this.muParam = muParam;
        setNormalization(normalize);
        this.useFrequencyModel=useFrequencyModel;

        if(!useFrequencyModel){
            for (int i=0; i<freqModel.getFrequencyCount();++i){
                freqModel.setFrequency(i,1.0);
            }
            dummyStationaryDistribution = new double[stateCount];
            Arrays.fill(dummyStationaryDistribution, 1);
        }

        generateStateVar();
        updateMatrix = true;
    }

    /*
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        updateMatrix = true;
    }
    */

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if ((model == freqModel) && useFrequencyModel) {
            updateMatrix = true;
            frequenciesChanged();
        } else { // Otherwise, freqModel does not affect likelihood calculation and the super deals with that
            super.handleModelChangedEvent(model, object, index);
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        updateMatrix = true;
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

        // DM we could do all the matrix at once. This is inherited from filling up a reversible transition matrix,
        // In which the upper triangle and the lower triangle are the same but transposed.

        //Upper triangle
        for (i = 0; i < stateCount; i++) {
            for (j = i+1; j < stateCount; j++) {
                amat[i][j] = amat[i][j] * pi[j];
            }
        }

        //Lower triangle
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
        //Update the frequencies if we are not using user-specified/sampled frequencies
        if (!useFrequencyModel){
            double [] frequencies = calculateEmpiricalFrequencies();
            for(int i=0;i<freqModel.getFrequencyCount();i++){
                //freqModel.setFrequency(i,frequencies[i]);
                freqModel.getFrequencyParameter().setParameterValueQuietly(i,frequencies[i]); // we do not need to fire anything
            }
        }

        super.computeStationaryDistribution(); //returns the frequencies
    }

    public double [] calculateEmpiricalFrequencies(){
        synchronized (this) {
            if (updateMatrix) {
                setupMatrix();
            }
        }

        if (!this.wellConditioned) {
            throw new RuntimeException("not well conditioned");
        } else {

            int eigenValPos = -1;
            for (int i = 0; i < stateCount; i++) {
                if (Eval[i] == 0.0) {
                    eigenValPos = i;
                }
            }

            double[] empFreq = new double[stateCount];

            for (int i = 0; i < stateCount; i++) {
                empFreq[i] = Evec[i][eigenValPos] * Ievc[eigenValPos][i];

            }
            return empFreq;
        }
    }

    protected double[] getPi() {
        if (useFrequencyModel) {
            return freqModel.getFrequencies();
        } else {
            return dummyStationaryDistribution;
        }
    }

    public double[][] getAmat(){return amat;}
}