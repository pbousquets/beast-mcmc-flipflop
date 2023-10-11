/*
 * FrequencyModel.java
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

import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;

/**
 * A model of equilibrium frequencies that depends only on a parameter
 *
 * @author Diego Mallo
 */
public class FlipFlopCenancestorFrequencyProportionModel extends AbstractGeneralFrequencyModel {


    public FlipFlopCenancestorFrequencyProportionModel(DataType dataType, Parameter methylatedProportionParam) {
        super(dataType);
        dim = dataType.getStateCount();

        this.methylatedProportionParameter = methylatedProportionParam;
        addVariable(methylatedProportionParam);

        this.dataType = dataType;
    }

    public void setFrequency(int i, double value) {
        if (i==0){
           methylatedProportionParameter.setParameterValue(0,1-value);
        }
        if (i==dim-1) {
            methylatedProportionParameter.setParameterValue(0,value);
        } else if(Math.abs(value) > Math.ulp(value)) { // These values can only be 0s
            throw new RuntimeException("");
        }
    }

    public double getFrequency(int i) {
        if(i==0){
            return 1-methylatedProportionParameter.getParameterValue(0);
        } else if (i==dim-1){
            return methylatedProportionParameter.getParameterValue(0);
        }
        return 0;
    }

    public int getFrequencyCount() {
        return dim;
    }

    public double[] getFrequencies() {
        double [] frequencies = new double [dim];
        frequencies[0]=1-methylatedProportionParameter.getParameterValue(0);
        frequencies[dim-1]= methylatedProportionParameter.getParameterValue(0);
        return frequencies;
    }

    private DataType dataType = null;
    private int dim;
    Parameter methylatedProportionParameter = null;
}
