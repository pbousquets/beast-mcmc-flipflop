/*
 * FlipFlopModelParser.java
 *
 * By DM, modified from HKYParser
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

package dr.evomodelxml.substmodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.inference.model.Variable;
import dr.inference.model.Parameter;
import dr.evolution.datatype.AFsequence;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.xml.*;

import java.text.NumberFormat;
import java.util.logging.Logger;
import java.util.regex.Pattern;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Diego Mallo
 */
public class FlipFlopModelParser extends AbstractXMLObjectParser {

    public static final String FLIPFLOP_MODEL = "flipflopModel";
    public static final String STEM_CELLS = "stemCells";
    public static final String GAMMA = "gamma";
    public static final String LAMBDA = "lambda";
    public static final String MU = "mu";
    public static final String NORMALIZE = "normalize";
    public static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return FLIPFLOP_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter stemCellParam = (Parameter) xo.getElementFirstChild(STEM_CELLS);
        Variable gammaParam = (Variable) xo.getElementFirstChild(GAMMA);
        Variable lambdaParam = (Variable) xo.getElementFirstChild(LAMBDA);
        Variable muParam = (Variable) xo.getElementFirstChild(MU);
        boolean normalize = xo.getAttribute(NORMALIZE, false);

        PatternList data = (PatternList) xo.getChild(PatternList.class);
        DataType dataType = data.getDataType();
        int stateCount = dataType.getStateCount();

        FrequencyModel freqModel=null;
        boolean useFrequencyModel = false;

        StringBuilder sb = new StringBuilder("\n---\n\nCreating FlipFlip model");
        sb.append("\n\t- Initial gamma = " + gammaParam.getValue(0));
        sb.append("\n\t- Initial lambda = " + lambdaParam.getValue(0));
        sb.append("\n\t- Initial mu = " + muParam.getValue(0));
        sb.append("\n\t- Normalization using flux of out-states: " + normalize);

        boolean updateDimension = false;

        if (xo.hasChildNamed(FREQUENCIES)) { //Given frequency parameter. We normalize them by default
            Parameter freqsParam = (Parameter) xo.getElementFirstChild(FREQUENCIES);

            double cumfreqs=0;

            if(freqsParam.getDimension()!=stateCount){
                freqsParam.setDimension(stateCount);
                updateDimension = true;
            } else {
                for (int i=0; i<freqsParam.getDimension();i++){
                    cumfreqs += freqsParam.getParameterValue(i);
                }
            }

            for (int i = 0; i < freqsParam.getDimension(); i++) {
                if (cumfreqs != 0)
                    freqsParam.setParameterValue(i, freqsParam.getParameterValue(i) / cumfreqs);
                else
                    freqsParam.setParameterValue(i, 1.0 / freqsParam.getDimension());
            }

            freqModel = new FrequencyModel(dataType, freqsParam);
            useFrequencyModel = true;
        } else { //No frequency parameter given. We just make one with equal freqs. It will not be used
            // to specify the transition matrix, but they will be used to calculate the stationary freqs
            double[] freqs = new double[stateCount];
            for (int i = 0; i < freqs.length; i++) {
                freqs[i] = 1.0 / stateCount;
            }
            freqModel = new FrequencyModel(dataType, freqs);
        }

        sb.append("\n\t- Using a frequency model to set the transition matrix: " + useFrequencyModel);
        if (updateDimension){
            sb.append("\n\t- WARNING: XML frequency parameter of improper dimension. Fixed.");
        }
        sb.append("\n\t- Initial frequency model: {");
        NumberFormat format = NumberFormat.getNumberInstance();
        format.setMaximumFractionDigits(5);
        sb.append(format.format(freqModel.getFrequencyParameter().getParameterValue(0)));
        for (int i = 1; i < stateCount; i++) {
            sb.append(", ");
            sb.append(format.format(freqModel.getFrequencyParameter().getParameterValue(i)));
        }
        sb.append("}");
        sb.append("\n---");
        Logger.getLogger("dr.evomodel").info(sb.toString());

        return new FlipFlopModel(xo.getId(), dataType, stemCellParam, gammaParam, lambdaParam, muParam, normalize, useFrequencyModel, freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents an instance of the flipflop model (Gabbutt, 2022), depicting the dynamics of fluctuating CpG sites.";
    }

    public Class getReturnType() {
        return FlipFlopModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GAMMA,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
            new ElementRule(LAMBDA,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
            new ElementRule(MU,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
            new ElementRule(STEM_CELLS,
                    new XMLSyntaxRule[]{new ElementRule(Variable.class)}),
            new ElementRule(PatternList.class,false),
            new ElementRule(FREQUENCIES, Parameter.class,"Equilibrium frequencies",true),
            AttributeRule.newBooleanRule(NORMALIZE, true)
    };
}
