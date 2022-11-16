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

import dr.inference.model.Variable;
import dr.inference.model.Parameter;
import dr.evolution.datatype.AFsequence;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Diego Mallo
 */
public class FlipFlopModelParser extends AbstractXMLObjectParser {

    public static final String FLIPFLOP_MODEL = "flipflopModel";
    public static final String STEM_CELLS = "stemCells";
    public static final String GAMMA = "gamma";
    public static final String LAMBDA= "lambda";
    public static final String MU= "mu";
    public static final String USEFREQMODEL= "useFrequencyModel";

    public String getParserName() {
        return FLIPFLOP_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter stemCellParam = (Parameter) xo.getElementFirstChild(STEM_CELLS);
        Variable gammaParam = (Variable) xo.getElementFirstChild(GAMMA);
        Variable lambdaParam = (Variable) xo.getElementFirstChild(LAMBDA);
        Variable muParam = (Variable) xo.getElementFirstChild(MU);
        boolean useFrequencyModel = xo.getAttribute(USEFREQMODEL, true);

        int S = (int) stemCellParam.getParameterValue(0);
        int stateCount = (int) (0.5*(S+1)*(S+2));
        double[] freqs = new double[stateCount];
        for(int i = 0; i < freqs.length; i++){
            freqs[i] = 1.0/stateCount;
        }

        AFsequence afseq = new AFsequence(stateCount);
        FrequencyModel freqModel = new FrequencyModel(afseq, freqs);

        Logger.getLogger("dr.evomodel").info("\n---\n\nCreating FlipFlip model");
        Logger.getLogger("dr.evomodel").info("  - Initial gamma = " + gammaParam.getValue(0));
        Logger.getLogger("dr.evomodel").info("  - Initial lambda = " + lambdaParam.getValue(0));
        Logger.getLogger("dr.evomodel").info("  - Initial mu = " + muParam.getValue(0));
        Logger.getLogger("dr.evomodel").info("\n---");

        return new FlipFlopModel(xo.getId(), stemCellParam, gammaParam, lambdaParam, muParam, useFrequencyModel, freqModel);
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
            AttributeRule.newBooleanRule(USEFREQMODEL, true)
    };
}
