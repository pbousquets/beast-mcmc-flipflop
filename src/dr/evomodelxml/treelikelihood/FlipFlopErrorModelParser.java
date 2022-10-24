/*
 * FlipFlopErrorModelParser.java
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


package dr.evomodelxml.treelikelihood;

import dr.evolution.util.TaxonList;
import dr.evomodel.treelikelihood.FlipFlopErrorModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 */
public class FlipFlopErrorModelParser extends AbstractXMLObjectParser {

    public static final String AFSEQUENCE_ERROR_MODEL = "AFsequenceErrorModel";
    public static final String STEM_CELLS = "stemCells";
    public static final String DELTA_OFFSET = "deltaOffset";
    public static final String ETA_OFFSET = "etaOffset";
    public static final String KAPPA_SCALE = "kappaScale";
    public static final String EXCLUDE = "exclude";
    public static final String INCLUDE = "include";

    public String getParserName() {
        return AFSEQUENCE_ERROR_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter deltaParameter = null;
        if (xo.hasChildNamed(DELTA_OFFSET)) {
            deltaParameter = (Parameter) xo.getElementFirstChild(DELTA_OFFSET);
        }

        Parameter etaParameter = null;
        if (xo.hasChildNamed(ETA_OFFSET)) {
            etaParameter = (Parameter) xo.getElementFirstChild(ETA_OFFSET);
        }

        Parameter kappaParameter = null;
        if (xo.hasChildNamed(KAPPA_SCALE)) {
            kappaParameter = (Parameter) xo.getElementFirstChild(KAPPA_SCALE);
        }

        Parameter stemCellParameter = null;
        if (xo.hasChildNamed(STEM_CELLS)) {
            stemCellParameter = (Parameter) xo.getElementFirstChild(STEM_CELLS);
            if ((stemCellParameter.getParameterValue(0)) % 1 != 0){
                throw new XMLParseException("The stem cell parameter must be an integer. Current value: " + stemCellParameter.getParameterValue(0));
            }
        }

        TaxonList includeTaxa = null;
        TaxonList excludeTaxa = null;

        if (xo.hasChildNamed(INCLUDE)) {
            includeTaxa = (TaxonList) xo.getElementFirstChild(INCLUDE);
        }

        if (xo.hasChildNamed(EXCLUDE)) {
            excludeTaxa = (TaxonList) xo.getElementFirstChild(EXCLUDE);
        }

        FlipFlopErrorModel afModel = new FlipFlopErrorModel(includeTaxa, excludeTaxa, stemCellParameter, deltaParameter, etaParameter, kappaParameter);

        Logger.getLogger("dr.evomodel").info("Using allele frequency sequence error model with " + (int) stemCellParameter.getParameterValue(0) + " stem cells");

        return afModel;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element returns a model that allows for allele frequency rather than discrete states.";
    }

    public Class getReturnType() {
        return FlipFlopErrorModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(STEM_CELLS, Parameter.class, "Number of stem cells"),
            new ElementRule(DELTA_OFFSET, Parameter.class, "Delta prior, regarding the error model", true),
            new ElementRule(ETA_OFFSET, Parameter.class, "Eta prior, regarding the error model", true),
            new ElementRule(KAPPA_SCALE, Parameter.class, "Kappa prior, regarding the scale parameter of the AF peak", true),
            new XORRule(
                    new ElementRule(INCLUDE, TaxonList.class, "A set of taxa to which to apply the damage model to"),
                    new ElementRule(EXCLUDE, TaxonList.class, "A set of taxa to which to not apply the damage model to"),
                    true)
    };
}
