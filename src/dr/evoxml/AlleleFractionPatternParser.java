/*
 * AlleleFractionPatternParser.java
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

package dr.evoxml;

import dr.xml.*;
import dr.evolution.util.Taxa;
import dr.evolution.datatype.AlleleFraction;
import dr.evolution.alignment.AFalignment;
import dr.evolution.datatype.DataType;
import dr.evoxml.util.DataTypeUtils;
import dr.evolution.sequence.AFsequence;
import java.util.logging.Logger;

/**
 * @author Pablo Bousquets
 *
 * Parser that returns a list of allele frequency patterns
 *
 */
public class AlleleFractionPatternParser extends AbstractXMLObjectParser {

    public static final String ALLELEFRACTIONS = "afalignment"; //Reference to main XML class
    public static final String AFSEQUENCE = "afseq";

    public String getParserName() {
        return ALLELEFRACTIONS;
    }

    public static final String ID = "id"; //Reference to the ALLELEFRACTION ID

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        AFalignment alignment = new AFalignment();

        final DataType dataType = DataTypeUtils.getDataType(xo);

        if (dataType == null) {
            throw new XMLParseException("dataType attribute expected for alignment element");
        }

        //alignment.setDataType(dataType); //TODO

        for (int i = 0; i < xo.getChildCount(); i++) {
            final Object child = xo.getChild(i);

            if (child instanceof AFsequence) {
                alignment.addSequence((AFsequence) child); //
            } else if (child instanceof DataType) {
                // already dealt with (the exception above)
            } else {
                throw new XMLParseException("Unknown child element found in alignment");
            }
        }

        final Logger logger = Logger.getLogger("dr.evoxml");
        logger.info("Read alignment" + (xo.hasAttribute(XMLParser.ID) ? ": " + xo.getId() : "") +
                "\n  Sequences = " + alignment.getAFSequenceCount() +
                "\n      Sites = " + alignment.getAFSiteCount());

        return alignment;
    }


    /**
     * @return the taxon for this sequences.
     */

    public String getParserDescription() {
        return "This element represents an allele frequency alignment.";
    }

    public Class getReturnType() {
        return AFalignment.class;
    }

    public String getExample() {

        return //TODO: check the datatype
                "<!-- An alignment of three short DNA sequences -->\n" +
                        "<afalignment missing=\"-?\" dataType=\"AlleleFrequencies\">\n" +
                        "  <afsequence>\n" +
                        "    <taxon idref=\"taxon1\"/>\n" +
                        "    0.1,0.2,0,1\n" +
                        "  </afsequence>\n" +
                        "  <afsequence>\n" +
                        "    <taxon idref=\"taxon2\"/>\n" +
                        "    0.11,0.2,0,1\n" +
                        "  </afsequence>\n" +
                        "  <afsequence>\n" +
                        "    <taxon idref=\"taxon3\"/>\n" +
                        "    0.24,0.2,0,1\n" +
                        "  </afsequence>\n" +
                        "</afalignment>\n";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new StringAttributeRule(
                            DataType.DATA_TYPE,
                            "The data type",
                            DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class)
            ),
            new ElementRule(AFsequence.class,
                    "A string of numbers representing the AFs",
                    1, Integer.MAX_VALUE)
    };
}

