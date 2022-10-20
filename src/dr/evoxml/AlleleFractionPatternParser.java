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
import dr.evolution.datatype.AFsequence;
import dr.evolution.datatype.DataType;
import dr.evoxml.util.DataTypeUtils;
import dr.evolution.alignment.DoublePatterns;
import dr.evolution.alignment.AFPatternList;
import java.util.logging.Logger;
import java.util.ArrayList;
import java.util.Arrays;

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

    public static final int COUNT_INCREMENT = 100;

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DoublePatterns patterns = new DoublePatterns();
        Taxa taxonList = new Taxa();
        ArrayList<AFsequence> seqlist = new ArrayList<AFsequence>();

        for (int i = 0; i < xo.getChildCount(); i++) {
            final Object child = xo.getChild(i);

            if (child instanceof AFsequence) {
                taxonList.addTaxon(((AFsequence) child).getTaxon());
                patterns.addPattern(((AFsequence) child).getSequence());
            } else {
                throw new XMLParseException("Unknown child element found in alignment");
            }
        }

        patterns.setTaxon(taxonList);

        final Logger logger = Logger.getLogger("dr.evoxml");
        logger.info("Site patterns '" + xo.getId() + "' created:");
        logger.info("  - Taxa count = " + patterns.getTaxonCount());
        logger.info("  - Site count = " + patterns.getPatternLength());

        return patterns;
    }

    public String getParserDescription() {
        return "This element represents an allele frequency alignment.";
    }

    public Class getReturnType() {
        return AFPatternList.class;
    }

    public String getExample() {

        return
                "<!-- An alignment of three short DNA sequences -->\n" +
                        "<afalignment id = \"alignment\">" +
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
            new ElementRule(AFsequence.class,
                    "A string of numbers representing the AFs",
                    1, Integer.MAX_VALUE)
    };
}

