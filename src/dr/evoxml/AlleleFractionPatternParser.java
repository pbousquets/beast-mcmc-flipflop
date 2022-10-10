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
import dr.evolution.alignment.Patterns;
import java.util.logging.Logger;

/**
 * @author Pablo Bousquets
 *
 * Parser that returns a list of allele frequency patterns
 *
 */
public class AlleleFractionPatternParser extends AbstractXMLObjectParser {

    public static final String ALLELEFRACTIONS = "allelefractions"; //Rerference to main XML class
    public static final String AF_SEQ = "afsequence"; //Reference to the sequences
    public static final String PRINT_DETAILS = "printDetails";
    public static final String PRINT_AF_PATTERN_CONTENT = "printAFPatContent";

    public static final String ID ="id"; //Reference to the ALLELEFRACTION set ID

    public static final int COUNT_INCREMENT = 100; // TODO: Dunno what it does

    public String getParserName() {
        return ALLELEFRACTIONS;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Taxa taxonList = (Taxa)xo.getChild(Taxa.class); // Retrieve the taxon id

        AlleleFraction allelefraction = (AlleleFraction)xo.getChild(AlleleFraction.class);

        String[] sites = ((String) xo.getElementFirstChild(AF_SEQ)).split(","); // split the sequences
        int[] pattern = new int [sites.length];
        try {
            for (int i = 0; i < sites.length; i++) {
                pattern[i] = Integer.parseInt(sites[i]);
                checkRange(pattern[i], 0, 100);
            }
        } catch (NumberFormatException nfe) {
            throw new XMLParseException("Unable to parse AF data: " + nfe.getMessage()); //TODO: change these excepctions?
        } catch (IllegalArgumentException iae) {
            throw new XMLParseException("Unable to parse AF data: " + iae.getMessage());
        }

        Patterns allelefractionPat = new Patterns(allelefraction, taxonList); // Create Pattern object
            allelefractionPat.addPattern(pattern); // Feeds the pattern with the AF sequences
            allelefractionPat.setId((String)xo.getAttribute(ID));  //Provide the Taxon ID to this object

        if(xo.getAttribute(PRINT_DETAILS,true)){
            printDetails(allelefractionPat);
        }

        if(xo.getAttribute(PRINT_AF_PATTERN_CONTENT,true)){
            printAFContent(allelefractionPat);
        }

        return allelefractionPat;

    }

    public static void printDetails(Patterns allelefractionPat){
        Logger.getLogger("dr.evoxml").info(
                "    Locus name: "+allelefractionPat.getId()+
                        "\n    Number of Taxa: "+allelefractionPat.getPattern(0).length);
    }

    public static void printAFContent(Patterns allelefractionPat){
        Logger.getLogger("dr.evoxml").info(
                "    Locus name: "+ allelefractionPat.getId());
        int[] pat = allelefractionPat.getPattern(0);
        for(int i = 0; i < pat.length; i++){
            Logger.getLogger("dr.evoxml").info("    Taxon: "+allelefractionPat.getTaxon(i)+" "+"state: "+pat[i]);
        }
        Logger.getLogger("dr.evoxml").info("\n");
    }


    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(Taxa.class),
            new ElementRule(AlleleFraction.class),
            new ElementRule(AF_SEQ,new XMLSyntaxRule[]{
                    new ElementRule(String.class,
                            "A string of numbers representing the allele frequencies in a taxa",
                            "0,10,40,100")},false),
            new StringAttributeRule(ID, "the name of the character"),
            AttributeRule.newBooleanRule(PRINT_DETAILS, true),
            AttributeRule.newBooleanRule(PRINT_AF_PATTERN_CONTENT, true)
    };

    public String getParserDescription() {
        return "This element represents an allele frequency pattern.";
    }

    public static void checkRange(int value, int min, int max) {
        if (value < min || value > max) {
            throw new java.lang.RuntimeException(String.format("Allele frequencies are expected to be in range [%d-%d]", min, max));
        }
    }

    public Class getReturnType() {
        return Patterns.class;
    }


}