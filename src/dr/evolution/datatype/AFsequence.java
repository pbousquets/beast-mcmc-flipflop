/*
 * AFsequence.java
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

package dr.evolution.datatype;
import java.util.ArrayList;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.util.FlipFlopUtils;

/**
 * @author Pablo Bousquets
 *
 * Allele frequency data type
 */
public class AFsequence extends DataType {

    public static final String DESCRIPTION = "afsequence";
    public static int UNKNOWN_STATE_LENGTH = -1;
    public int[] sequence;
    public String sequenceString;
    public Taxon taxon;

    private int sequenceLength = 0;

    private String name;

    public static final AFsequence INSTANCE = new AFsequence();

    public AFsequence() {}

    public AFsequence(int nStates) { this.stateCount=nStates;}; // DM This is for the test only. We'll need to either require a parameter that indicates this or to add a method to modify this. A proper state count it is needed for the initialization of the tipStatesModel, and probably more things down the road

    public AFsequence(String sequenceString){ // Create the sequence object
        setSequenceString(sequenceString);
        this.taxon = getTaxon();
        this.sequenceLength = sequence.length;
    }

    public void setTaxon(Taxon taxon) {
        this.taxon = taxon;
    }

    public void setSequenceString(String sequenceString) {
        this.sequenceString = sequenceString;
        String[] split_sequence = sequenceString.split(",");
        double[] double_seq = new double[split_sequence.length];
        for (int i = 0; i < split_sequence.length; i++) {
            double value = Double.parseDouble(split_sequence[i]);
            checkRange(value, 0, 1);
            double_seq[i] = value;
        }
        this.sequence = FlipFlopUtils.mapDoubleRangeToInt(double_seq, 0, 1, 200); //TODO: parametize the precision?
    }

    public int[] getSequence(){ return this.sequence; }

    public void checkRange(double value, int min, int max) {
        if (value < min || value > max) {
            throw new java.lang.RuntimeException(String.format("Allele frequencies are expected to be in range [%d-%d]", min, max));
        }
    }

    /**
     * @return the taxon for this sequences.
     */
    public Taxon getTaxon() {
        return taxon;
    }

    public int getLength() {return this.sequenceLength; }


    @Override
    public char[] getValidChars() {
        return null;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    /**
     * @return the description of the data type
     */
    public String getDescription() {
        return DESCRIPTION;
    }

    public int getType(){
        return AF_SEQ;
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    protected String id = null;

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

}




