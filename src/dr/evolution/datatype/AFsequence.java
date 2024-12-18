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
    private boolean hitLimit = false;
    private String name;

    public static final AFsequence INSTANCE = new AFsequence();

    public AFsequence() {}

    public AFsequence(int nStates) { this.stateCount=nStates;}; // DM TODO: I implemented this quick fix for the test only. A proper state count it is needed for the initialization of the tipStatesModel, and probably more things down the road.
    //DM TODO: We'll probably want to add a required parameter in AFSequence's XML parser with this information, and use a similar constructor. Alternatively, it would be better to modify the setStates method of the abstract super(), but it is final, so we can't. This would complicate things if we wanted to use the BEAST app without modification.

    public AFsequence(String sequenceString){
        setSequenceString(sequenceString);
        this.taxon = getTaxon();
        this.sequenceLength = sequence.length;
    }

    public AFsequence(String sequenceString, int nStates){
        setSequenceString(sequenceString);
        this.taxon = getTaxon();
        this.sequenceLength = sequence.length;
        this.stateCount=nStates;
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
            value = checkRange(value, 0, 1);
            double_seq[i] = value;
        }
        this.sequence = FlipFlopUtils.mapDoubleRangeToInt(double_seq, 0, 1, 200); //TODO: parametize the precision?
    }

    public void setStateCount(int nStates){
        this.stateCount=nStates;
    }
    public int[] getSequence(){ return this.sequence; }

    public double checkRange(double value, int min, int max) {
        if (value < min || value > max) {
            throw new java.lang.RuntimeException(String.format("Allele frequencies are expected to be in range [%d-%d]", min, max));
        }

        if (value == 0)  {
            setHitLimit(true);
            value = 0.0000000001;
        }

        if (value == 1)  {
            setHitLimit(true);
            value = 1-0.0000000001;
        }

        return value;
    }

    /**
     * @return the taxon for this sequences.
     */
    public Taxon getTaxon() {
        return taxon;
    }

    public int getLength() {return this.sequenceLength; }


    /**
     * @return the stateCount.
     */
    public int getStateCount(){
        return stateCount;
    }


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
    public boolean getHitLimit(){return this.hitLimit;}
    public void setHitLimit(boolean limit) {this.hitLimit = limit;}

}




