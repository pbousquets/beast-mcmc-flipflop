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

package dr.evolution.sequence;

import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;
import dr.util.Attributable;
import dr.util.Identifiable;

import java.util.Iterator;

/**
 * A simple class to use Allele frequencies sequences rather than nucleotides.
 *
 * @author Pablo Bousquets
 */
@SuppressWarnings("serial")
public class AFsequence implements Identifiable, Attributable {

    public double[] sequence;
    public String sequenceString;
    public Taxon taxon;

    private int sequenceLength = 0;

    public AFsequence(){
    }

    public AFsequence(String sequenceString){ // Create the sequence object
        this.sequenceString = sequenceString;
        String[] split_sequence = sequenceString.split(",");

        for (int i = 0; i < split_sequence.length; i++) {
            double value = Double.parseDouble(split_sequence[i]);
            checkRange(value,0,1);
            this.sequence[i] = value;
        }

        this.sequenceLength = split_sequence.length;
        this.taxon = getTaxon();
    }

    public int getLength() {return this.sequenceLength; }

    public void setTaxon(Taxon taxon) {
        this.taxon = taxon;
    }

    public double[] getSequence(){ return this.sequence; }

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

    /**
     * Set the DataType of the sequences.
     */
    public void setDataType(DataType dataType) {
        this.dataType = dataType;
    }

    /**
     * Set the DataType of the sequences.
     */
    public DataType guessDataType() {
        return DataType.guessDataType(sequenceString.toString());
    }


    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

    private Attributable.AttributeHelper attributes = null;

    /**
     * Sets an named attribute for this object.
     *
     * @param name  the name of the attribute.
     * @param value the new value of the attribute.
     */
    public void setAttribute(String name, Object value) {
        if (attributes == null)
            attributes = new Attributable.AttributeHelper();
        attributes.setAttribute(name, value);
    }

    /**
     * @param name the name of the attribute of interest.
     * @return an object representing the named attributed for this object.
     */
    public Object getAttribute(String name) {
        if (attributes == null)
            return null;
        else
            return attributes.getAttribute(name);
    }

    /**
     * @return an iterator of the attributes that this object has.
     */
    public Iterator<String> getAttributeNames() {
        if (attributes == null)
            return null;
        else
            return attributes.getAttributeNames();
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    protected String id = null;

    /**
     * @return the id.
     */
    public String getId() {
        return id;
    }

    /**
     * Sets the id.
     */
    public void setId(String id) {
        this.id = id;
    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    protected DataType dataType = null;

}// END: class
