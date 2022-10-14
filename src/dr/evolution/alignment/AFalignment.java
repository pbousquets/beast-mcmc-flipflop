/*
 * SimpleAlignment.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.DataType;
import dr.evolution.sequence.SequenceList;

import dr.evolution.sequence.AFsequence;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import java.util.ArrayList;


/**
 * A simple alignment class that implements Allele Frequencies instead of nucleotides.
 *
 * @author Pablo Bousquets
 */
@SuppressWarnings("serial")

public class AFalignment extends AFsequence {

    public ArrayList<AFsequence> alignment = new ArrayList<AFsequence>();

    public int ntaxa = 0;
    public int alignment_length = 0;

    public void addSequence(AFsequence sequence){
        alignment.add(sequence);

        if (ntaxa == 0)  {
            alignment_length = sequence.getLength();
        }

        if (alignment_length != sequence.getLength()) {
            throw new java.lang.RuntimeException(String.format("Different sites found among the taxa. Please, check that all have the same amount of sites"));
        }

        ntaxa++;
    }

    public double[] getTaxonSequence(int taxonIndex){
        AFsequence s = this.alignment.get(taxonIndex);
        return s.getSequence(); //returns the actual AF values
    }

    public int getAFSequenceCount(){
        return this.ntaxa; //returns the number of sequences
    }

    public int getAFSiteCount(){
        return this.alignment_length; //returns the number of sites
    }

}// END: class
