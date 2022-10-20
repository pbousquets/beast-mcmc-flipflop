/*
 * DoublePatterns.java
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
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 * A concrete implementation of PatternList. Patterns can be added and
 * removed from the list individually or in bulk from an alignment.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: Patterns.java,v 1.10 2005/07/08 11:27:53 rambaut Exp $
 */
public class DoublePatterns implements AFPatternList {

    public static final int COUNT_INCREMENT = 100;

    /**
     * number of patterns
     */
    protected int patternCount = 0;

    /**
     * length of patterns
     */
    protected int patternLength = 0;

    /**
     * weights of each pattern
     */
    protected double[] weights = new double[COUNT_INCREMENT];

    /**
     * site patterns [pattern][taxon]
     */
    protected double[][] patterns = new double[COUNT_INCREMENT][];

    protected DataType dataType = null;

    protected TaxonList taxonList = null;

    /**
     * Constructor
     */
    public DoublePatterns() {
    }

    /**
     * Constructor
     */
    public DoublePatterns(DataType dataType) {
        this.dataType = dataType;
    }

    /**
     * Constructor
     */
    public DoublePatterns(DataType dataType, TaxonList taxonList) {
        this.dataType = dataType;
        this.taxonList = taxonList;
        patternLength = taxonList.getTaxonCount();
    }

    /**
     * Constructor
     */
    public DoublePatterns(AFPatternList patternList) {
        addPatterns(patternList);
    }

    /**
     * Constructor
     */

    private void subSetPatterns(int subSet, int subSetCount) {
        if (subSetCount > 0) {
            // if we are using subSetCount then cut it down to only the subset we want...
            int div = patternCount / subSetCount;
            int rem = patternCount % subSetCount;

            int start = 0;
            for (int i = 0; i < subSet; i++) {
                start += div + (i < rem ? 1 : 0);
            }

            int newPatternCount = div;
            if (subSet < rem) {
                newPatternCount++;
            }

            double[][] newPatterns = new double[newPatternCount][];
            double[] newWeights = new double[newPatternCount];
            for (int i = 0; i < newPatternCount; i++) {
                newPatterns[i] = patterns[start + i];
                newWeights[i] = weights[start + i];
            }
            patterns = newPatterns;
            weights = newWeights;

            patternCount = newPatternCount;
        }
    }

    /**
     * adds patterns to the list from a SiteList
     */
    public void addPatterns(AFPatternList patternList) {

        if (patternList == null) {
            return;
        }

        if (taxonList == null) {
            taxonList = patternList;
            patternLength = taxonList.getTaxonCount();
        }

        if (dataType == null) {
            dataType = patternList.getDataType();
        } else if (dataType != patternList.getDataType()) {
            throw new IllegalArgumentException("Patterns' existing DataType does not match that of added PatternList" + dataType + " vs " + patternList.getDataType());
        }

        for (int i = 0; i < patternList.getPatternCount(); i++) {
            double[] pattern = patternList.getPattern(i);

            // don't add patterns that are all gaps or all ambiguous
                addPattern(pattern, patternList.getPatternWeight(i));

        }

    }

    /**
     * adds a pattern to the pattern list with a default weight of 1
     */
    public void addPattern(double[] pattern) {
        addPattern(pattern, 1.0);
    }

    /**
     * adds a pattern to the pattern list
     */
    public void addPattern(double[] pattern, double weight) {

        if (patternLength == 0) {
            patternLength = pattern.length;
        }

        if (patternLength != 0 && pattern.length != patternLength) {
            throw new IllegalArgumentException("Added pattern's length (" + pattern.length + ") does not match those of existing patterns (" + patternLength + ")");
        }

        for (int i = 0; i < patternCount; i++) {

            if (comparePatterns(patterns[i], pattern)) {

                weights[i] += weight;
                return;
            }
        }

        if (patternCount == patterns.length) {
            double[][] newPatterns = new double [patternCount + COUNT_INCREMENT][];
            double[] newWeights = new double[patternCount + COUNT_INCREMENT];
            for (int i = 0; i < patternCount; i++) {
                newPatterns[i] = patterns[i];
                newWeights[i] = weights[i];
            }
            patterns = newPatterns;
            weights = newWeights;
        }

        patterns[patternCount] = pattern;
        weights[patternCount] = weight;
        patternCount++;
    }

    /**
     * removes a pattern from the pattern list
     */
    public void removePattern(double[] pattern) {

        int index = -1;
        for (int i = 0; i < patternCount; i++) {

            if (comparePatterns(patterns[i], pattern)) {
                index = i;
                break;
            }
        }

        if (index == -1) throw new IllegalArgumentException("Pattern not found");

        weights[index] -= 1;
        if (weights[index] == 0 && patternCount > 1) {
            patterns[index] = patterns[patternCount - 1];
            patterns[patternCount - 1] = null;
            weights[index] = weights[patternCount - 1];
            patternCount--;
        }

    }

    /**
     * removes all patterns from the pattern list
     */
    public void removeAllPatterns() {
        patternCount = 0;
        for (int i = 0; i < patterns.length; i++) patterns[i] = null;
    }


    protected boolean comparePatterns(double[] pattern1, double[] pattern2) {

        int len = pattern1.length;
        for (int i = 0; i < len; i++) {
            if (pattern1[i] != pattern2[i]) {
                return false;
            }
        }

        return true;
    }

    // **************************************************************
    // PatternList IMPLEMENTATION
    // **************************************************************

    /**
     * @return number of patterns
     */
    public int getPatternCount() {
        return patternCount;
    }

    /**
     * @return number of states for this siteList
     */
    public int getStateCount() {
        return dataType.getStateCount();
    }

    /**
     * Gets the length of the pattern strings which will usually be the
     * same as the number of taxa
     *
     * @return the length of patterns
     */
    public int getPatternLength() {
        return patternLength;
    }

    /**
     * Gets the pattern as an array of state numbers (one per sequences)
     *
     * @return the pattern at patternIndex
     */
    public double[] getPattern(int patternIndex) {
        return patterns[patternIndex];
    }

    /**
     * @return state at (taxonIndex, patternIndex)
     */
    public double getPatternState(int taxonIndex, int patternIndex) {
        return patterns[patternIndex][taxonIndex];
    }

    /**
     * Gets the weight of a site pattern
     */
    public double getPatternWeight(int patternIndex) {
        return weights[patternIndex];
    }

    /**
     * @return the array of pattern weights
     */
    public double[] getPatternWeights() {
        double[] w = new double[weights.length];
        for (int i = 0; i < weights.length; i++) w[i] = weights[i];
        return w;
    }

    /**
     * @return the DataType of this siteList
     */
    public DataType getDataType() {
        return dataType;
    }

    /**
     * @return the frequency of each state
     */

    /*
    public double[] getStateFrequencies() {
        return AFPatternList.Utils.empiricalStateFrequencies(this);
    }*/

    @Override
    public boolean areUnique() {
        return true;
    }

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

    /**
     * @return a count of the number of taxa in the list.
     */
    public int getTaxonCount() {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxonList.getTaxonCount();
    }

    /**
     * @return the ith taxon.
     */
    public Taxon getTaxon(int taxonIndex) {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxonList.getTaxon(taxonIndex);
    }

    public void setTaxon(TaxonList taxonList) {
        this.taxonList = taxonList;
    }

    /**
     * @return the ID of the ith taxon.
     */
    public String getTaxonId(int taxonIndex) {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxonList.getTaxonId(taxonIndex);
    }

    /**
     * returns the index of the taxon with the given id.
     */
    public int getTaxonIndex(String id) {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxonList.getTaxonIndex(id);
    }

    /**
     * returns the index of the given taxon.
     */
    public int getTaxonIndex(Taxon taxon) {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxonList.getTaxonIndex(taxon);
    }

    public List<Taxon> asList() {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxonList.asList();
    }

    public Iterator<Taxon> iterator() {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxonList.iterator();
    }


    /**
     * @param taxonIndex the index of the taxon whose attribute is being fetched.
     * @param name       the name of the attribute of interest.
     * @return an object representing the named attributed for the given taxon.
     */
    public Object getTaxonAttribute(int taxonIndex, String name) {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        return taxonList.getTaxonAttribute(taxonIndex, name);
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

    // ========= Mask =========
    // indexes to mask sth., e.g. taxon index whose state is unknown character in microsatellite
    protected Set<Integer> maskSet = new HashSet<Integer>();

    // no duplication, if duplicate, not add
    public boolean addMask(int index) {
        return maskSet.add(index);
    }

    public boolean isMasked(int index) {
        return maskSet.contains(index);
    }

    public boolean hasMask() {
        return maskSet.size() > 0;
    }

    public void clearMask() {
        maskSet.clear();
    }

    public Set<Integer> getMaskSet() {
        return maskSet;
    }

    /**
     * @return the ith taxon not masked.
     */
    public Taxon getTaxonMasked(int taxonIndex) {
        if (taxonList == null) throw new RuntimeException("Patterns has no TaxonList");
        if (isMasked(taxonIndex)) {
            return null;
        }
        return taxonList.getTaxon(taxonIndex);
    }

}
