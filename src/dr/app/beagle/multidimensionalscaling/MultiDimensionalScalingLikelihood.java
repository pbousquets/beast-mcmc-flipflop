/*
 * MultiDimensionalScalingLikelihood.java
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

package dr.app.beagle.multidimensionalscaling;

import dr.evomodel.antigenic.MultidimensionalScalingLikelihood;
import dr.inference.model.*;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class MultiDimensionalScalingLikelihood extends AbstractModelLikelihood implements Reportable {

    public static final String REQUIRED_FLAGS_PROPERTY = "mds.required.flags";

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();
        sb.append(getId() + ": " + getLogLikelihood());
        return sb.toString();
    }

    public enum ObservationType {
        POINT,
        UPPER_BOUND,
        LOWER_BOUND,
        MISSING
    }

    public final static String MULTIDIMENSIONAL_SCALING_LIKELIHOOD = "multiDimensionalScalingLikelihood";

//    public MultiDimensionalScalingLikelihood(
//            int mdsDimension,
//            Parameter mdsPrecision,
//            MatrixParameter locationsParameter,
//            DataTable<double[]> dataTable,
//            boolean reorderData) {
//        this(mdsDimension, mdsPrecision, locationsParameter, dataTable, false, reorderData);
//    }

    /**
     * A simple constructor for a fully specified symmetrical data matrix
     * @param mdsDimension
     * @param mdsPrecision
     * @param locationsParameter
     * @param dataTable
     * @param isLeftTruncated
     * @param reorderData
     */
    public MultiDimensionalScalingLikelihood(
            int mdsDimension,
            Parameter mdsPrecision,
            MatrixParameterInterface locationsParameter,
            DataTable<double[]> dataTable,
            boolean isLeftTruncated,
            boolean reorderData) {

        super(MULTIDIMENSIONAL_SCALING_LIKELIHOOD);

        this.mdsDimension = mdsDimension;
        this.isLeftTruncated = isLeftTruncated;

        // construct a compact data table
        String[] rowLabelsOriginal = dataTable.getRowLabels();
//        String[] columnLabels = dataTable.getRowLabels();

        int rowCount = dataTable.getRowCount();
        locationCount = rowCount;

        int[] permute = null;
        if (reorderData) {
            permute = getPermutation(rowLabelsOriginal, locationsParameter);
        } else {
            permute = new int[locationCount];
            for (int i = 0; i < locationCount; ++i) {
                permute[i] = i; // identity
            }
        }

        String[] rowLabels = new String[locationCount];

        int observationCount = rowCount * rowCount;
//        double[] observations = new double[observationCount];
        observations = new double[observationCount];
        ObservationType[] observationTypes = new ObservationType[observationCount];

        double[][] tmp = new double[rowCount][rowCount];

        for (int i = 0; i < rowCount; i++) {
            rowLabels[i] = rowLabelsOriginal[permute[i]];

            double[] dataRow = dataTable.getRow(permute[i]);

            for (int j = i + 1; j < rowCount; j++) {
                tmp[i][j] = tmp[j][i] = dataRow[permute[j]];
            }
        }

        int u = 0;
        for (int i = 0; i < rowCount; i++) {
            for (int j = 0; j < rowCount; j++) {
                observations[u] = (i == j ? 0 : tmp[i][j]);
                observationTypes[u] = ObservationType.POINT;
                u++;
            }
        }

        initialize(mdsDimension, mdsPrecision, isLeftTruncated, locationsParameter,
                rowLabels, observations, observationTypes);
    }

//    private class Data {
//        int observationCount;
//        double[] observations;
//        ObservationType[] observationTypes;
//
//        Data(int observationCount, double[] observations, ObservationType[] observationTypes) {
//            this.observationCount = observationCount;
//            this.observations = observations;
//            this.observationTypes = observationTypes;
//        }
//    }

    public double[] getObservations() { return observations; }    // TODO Grab from core when needed to save space

    public MatrixParameterInterface getMatrixParameter() { return locationsParameter; }

    private int[] getPermutation(String[] source, MatrixParameterInterface destination) {

        if (source.length != destination.getColumnDimension()) {
            throw new IllegalArgumentException("Dimension mismatch");
        }

        final int length = source.length;

        Map<String,Integer> map = new HashMap<String, Integer>(destination.getColumnDimension());
        for (int i = 0; i < length; ++i) {
            map.put(source[i],i);
        }

        int[] permute = new int[length];
        for (int i = 0; i < length; ++i) {
            Integer p = map.get(destination.getParameter(i).getParameterName());
            if (p == null) {
                throw new IllegalArgumentException("Missing label");
            }
            permute[i] = p;
        }

        return permute;
    }

    private MultiDimensionalScalingCore getCore() {
        long computeMode = 0;
        String r = System.getProperty(REQUIRED_FLAGS_PROPERTY);
        if (r != null) {
            computeMode = Long.parseLong(r.trim());
        }

        MultiDimensionalScalingCore core;
        if (computeMode > 0) {
            System.err.println("Attempting to use a native MDS core with flag: " + computeMode + "; may the force be with you ....");
            core = new MassivelyParallelMDSImpl();
            flags = computeMode;
        } else {
            core = new MultiDimensionalScalingCoreImpl2();
        }
//        switch (computeMode) {
//            case 1:
//                System.err.println("Attempting to use a native MDS core; may the force be with you ....");
//                core = new MassivelyParallelMDSImpl();
//                break;
//            default:
//                core = new MultiDimensionalScalingCoreImpl2();
//        }
        return core;
    }

    public int getMdsDimension() { return mdsDimension; }

    public int getLocationCount() { return locationCount; }

    protected void initialize(
            final int mdsDimension,
            final Parameter mdsPrecision,
            final boolean isLeftTruncated,
            final MatrixParameterInterface locationsParameter,
            final String[] locationLabels,
            final double[] observations,
            final ObservationType[] observationTypes) {

        this.mdsCore = getCore();

        if (isLeftTruncated) {
            flags |= MultiDimensionalScalingCore.LEFT_TRUNCATION;
        }

        System.err.println("Initializing with flags: " + flags);

        this.mdsCore.initialize(mdsDimension, locationCount, flags);
        this.locationLabels = locationLabels;

        this.locationsParameter = locationsParameter;
        setupLocationsParameter(this.locationsParameter);
        addVariable(locationsParameter);

        this.mdsPrecisionParameter = mdsPrecision;
        addVariable(mdsPrecision);

        mdsCore.setParameters(mdsPrecisionParameter.getParameterValues());
        mdsCore.setPairwiseData(observations);
//        for (int i = 0; i < locationCount; i++) {
//            mdsCore.updateLocation(i, locationsParameter.getColumnValues(i));
//        }
        mdsCore.updateLocation(-1, locationsParameter.getParameterValues());

        // make sure everything is calculated on first evaluation
        makeDirty();
    }

    protected void setupLocationsParameter(MatrixParameterInterface locationsParameter) {
        final boolean exisitingParameter = locationsParameter.getColumnDimension() > 0;

        if (exisitingParameter){
            if (locationsParameter.getColumnDimension() != locationCount){
                throw new RuntimeException("locationsParameter column dimension ("+locationsParameter.getColumnDimension()+") is not equal to the locationCount ("+locationCount+")");
            }
            if (locationsParameter.getRowDimension() != mdsDimension){
                throw new RuntimeException("locationsParameter row dimension ("+locationsParameter.getRowDimension()+") is not equal to the mdsDimension ("+mdsDimension+")");
            }
        } else {
//            locationsParameter.setColumnDimension(mdsDimension);
//            locationsParameter.setRowDimension(locationCount);
            throw new IllegalArgumentException("Dimensions on matrix must be set");
        }

        for (int i = 0; i < locationLabels.length; i++) {
            if (exisitingParameter) {
                if (locationsParameter.getParameter(i).getParameterName().compareTo(locationLabels[i]) != 0) {
                    throw new RuntimeException("Mismatched trait parameter name (" + locationsParameter.getParameter(i).getParameterName() +
                            ") and data dimension name (" + locationLabels[i] + ")");
                }
            } else {
                locationsParameter.getParameter(i).setId(locationLabels[i]);
            }
        }

        for (int i = 0; i < locationsParameter.getColumnDimension(); ++i) {
            Parameter param = locationsParameter.getParameter(i);
            try {
                if (param.getBounds() != null) {
                    // Do nothing
                }
            } catch (NullPointerException exception) {
                param.addBounds(new Parameter.DefaultBounds(
                        Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, param.getDimension()));
            }
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) { }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // TODO Flag which cachedDistances or mdsPrecision need to be updated

        if (variable == locationsParameter) {

            if (index == -1) {

                mdsCore.updateLocation(-1, locationsParameter.getParameterValues());
            } else {

                int locationIndex = index / mdsDimension;
                mdsCore.updateLocation(locationIndex, locationsParameter.getColumnValues(locationIndex));
            }
        } else if (variable == mdsPrecisionParameter) {
            mdsCore.setParameters(mdsPrecisionParameter.getParameterValues());
        } else {
            // could be a derived class's parameter
//            throw new IllegalArgumentException("Unknown parameter");
        }
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        mdsCore.storeState();
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
        mdsCore.restoreState();
    }

    @Override
    protected void acceptState() {
        mdsCore.acceptState();
        // do nothing
    }

    public void makeDirty() {
        likelihoodKnown = false;
        mdsCore.makeDirty();
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = mdsCore.calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
        public final static String FILE_NAME = "fileName";

        public final static String TIP_TRAIT = "tipTrait";
        public final static String LOCATIONS = "locations";
        public static final String MDS_DIMENSION = "mdsDimension";
        public static final String MDS_PRECISION = "mdsPrecision";
        public static final String INCLUDE_TRUNCATION = "includeTruncation";
        public static final String USE_OLD = "useOld";
        public static final String FORCE_REORDER = "forceReorder";

        public String getParserName() {
            return MULTIDIMENSIONAL_SCALING_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String fileName = xo.getStringAttribute(FILE_NAME);
            DataTable<double[]> distanceTable;
            try {
                distanceTable = DataTable.Double.parse(new FileReader(fileName));
            } catch (IOException e) {
                throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
            }

            if (distanceTable.getRowCount() != distanceTable.getColumnCount()) {
                throw new XMLParseException("Data table is not symmetrical.");
            }

            int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);

            MatrixParameterInterface locationsParameter = (MatrixParameterInterface) xo.getElementFirstChild(LOCATIONS);

            Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

            boolean useOld = xo.getAttribute(USE_OLD, false);

            boolean includeTrauncation = xo.getAttribute(INCLUDE_TRUNCATION, false);

            boolean forceReorder = xo.getAttribute(FORCE_REORDER, false);

            if (useOld) {
                System.err.println("USE OLD");
                return new MultidimensionalScalingLikelihood(mdsDimension, includeTrauncation, mdsPrecision, (MatrixParameter)locationsParameter, distanceTable);
            } else {
                return new MultiDimensionalScalingLikelihood(mdsDimension, mdsPrecision, locationsParameter,
                        distanceTable, includeTrauncation, forceReorder);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of pairwise distance given vectors of coordinates" +
                    "for points according to the multidimensional scaling scheme of XXX & Rafferty (xxxx).";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
                AttributeRule.newIntegerRule(MDS_DIMENSION, false, "The dimension of the space for MDS"),
                new ElementRule(LOCATIONS, MatrixParameterInterface.class),
                AttributeRule.newBooleanRule(USE_OLD, true),
                AttributeRule.newBooleanRule(INCLUDE_TRUNCATION, true),
                AttributeRule.newBooleanRule(FORCE_REORDER, true),
                new ElementRule(MDS_PRECISION, Parameter.class)
        };

        public Class getReturnType() {
            return MultiDimensionalScalingLikelihood.class;
        }
    };

    private final int mdsDimension;
    private final int locationCount;
    private final boolean isLeftTruncated;

    private MultiDimensionalScalingCore mdsCore;

    private String[] locationLabels;

    private Parameter mdsPrecisionParameter;
    private MatrixParameterInterface locationsParameter;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private long flags = 0;

    private double[] observations;
}
