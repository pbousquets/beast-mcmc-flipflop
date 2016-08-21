/*
 * DebugUtils.java
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

package dr.inference.mcmc;

/**
 * ${CLASS_NAME}
 *
 * @author Andrew Rambaut
 * @version $Id$
 *
 * $HeadURL$
 *
 * $LastChangedBy$
 * $LastChangedDate$
 * $LastChangedRevision$
 */

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.reflect.TypeToken;
import com.google.gson.stream.JsonWriter;
import dr.inference.model.Collectable;
import dr.inference.model.Likelihood;
import dr.math.MathUtils;

import java.io.*;
import java.lang.reflect.Type;
import java.util.*;

public class DebugUtils {
    private final static boolean FIXED_FILE_NAME = true;

//    /**
//     * Writes out the current state in a human readable format to help debugging.
//     * If it fails, then returns false but does not stop.
//     * @param file the file
//     * @param state the current state number
//     * @return success
//     */
//    public static boolean writeStateToFile(File file, long state, double lnL) {
//        OutputStream fileOut = null;
//        try {
//            fileOut = new FileOutputStream(file);
//            PrintStream out = new PrintStream(fileOut);
//
//            int[] rngState = MathUtils.getRandomState();
//            out.print("rng");
//            for (int i = 0; i < rngState.length; i++) {
//                out.print("\t");
//                out.print(rngState[i]);
//            }
//            out.println();
//
//            out.print("\nstate\t");
//            out.println(state);
//
//            out.print("lnL\t");
//            out.println(lnL);
//
//            for (Parameter parameter : Parameter.CONNECTED_SET) {
//                out.print(parameter.getParameterName());
//                out.print("\t");
//                out.print(parameter.getDimension());
//                for (int dim = 0; dim < parameter.getDimension(); dim++) {
//                    out.print("\t");
//                    out.print(parameter.getParameterValue(dim));
//                }
//                out.println();
//            }
//
//            for (Model model : Model.CONNECTED_SET) {
//                if (model instanceof TreeModel) {
//                    out.print(model.getModelName());
//                    out.print("\t");
//                    out.println(((TreeModel) model).getNewick());
//                }
//            }
//
//            out.close();
//            fileOut.close();
//        } catch (IOException ioe) {
//            System.err.println("Unable to write file: " + ioe.getMessage());
//            return false;
//        }
//
////        for (Likelihood likelihood : Likelihood.CONNECTED_SET) {
////            System.err.println(likelihood.getId() + ": " + likelihood.getLogLikelihood());
////        }
//
//        return true;
//    }
//
//    /**
//     * Attempts to read the current state from a state dump file. This should be a state
//     * dump created using the same XML file (some rudimentary checking of this is done).
//     * If it fails then it will throw a RuntimeException. If successful it will return the
//     * current state number.
//     * @param file the file
//     * @return the state number
//     */
//    public static long readStateFromFile(File file, double[] lnL) {
//        long state = -1;
//
//        try {
//            FileReader fileIn = new FileReader(file);
//            BufferedReader in = new BufferedReader(fileIn);
//
//            int[] rngState = null;
//
//            String line = in.readLine();
//            String[] fields = line.split("\t");
//            if (fields[0].equals("rng")) {
//                // if there is a random number generator state present then load it...
//                try {
//                    rngState = new int[fields.length - 1];
//                    for (int i = 0; i < rngState.length; i++) {
//                        rngState[i] = Integer.parseInt(fields[i + 1]);
//                    }
//
//                } catch (NumberFormatException nfe) {
//                    throw new RuntimeException("Unable to read state number from state file");
//                }
//
//                line = in.readLine();
//                fields = line.split("\t");
//            }
//
//            try {
//                if (!fields[0].equals("state")) {
//                    throw new RuntimeException("Unable to read state number from state file");
//                }
//                state = Long.parseLong(fields[1]);
//            } catch (NumberFormatException nfe) {
//                throw new RuntimeException("Unable to read state number from state file");
//            }
//
//            line = in.readLine();
//            fields = line.split("\t");
//            try {
//                if (!fields[0].equals("lnL")) {
//                    throw new RuntimeException("Unable to read lnL from state file");
//                }
//                if (lnL != null) {
//                    lnL[0] = Double.parseDouble(fields[1]);
//                }
//            } catch (NumberFormatException nfe) {
//                throw new RuntimeException("Unable to read lnL from state file");
//            }
//
//            for (Parameter parameter : Parameter.FULL_SET) {
//                line = in.readLine();
//                fields = line.split("\t");
////                if (!fields[0].equals(parameter.getParameterName())) {
////                    System.err.println("Unable to match state parameter: " + fields[0] + ", expecting " + parameter.getParameterName());
////                }
//                int dimension = Integer.parseInt(fields[1]);
//
//                if (dimension != parameter.getDimension()) {
//                    System.err.println("Unable to match state parameter dimension: " + dimension + ", expecting " + parameter.getDimension());
//                }
//
//                if (fields[0].equals("branchRates.categories.rootNodeNumber")) {
//                    System.out.println("eek");
//                    double value = Double.parseDouble(fields[2]);
//                    parameter.setParameterValue(0, 160.0);
//                } else {
//                    for (int dim = 0; dim < parameter.getDimension(); dim++) {
//                        parameter.setParameterValue(dim, Double.parseDouble(fields[dim + 2]));
//                    }
//                }
//
//            }
//
//            // load the tree models last as we get the node heights from the tree (not the parameters which
//            // which may not be associated with the right node
//            Set<String> expectedTreeModelNames = new HashSet<String>();
//            for (Model model : Model.FULL_SET) {
//                if (model instanceof TreeModel) {
//                    expectedTreeModelNames.add(model.getModelName());
//                }
//            }
//
//            // Read in all (possibly more than one) tree
//            while((line = in.readLine()) != null) {
//                fields = line.split("\t");
//                boolean treeFound = false;
//
//                for (Model model : Model.FULL_SET) {
//                    if (model instanceof TreeModel && fields[0].equals(model.getModelName())) {
//                        treeFound = true;
//                        NewickImporter importer = new NewickImporter(fields[1]);
//                        Tree tree = importer.importNextTree();
//                        ((TreeModel) model).beginTreeEdit();
//                        ((TreeModel) model).adoptTreeStructure(tree);
//                        ((TreeModel) model).endTreeEdit();
//
//                        expectedTreeModelNames.remove(model.getModelName());
//                    }
//                }
//
//                if (!treeFound) {
//                    throw new RuntimeException("Unable to match state parameter: " + fields[0]);
//                }
//            }
//
//            if (expectedTreeModelNames.size() > 0) {
//                StringBuilder sb = new StringBuilder();
//                for (String notFoundName : expectedTreeModelNames) {
//                    sb.append("Expecting, but unable to match state parameter:" + notFoundName + "\n");
//                }
//                throw new RuntimeException(sb.toString());
//            }
//
//            if (rngState != null) {
//                MathUtils.setRandomState(rngState);
//            }
//
//            in.close();
//            fileIn.close();
//            for (Likelihood likelihood : Likelihood.CONNECTED_SET) {
//                likelihood.makeDirty();
//            }
//        } catch (IOException ioe) {
//            throw new RuntimeException("Unable to read file: " + ioe.getMessage());
//        } catch (Importer.ImportException ie) {
//            throw new RuntimeException("Unable to import tree: " + ie.getMessage());
//        }
//
//        return state;
//    }

    /**
     * Writes out the current state in a human readable format to help debugging.
     * If it fails, then returns false but does not stop.
     * @param file the file
     * @param state the current state number
     * @return success
     */
    public static boolean writeStateToFile(File file, long state, double lnL) {

        if (FIXED_FILE_NAME) {
            file = new File("debug.json");
        }

        try {
            JsonWriter writer = new JsonWriter(new FileWriter(file));
            writer.setIndent("  ");

            Map<String, Map<String, Object>> stateMap = new LinkedHashMap<String, Map<String, Object>>();

            Map<String, Object> settingsMap = new HashMap<String, Object>();

            int[] rngState = MathUtils.getRandomState();
            List<Number> rngList = new ArrayList<Number>();
            for (int i : rngState) {
                rngList.add(i);
            }

            settingsMap.put("rng", rngList);
            settingsMap.put("state", state);
            settingsMap.put("lnL", lnL);

            stateMap.put("beast", settingsMap);

            for (Collectable collectable: Collectable.FULL_SET) {
                collectable.saveModelState(stateMap);
            }
            Gson gson = new GsonBuilder().setPrettyPrinting().create();
            Type collectionType = new TypeToken<Map<String, Map<String, Object>>>(){}.getType();
            gson.toJson(stateMap, collectionType, writer);

            writer.close();
        } catch (IOException ioe) {
            System.err.println("Unable to write file: " + ioe.getMessage());
            return false;
        }

        return true;
    }

    /**
     * Attempts to read the current state from a state dump file. This should be a state
     * dump created using the same XML file (some rudimentary checking of this is done).
     * If it fails then it will throw a RuntimeException. If successful it will return the
     * current state number.
     * @param file the file
     * @return the state number
     */
    public static long readStateFromFile(File file, double[] lnL) {

        if (FIXED_FILE_NAME) {
            file = new File("debug.json");
        }

        long state;

        try {
            FileReader fileIn = new FileReader(file);
            BufferedReader in = new BufferedReader(fileIn);

            // Deserialization
            Type collectionType = new TypeToken<Map<String, Map<String, Object>>>(){}.getType();
            Gson gson = new Gson();
            Map<String, Map<String, Object>> stateMap = gson.fromJson(in, collectionType);

            in.close();
            fileIn.close();

            Map<String, Object> settingsMap = stateMap.get("beast");

            int[] rngState = null;
            if (settingsMap.containsKey("rng")) {
                List<Integer> rngStateList = ((List<Integer>) settingsMap.get("rng"));
                rngState = new int[rngStateList.size()];
                int i = 0;
                for (Object e : rngStateList) {
                    rngState[i++] = ((Number)e).intValue();
                }
            }

            if (!settingsMap.containsKey("state")) {
                throw new RuntimeException("Unable to read state number from state file");
            }
            state = ((Number)settingsMap.get("state")).longValue();

            if (lnL != null) {
                if (!settingsMap.containsKey("lnL")) {
                    throw new RuntimeException("Unable to read lnL from state file");
                }
                lnL[0] = ((Number) settingsMap.get("lnL")).doubleValue();
            }

            if (rngState != null) {
                MathUtils.setRandomState(rngState);
            }

            for (Collectable collectable : Collectable.FULL_SET) {
                collectable.loadModelState(stateMap);
            }

            for (Likelihood likelihood : Likelihood.CONNECTED_SET) {
                likelihood.makeDirty();
            }
        } catch (IOException ioe) {
            throw new RuntimeException("Unable to read state file: " + ioe.getMessage());
        }

        return state;
    }

}
