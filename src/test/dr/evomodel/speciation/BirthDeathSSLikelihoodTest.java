/*
 * BirthDeathSSLikelihoodTest.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.evomodel.speciation;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.speciation.BirthDeathSerialSamplingModel;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * YuleModel Tester.
 *
 * @author Alexei Drummond
 * @version 1.0
 * @since <pre>08/26/2007</pre>
 */
public class BirthDeathSSLikelihoodTest extends TestCase {

    private FlexibleTree tree;
    private FlexibleTree tree2;

    public BirthDeathSSLikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        NewickImporter importer = new NewickImporter("((1:1.0,2:1.0):1.0,3:2.0);");
        tree = (FlexibleTree) importer.importTree(null);

        importer = new NewickImporter("((3 : 1.5, 4 : 0.5) : 1 , (1 : 2, 2 : 1) : 3);");
        tree2 = (FlexibleTree) importer.importTree(null);
    }

    public void testBirthDeathLikelihoodBEAST2() {
        System.out.println("RootHeight = " + tree2.getRootHeight());
        Parameter origin = new Parameter.Default("origin", 6.0);

        final double birthRate = 2.0;
        final double deathRate = 1.0;
        final double psiRate = 0.5; // rate of sampling taxa through time
        final double sampleProbability = 0.0; // the proportion of taxa sampled, default to fix to 0
        final boolean hasFinalSample = false;
        Parameter b = new Parameter.Default("b", birthRate);
        Parameter d = new Parameter.Default("d", deathRate);
        Parameter psi = new Parameter.Default("psi", psiRate);
        Parameter p = new Parameter.Default("p", sampleProbability);
        Parameter r = new Parameter.Default("r", 0.0); // sampleBecomesNonInfectiousProb

        SpeciationModel speciationModel = new BirthDeathSerialSamplingModel(b, d, psi, p, false, r, hasFinalSample, origin, Units.Type.YEARS);
        Likelihood likelihood = new SpeciationLikelihood(tree2, speciationModel, "bdss.like");

        assertEquals(-19.0198, likelihood.getLogLikelihood(), 1e-5);
    }

    final double birthRate = 1.0;
    final double deathRate = 0.5;
    final double p = 1.0;
    final double psi = 0.0;

    public void testBirthDeathLikelihoodP0() {
        assertEquals(BirthDeathSerialSamplingModel.p0(birthRate, deathRate, p, psi, 1.0), 0.28236670080320814);
    }

//    public void testP1() {
//        assertEquals(BirthDeathSerialSamplingModel.p1(birthRate, deathRate, p, psi, 1.0), 0.31236180503535266);
//    }

//    public void testPureBirthLikelihood() {
//        likelihoodTester(tree, 1, 0, -2.8219461696520542);
//    }

//    public void testBirthDeathLikelihood() {
//        likelihoodTester(tree, birthRate, deathRate, null, -4.633233508436623);
//    }

    public void testBirthDeathLikelihoodOrigin() {
        System.out.println("RootHeight = " + tree.getRootHeight());
        Parameter origin = new Parameter.Default("origin", 50);
        likelihoodTester(tree, birthRate, deathRate, origin, -29.529647743897872);
    }

    private void likelihoodTester(Tree tree, double birthRate, double deathRate, Parameter origin, double logL) {

        Parameter b = new Parameter.Default("b", birthRate);
        Parameter d = new Parameter.Default("d", deathRate);
        Parameter psi = new Parameter.Default("psi", this.psi);
        Parameter p = new Parameter.Default("p", this.p);
        Parameter r = new Parameter.Default("r", 0.5);
        Parameter fTime = new Parameter.Default("time", 0.0);

        SpeciationModel speciationModel = new BirthDeathSerialSamplingModel(b, d, psi, p, false, r, true, origin, Units.Type.YEARS);
        Likelihood likelihood = new SpeciationLikelihood(tree, speciationModel, "bdss.like");

        assertEquals(logL, likelihood.getLogLikelihood());
    }

    public static Test suite() {
        return new TestSuite(BirthDeathSSLikelihoodTest.class);
    }
}