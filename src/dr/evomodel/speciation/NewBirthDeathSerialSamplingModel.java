/*
 * NewBirthDeathSerialSamplingModel.java
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

package dr.evomodel.speciation;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;
import java.util.Set;

/**
 * A phylogenetic birth-death-sampling model which includes serial sampling, sampling at present, and the possibility of treatmentProbability.
 */
public class NewBirthDeathSerialSamplingModel extends MaskableSpeciationModel implements Citable {

    // extant sampling proportion
    Parameter samplingFractionAtPresent;

    // birth rate
    Parameter birthRate;

    // death rate
    Parameter deathRate;

    // serial sampling rate
    Parameter serialSamplingRate;

    // "treatmentProbability" parameter aka r aka Pr(death | lineage is sampled)
    Parameter treatmentProbability;

    // the originTime of the infection, origin > tree.getRoot();
    Parameter originTime;

    private boolean conditionOnSurvival;

    // useful constants we don't want to compute nTaxa times
    private double storedC1 = Double.NEGATIVE_INFINITY;
    private double storedC2 = Double.NEGATIVE_INFINITY;

    public NewBirthDeathSerialSamplingModel(
            Parameter birthRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter treatmentProbability,
            Parameter samplingFractionAtPresent,
            Parameter originTime,
            boolean condition,
            Type units) {

        this("NewBirthDeathSerialSamplingModel", birthRate, deathRate, serialSamplingRate, treatmentProbability, samplingFractionAtPresent, originTime, condition, units);
    }

    public SpeciationModelGradientProvider getProvider() { // This is less INTRUSIVE to the exisiting file
        return new NewBirthDeathSerialSamplingModelGradient(this);
    }

    public NewBirthDeathSerialSamplingModel(
            String modelName,
            Parameter birthRate,
            Parameter deathRate,
            Parameter serialSamplingRate,
            Parameter treatmentProbability,
            Parameter samplingFractionAtPresent,
            Parameter originTime,
            boolean condition,
            Type units) {

        super(modelName, units);

        this.birthRate = birthRate;
        addVariable(birthRate);
        birthRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.deathRate = deathRate;
        addVariable(deathRate);
        deathRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.serialSamplingRate = serialSamplingRate;
        addVariable(serialSamplingRate);
        serialSamplingRate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.samplingFractionAtPresent = samplingFractionAtPresent;
        addVariable(samplingFractionAtPresent);
        samplingFractionAtPresent.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        this.treatmentProbability = treatmentProbability;
        addVariable(treatmentProbability);
        treatmentProbability.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        this.conditionOnSurvival = condition;

        this.originTime = originTime;
        if (originTime != null) {
            addVariable(originTime);
            originTime.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }
    }

    /**
     * @param lambda   birth rate
     * @param mu   death rate
     * @param psi   proportion sampled at final time point
     * @param rho rate of sampling per lineage per unit time
     * @param t   time
     * @return the probability of no sampled descendants after time, t
     */
    public static double p0(double lambda, double mu, double psi, double rho, double c1, double c2, double t) {

        double expc1trc2 = Math.exp(-c1 * t) * (1.0 - c2);

        // Stadler 2011 p 349
        return (lambda + mu + psi + c1 * ((expc1trc2 - (1.0 + c2)) / (expc1trc2 + (1.0 + c2)))) / (2.0 * lambda);
    }

    /**
     * @param t   time
     * @return the probability of no sampled descendants after time, t
     */
    public static double logq(double c1, double c2, double t) {
        // TODO off by a factor of 4 from Gavryushkina et al (2014) and Magee and Hoehna (2021)
        // Should only be a problem if tree is being inferred and sampled ancestors are allowed
        double res = c1 * t + 2.0 * Math.log( Math.exp(-c1 * t) * (1.0 - c2) + (1.0 + c2) ); // operate directly in logspace, c1 * t too big
        return res;
    }

    private static double c1(double lambda, double mu, double psi) {
        return Math.abs(Math.sqrt(Math.pow(lambda - mu - psi, 2.0) + 4.0 * lambda * psi));
    }

    private static double c2(double lambda, double mu, double psi, double rho) {
        return -(lambda - 2.0 * rho * lambda - mu - psi)/c1(lambda, mu, psi);
    }

    private void precomputeConstants() {
        this.storedC1 = c1(lambda(), mu(), psi());
        this.storedC2 = c2(lambda(), mu(), psi(), rho());
    }

    public double p0(double t) {
        return p0(lambda(), mu(), psi(), rho(), storedC1, storedC2, t);
    }

    public double logq(double t) {
        return logq(storedC1, storedC2, t);
    }

    public double lambda() {
        if (mask != null) return mask.lambda();
        else {
            return birthRate.getValue(0);
        }
    }

    public double mu() {
        if (mask != null) return mask.mu();
        else {
            return deathRate.getValue(0);
        }
    }

    public double psi() {
        if (mask != null) return mask.psi();
        else {
            return serialSamplingRate.getValue(0);
        }
    }

    public double r() {
        if (mask != null) return mask.r();
        return treatmentProbability.getValue(0);
    }

    public double rho() {
        if (mask != null) return mask.rho();
        return samplingFractionAtPresent.getValue(0);
    }


    /**
     * Generic likelihood calculation
     *
     * @param tree the tree to calculate likelihood of
     * @return log-likelihood of density
     */
    public final double calculateTreeLogLikelihood(Tree tree) {
        precomputeConstants();

        double logL = calculateUnconditionedLogLikelihoodOverIntervals(tree);
//        double logL = calculateUnconditionedTreeLogLikelihood(tree);

        double origin = originTime.getValue(0);
        if (origin < tree.getNodeHeight(tree.getRoot())) {
            return Double.NEGATIVE_INFINITY;
        }

        if ( conditionOnSurvival ) {
            logL -= Math.log(1.0 - p0(origin));
        }

        return logL;
    }

    // Log-likelihood of tree without conditioning on anything
    public final double calculateUnconditionedTreeLogLikelihood(Tree tree) {

        double lambda = lambda();
        double mu = mu();
        double psi = psi();
        double r = r();
        double rho = rho();

        double timeZeroTolerance = Double.MIN_VALUE;
        boolean noSamplingAtPresent = rho < Double.MIN_VALUE;

        double origin = originTime.getValue(0);
        if (origin < tree.getNodeHeight(tree.getRoot())) {
            return Double.NEGATIVE_INFINITY;
        }

        // extant leaves
        int n = 0;
        // extinct leaves
        int m = 0;
        // sampled ancestors
        int k = 0;

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            if (noSamplingAtPresent || tree.getNodeHeight(node) > timeZeroTolerance) {
                m += 1;
            } else {
                n += 1;
            }
        }

        if ((!noSamplingAtPresent) && n < 1) {
            throw new RuntimeException(
                    "Sampling fraction at time zero (rho) is >0 but there are no samples at time zero"
            );
        }

        double logL = 0.0;

        logL += (double)(n + m - 1) * Math.log(lambda);

        if ( k > 0 ) {
            logL += (double)(k) * Math.log(psi * (1.0 - r));
        }

        if (!noSamplingAtPresent) {
            logL += n * Math.log(rho);
        }

        logL -= logq(origin);

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            double x = tree.getNodeHeight(tree.getInternalNode(i));
            logL -= logq(x);
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double y = tree.getNodeHeight(tree.getExternalNode(i));
            if (noSamplingAtPresent || y > timeZeroTolerance) {
                logL += Math.log(psi * (r + (1.0 - r) * p0(y))) + logq(y);
//                System.err.println("logq(y) = " + logq(y));
            }
        }

        return logL;
    }

    @Override
    public double[] getBreakPoints() {
        return new double[] { Double.POSITIVE_INFINITY };
    }

    @Override
    public double processInterval(int model, double tYoung, double tOld, int nLineages) {
        return nLineages * (Math.log(tYoung) - Math.log(tOld));
    }

    @Override
    public double processCoalescence(int model) {
        return Math.log(lambda()); // TODO Notice the natural parameterization is `log lambda`
    }

    @Override
    public double processSampling(int model, double tOld) {

        double logPsi = Math.log(psi()); // TODO Notice the natural parameterization is `log psi`
        double r = r();
        double logRho = Math.log(rho()); // TODO Notice the natural parameterization is `log rho`

        double timeZeroTolerance = Double.MIN_VALUE;
        boolean noSamplingAtPresent = rho() < Double.MIN_VALUE;

        if (noSamplingAtPresent || tOld > timeZeroTolerance) {
            return logPsi + Math.log(r + (1.0 - r) * p0(tOld));
        } else {
            return logRho;
        }
    }

    // Log-likelihood of tree without conditioning on anything
    private double calculateUnconditionedLogLikelihoodOverIntervals(Tree tree) {

        if (!(tree instanceof TreeModel)) {
            throw new IllegalArgumentException("Failed test");
        }

        double logLambda = Math.log(lambda());
        double mu = mu();
        double logPsi = Math.log(psi());
        double r = r();
        double logRho = Math.log(rho());

        double timeZeroTolerance = Double.MIN_VALUE;
        boolean noSamplingAtPresent = rho() < Double.MIN_VALUE;

        double origin = originTime.getValue(0);

        // TODO Make cached class-object
        BigFastTreeIntervals treeIntervals = new BigFastTreeIntervals((TreeModel)tree);

        double logL = 0.0;

        // TODO might be able to avoid O(nnodes) computations of q
        // nominally in every loop, we can set tYoung <- tOld and logqYoung <- logqOld
        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {
            double tYoung = treeIntervals.getIntervalTime(i);
            double tOld = tYoung + treeIntervals.getInterval(i);

            double logqYoung = logq(tYoung);
            double logqOld = logq(tOld);

            // Prob of 1 lineage surviving from interval b/n tYoung and tOld is q(tYoung)/q(tOld)
            int nLineages = treeIntervals.getLineageCount(i);
            logL += nLineages * (logqYoung - logqOld);

            // Interval ends with a coalescent or sampling event at time tOld
            if (treeIntervals.getIntervalType(i) == IntervalType.SAMPLE) {
                if (noSamplingAtPresent || tOld > timeZeroTolerance) {
                    logL += logPsi + Math.log(r + (1.0 - r) * p0(tOld));
                } else {
                    logL += logRho;
                }
            } else if (treeIntervals.getIntervalType(i) == IntervalType.COALESCENT) {
                logL += logLambda;
            } else {
                throw new RuntimeException("Birth-death tree includes non birth/death/sampling event.");
            }
        }

        // We've missed the first sample and need to add it back
        double t0 = treeIntervals.getStartTime();
        if (noSamplingAtPresent || t0 > timeZeroTolerance) {
            logL += logPsi + Math.log(r + (1.0 - r) * p0(t0));
        } else {
            logL += logRho;
        }

        // origin branch is a fake branch that doesn't exist in the tree, now compute its contribution
        double tYoung = treeIntervals.getTotalDuration();
        double tOld = origin;

        double logqYoung = logq(tYoung);
        double logqOld = logq(tOld);

        logL += (logqYoung - logqOld);

//        System.err.println(">>>>>> Intervals <<<<<<");
//        System.err.println(treeIntervals.getIntervalCount() + " " + treeIntervals.getSampleCount() + " " + treeIntervals.getTotalDuration());
//        for (int i = 0; i < treeIntervals.getIntervalCount(); ++i) {
//            System.err.println(treeIntervals.getInterval(i) + " " + treeIntervals.getLineageCount(i) + " " +
//                    treeIntervals.getIntervalTime(i) + " " + treeIntervals.getIntervalType(i));
//        }
//        System.err.println("<<<<<< Intervals >>>>>>");
        return logL;
    }



    public double calculateTreeLogLikelihood(Tree tree, Set<Taxon> exclude) {
        if (exclude.size() == 0) return calculateTreeLogLikelihood(tree);
        throw new RuntimeException("Not implemented!");
    }

    public void mask(SpeciationModel mask) {
        if (mask instanceof NewBirthDeathSerialSamplingModel) {
            this.mask = (NewBirthDeathSerialSamplingModel) mask;
        } else {
            throw new IllegalArgumentException();
        }
    }

    public void unmask() {
        mask = null;
    }

    // if a mask exists then use the mask's parameters instead (except for originTime and finalTimeInterval)
    NewBirthDeathSerialSamplingModel mask = null;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "A serially-sampled birth-death model with the possibility of treatment and sampling at present.";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(new Citation(
                new Author[]{
                        new Author("T", "Gernhard"),
                },
                "The conditioned reconstructed process",
                2008,
                "Journal of Theoretical Biology",
                253,
                769, 778,
                "10.1016/j.jtbi.2008.04.005"
        ));
    }
}