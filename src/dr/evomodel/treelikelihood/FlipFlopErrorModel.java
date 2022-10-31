/*
 * SequenceErrorModel.java
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

package dr.evomodel.treelikelihood;

import dr.evolution.util.TaxonList;
import dr.math.GammaFunction;
import dr.evomodelxml.treelikelihood.FlipFlopErrorModelParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.FlipFlopUtils;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 */
public class FlipFlopErrorModel extends TipStatesModel implements Citable {

    public FlipFlopErrorModel(TaxonList includeTaxa, TaxonList excludeTaxa,
                              Parameter stemCellParameter,
                              Parameter deltaParameter,
                              Parameter etaParameter,
                              Parameter kappaParameter) {
        super(FlipFlopErrorModelParser.AFSEQUENCE_ERROR_MODEL, includeTaxa, excludeTaxa);

        this.stemCellParameter = stemCellParameter; //No need to evaluate if is null, already done in parser
        addVariable(this.stemCellParameter);

        if (deltaParameter != null) {
            this.deltaParameter = deltaParameter;
            addVariable(deltaParameter);
        } else {
            this.deltaParameter = null;
        }

        if (etaParameter != null) {
            this.etaParameter = etaParameter;
            addVariable(etaParameter);
        } else {
            this.etaParameter = null;
        }

        if (kappaParameter != null) {
            this.kappaParameter = kappaParameter;
            addVariable(kappaParameter);
        } else {
            this.kappaParameter = null;
        }

        this.stateCount = 2 * (int) stemCellParameter.getParameterValue(0) + 1;
    }

    /*
    protected void taxaChanged() {
        if (indicatorParameter != null && indicatorParameter.getDimension() <= 1) {
            this.indicatorParameter.setDimension(tree.getExternalNodeCount());
        }
    }
*/

    @Override
    public Type getModelType() {
        return Type.PARTIALS;
    }

    @Override
    public void getTipStates(int nodeIndex, int[] tipStates) {
        throw new IllegalArgumentException("This model emits only tip partials");
    }

    @Override
    public void getTipPartials(int nodeIndex, double[] partials) {
        /* Note: the partials are generated in the abstracttreelikelihood class as:
            double[] partials = new double[patternCount * stateCount];
         */

        int[] states = this.states[nodeIndex];
        double[] double_states = FlipFlopUtils.remapDigitized(states, 200); // Turn the states into the 0,1 original range

        double delta = deltaParameter.getParameterValue(0);
        double eta = etaParameter.getParameterValue(0);
        double kappa = kappaParameter.getParameterValue(0);

        int total_states = this.stateCount;
        double[] ideal_beta = new double [total_states];
        double[] beta = new double [total_states];
        double[] transformed_alpha = new double [total_states];
        double[] transformed_beta = new double [total_states];

        for (int state = 0; state < total_states; state ++){
            ideal_beta[state] = state/total_states; // We'll have 0/N methylated alleles, 1/N allele, 2/N alleles, ... N/N methylated alleles (considering here N = total alleles)
            beta[state] = (eta - delta) * ideal_beta[state] + delta; // Rescale the expected peaks to linear scale

            //Convert mean/dispersion parameterization of a beta distribution
            transformed_alpha[state] = beta[state] * kappa;
            transformed_beta[state] =  (1 - beta[state]) * kappa;
        }

        //Loop through each of the possible Z=2S+1 peaks and calculate the log-pdf of the nth beta value in y[N] for the zth beta distribution
        int N = double_states.length; //number of sites: y.len in the original function
        int Z = transformed_alpha.length; //number of states: alpha.len in the original function
        double lgamma_alpha;
        double lgamma_beta;
        double lgamma_alphaplusbeta;

        for (int z = 0; z < Z; z++) {  //Precompute the log-gamma constants. Iterate through states
            lgamma_alpha = GammaFunction.lnGamma(transformed_alpha[z]);
            lgamma_beta = GammaFunction.lnGamma(transformed_beta[z]);
            lgamma_alphaplusbeta = GammaFunction.lnGamma(transformed_alpha[z] + transformed_beta[z]);

            for (int n = 0; n < N; n++) { //Calculate the log pdf of the beta distribution for the nth dataproint drawn from the zth peak. Iterate through sites
                //Partials format is partials per site per state: [site1state1,site1state2,site1state3,...,site2state1,site2state2,site2state3,...]
                partials[n * Z + z] = ((transformed_alpha[z] - 1) * Math.log(double_states[n]) + (transformed_beta[z] - 1) * StrictMath.log1p(-double_states[n]) - lgamma_alpha - lgamma_beta + lgamma_alphaplusbeta);
            }
        }
    }

    private final Parameter stemCellParameter;
    private final Parameter deltaParameter;
    private final Parameter etaParameter;
    private final Parameter kappaParameter;

    @Override
    protected void taxaChanged(){
        throw new RuntimeException("Unexpected behaviour. taxaChanged method has been activated, but it shouldn't.");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Sequence error model";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                        new Author[]{
                                new Author("A", "Rambaut"),
                                new Author("SYW", "Ho"),
                                new Author("AJ", "Drummond"),
                                new Author("B", "Shapiro"),
                        },
                        "Accommodating the effect of ancient DNA damage on inferences of demographic histories",
                        2008,
                        "Mol Biol Evol",
                        26,
                        245, 248,
                        "10.1093/molbev/msn256"
                ),
                new Citation(
                        new Author[]{
                                new Author("J", "Felsenstein"),
                        },
                        "Inferring Phylogenies",
                        2004,
                        "Sinauer Associates",
                        ""
                ));
    }
}