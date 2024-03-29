/*
 * PartitionClockModelTreeModelLink.java
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

package dr.app.beauti.options;

import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxa;
import dr.evomodelxml.tree.RateStatisticParser;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public class PartitionClockModelTreeModelLink extends PartitionOptions {

    private static final long serialVersionUID = -5594175190415743674L;

    private final PartitionClockModel model;
    private final PartitionTreeModel tree;

    public PartitionClockModelTreeModelLink(BeautiOptions options, PartitionClockModel model, PartitionTreeModel tree) {
//        super(options, model.getName() + "." + tree.getName());
        // clockModel and substModel have to be assigned before initModelParametersAndOpererators()
        super(options);
        this.model = model;
        this.tree = tree;
        initModelParametersAndOpererators();
    }

    protected void initModelParametersAndOpererators() {
//        {
//            final Parameter p = createParameter("branchRates.var", "autocorrelated lognormal relaxed clock rate variance ", PriorScaleType.LOG_VAR_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
//            p.priorType = PriorType.GAMMA_PRIOR;
//            p.shape = 1;
//            p.scale = 0.0001;
//        }
        createParameterGammaPrior("branchRates.var", "autocorrelated lognormal relaxed clock rate variance",
                PriorScaleType.LOG_VAR_SCALE, 0.1, 1, 0.0001, false);
        createParameter("branchRates.categories", "relaxed clock branch rate categories");
        createZeroOneParameter("branchRates.quantiles", "relaxed clock branch rate quantiles", 0.5);

//        {
//            final Parameter p = createParameter("treeModel.rootRate", "autocorrelated lognormal relaxed clock root rate", PriorScaleType.ROOT_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//            p.priorType = PriorType.GAMMA_PRIOR;
//            p.shape = 1;
//            p.scale = 0.0001;
//        }
        createParameterGammaPrior("treeModel.rootRate", "autocorrelated lognormal relaxed clock root rate",
                PriorScaleType.ROOT_RATE_SCALE, 1.0, 1, 0.0001, false);
        createParameterUniformPrior("treeModel.nodeRates", "autocorrelated lognormal relaxed clock non-root rates",
                PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);
        createParameterUniformPrior("treeModel.allRates", "autocorrelated lognormal relaxed clock all rates",
                PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Parameter.UNIFORM_MAX_BOUND);


        createScaleOperator("branchRates.var", demoTuning, rateWeights);

        createOperator("scaleRootRate", "treeModel.rootRate", "Scales root rate", "treeModel.rootRate", OperatorType.SCALE, demoTuning, rateWeights);
        createOperator("scaleOneRate", "treeModel.nodeRates", "Scales one non-root rate", "treeModel.nodeRates",
                OperatorType.SCALE, demoTuning, branchWeights);
        createOperator("scaleAllRates", "treeModel.allRates", "Scales all rates simultaneously", "treeModel.allRates",
                OperatorType.SCALE_ALL, demoTuning, rateWeights);
        createOperator("scaleAllRatesIndependently", "treeModel.nodeRates", "Scales all non-root rates independently", "treeModel.nodeRates",
                OperatorType.SCALE_INDEPENDENTLY, demoTuning, rateWeights);

        createOperator("swapBranchRateCategories", "branchRates.categories", "Performs a swap of branch rate categories",
                "branchRates.categories", OperatorType.SWAP, 1, branchWeights / 3);
//        createOperator("randomWalkBranchRateCategories", "branchRates.categories", "Performs an integer random walk of branch rate categories",
//                "branchRates.categories", OperatorType.INTEGER_RANDOM_WALK, 1, branchWeights / 3);
        createOperator("uniformBranchRateCategories", "branchRates.categories", "Performs an integer uniform draw of branch rate categories",
                "branchRates.categories", OperatorType.INTEGER_UNIFORM, 1, branchWeights / 3);

        createOperator("uniformBranchRateQuantiles", "branchRates.quantiles", "Performs an uniform draw of branch rate quantiles",
                "branchRates.quantiles", OperatorType.UNIFORM, 0, branchWeights);

        createUpDownOperator("upDownRateHeights", "Substitution rate and heights",
                "Scales substitution rates inversely to node heights of the tree", model.getParameter("clock.rate"),
                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, true, demoTuning, rateWeights);
        createUpDownOperator("upDownUCEDMeanHeights", "UCED mean and heights",
                "Scales UCED mean inversely to node heights of the tree", model.getParameter(ClockType.UCED_MEAN),
                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, true, demoTuning, rateWeights);
        createUpDownOperator("upDownUCLDMeanHeights", "UCLD mean and heights",
                "Scales UCLD mean inversely to node heights of the tree", model.getParameter(ClockType.UCLD_MEAN),
                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, true, demoTuning, rateWeights);
        createUpDownOperator("upDownUCGDMeanHeights", "UCGD mean and heights",
                "Scales UCGD mean inversely to node heights of the tree", model.getParameter(ClockType.UCGD_MEAN),
                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.UP_DOWN, true, demoTuning, rateWeights);

        // These should not have priors on as there will be priors on the clock model parameters already..

        // These are statistics which could have priors on...
        // #meanRate = #Relaxed Clock Model * #Tree Model
//        createNonNegativeStatistic("meanRate", "The mean rate of evolution over the whole tree");
        // #covariance = #Relaxed Clock Model * #Tree Model
//        createStatistic("covariance", "The covariance in rates of evolution on each lineage with their ancestral lineages");
        // #COEFFICIENT_OF_VARIATION = #Uncorrelated Clock Model
//        createNonNegativeStatistic(RateStatisticParser.COEFFICIENT_OF_VARIATION, "The variation in rate of evolution over the whole tree");

        createUpDownOperator("microsatUpDownRateHeights", "Substitution rate and heights",
                "Scales substitution rates inversely to node heights of the tree", model.getParameter("clock.rate"),
                tree.getParameter("treeModel.allInternalNodeHeights"), OperatorType.MICROSAT_UP_DOWN, true, demoTuning, branchWeights);
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectParameters(List<Parameter> params) {
//        setAvgRootAndRate();
//        getParameter("branchRates.categories");
//        getParameter("treeModel.rootRate");
//        getParameter("treeModel.nodeRates");
//        getParameter("treeModel.allRates");
//
//        if (options.hasData()) {
//            // if not fixed then do mutation rate move and up/down move
//            boolean fixed = !model.isEstimatedRate();
//
//            Parameter rateParam;
//
//            switch (model.getClockType()) {
//                case AUTOCORRELATED:
//                    rateParam = getParameter("treeModel.rootRate");
//                    rateParam.isFixed = fixed;
//                    if (!fixed) params.add(rateParam);
//
//                    params.add(getParameter("branchRates.var"));
//                    break;
//            }
//        }
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    public void selectOperators(List<Operator> ops) {
        if (options.hasData()) {
            // always have upDown(rate, allInternalNodeHeights), but when isEstimatedRate() = false, write nothing on up part (rate)
            Operator op;
            if (model.getDataType().getType() == DataType.MICRO_SAT) {
                if (model.getClockType() == ClockType.STRICT_CLOCK) {
                    op = getOperator("microsatUpDownRateHeights");
                    ops.add(op);
                } else {
                    throw new UnsupportedOperationException("Microsatellite only supports strict clock model");
                }

            } else {

                switch (model.getClockType()) {
                    case STRICT_CLOCK:
                        op = getOperator("upDownRateHeights");
                        ops.add(op);
                        break;

                    case UNCORRELATED:
                        switch (model.getClockDistributionType()) {

                            case LOGNORMAL:
                                op = getOperator("upDownUCLDMeanHeights");
                                ops.add(op);
                                break;
                            case GAMMA:
//                                throw new UnsupportedOperationException("Uncorrelated gamma model not implemented yet");
                                op = getOperator("upDownUCGDMeanHeights");
                                ops.add(op);
                            break;
                            case CAUCHY:
                                throw new UnsupportedOperationException("Uncorrelated Cauchy model not implemented yet");
//                            break;
                            case EXPONENTIAL:
                                op = getOperator("upDownUCEDMeanHeights");
                                ops.add(op);
                                break;
                        }
                        if (model.isContinuousQuantile()) {
                            ops.add(getOperator("uniformBranchRateQuantiles"));
                        } else {
                            ops.add(getOperator("swapBranchRateCategories"));
                            ops.add(getOperator("uniformBranchRateCategories"));
                        }
                        break;

                    case AUTOCORRELATED:

                        throw new UnsupportedOperationException("Autocorrelated relaxed clock model not implemented yet");

//                	if (options.clockModelOptions.getRateOptionClockModel() == FixRateType.RElATIVE_TO) {//&& model.isEstimatedRate()) {
//	                	ops.add(getOperator("scaleRootRate"));
//	                    ops.add(getOperator("scaleOneRate"));
//	                    ops.add(getOperator("scaleAllRates"));
//	                    ops.add(getOperator("scaleAllRatesIndependently"));
//	                    ops.add(getOperator("branchRates.var"));
////	                	ops.add(getOperator("upDownAllRatesHeights"));
//                	} else {
//                		ops.add(getOperator("scaleOneRate"));
//	                    ops.add(getOperator("scaleAllRatesIndependently"));
//	                    ops.add(getOperator("branchRates.var"));
//                	}
//                    break;

                    case RANDOM_LOCAL_CLOCK:
                    case FIXED_LOCAL_CLOCK:
                        op = getOperator("upDownRateHeights");
                        ops.add(op);

                        break;

                    default:
                        throw new IllegalArgumentException("Unknown clock model");
                }
            }
        }
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    public void selectStatistics(List<Parameter> params) {
    }

    /////////////////////////////////////////////////////////////
    public PartitionClockModel getPartitionClockModel() {
        return model;
    }

    public PartitionTreeModel getPartitionTreeTree() {
        return tree;
    }

    public String getPrefix() {
        return noDuplicatedPrefix(model.getPrefix(), tree.getPrefix());
    }

}