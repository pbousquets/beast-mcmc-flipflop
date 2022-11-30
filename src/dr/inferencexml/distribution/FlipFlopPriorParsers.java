/*
 * PriorParsers.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.distribution;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Statistic;
import dr.math.distributions.*;
import dr.xml.*;

/**
 */
public class FlipFlopPriorParsers extends PriorParsers{
    public final static boolean DEBUG = false;

    public static final String HALF_NORMAL_PRIOR = "halfNormalPrior";

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser HALF_NORMAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return HALF_NORMAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double stdev = xo.getDoubleAttribute(STDEV);

            DistributionLikelihood likelihood = new DistributionLikelihood(new HalfNormalDistribution(mean, stdev));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MEAN),
                AttributeRule.newDoubleRule(STDEV),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE)
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given half-normal distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };
}
