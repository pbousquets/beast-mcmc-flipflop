/*
 * MicroSatModelType.java
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

package dr.app.beauti.types;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class MicroSatModelType {

    public enum RateProportionality {
        EQUAL_RATE("Equal Rate"),
        PROPORTIONAL_RATE("Proportional Rate"),
        ASYM_QUAD("Asymetric Quadratic Rate");

        RateProportionality(String displayName) {
            this.displayName = displayName;
        }

        public String toString() {
            return displayName;
        }

        private String displayName;
    }

    public enum MutationalBias {
        UNBIASED("Unbiased"),
        CONSTANT_BIAS("Constant Bias"),
        LINEAR_BIAS("Linear Bias");

        MutationalBias(String displayName) {
            this.displayName = displayName;
        }

        public String toString() {
            return displayName;
        }

        private String displayName;
    }

    public enum Phase {
        ONE_PHASE("One Phase"),
        TWO_PHASE("Two Phase"),
        TWO_PHASE_STAR("Two Phase *");

        Phase(String displayName) {
            this.displayName = displayName;
        }

        public String toString() {
            return displayName;
        }

        private String displayName;
    }

}
