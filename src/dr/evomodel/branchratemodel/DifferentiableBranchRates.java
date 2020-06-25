/*
 * BranchRates.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 * @author Alexander Fisher
 */
public interface DifferentiableBranchRates  extends BranchRates {

    double getBranchRateDifferential(Tree tree, NodeRef node);

    double getBranchRateSecondDifferential(Tree tree, NodeRef node);

    Parameter getRateParameter();

    int getParameterIndexFromNode(NodeRef node);

    ArbitraryBranchRates.BranchRateTransform getTransform(); // TODO Should remove
}
