/*
 * ProductStatistic.java
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

package dr.inference.model;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ProductStatistic.java,v 1.2 2005/05/24 20:26:00 rambaut Exp $
 */
public class ProductStatistic extends Statistic.Abstract {

    private int dimension = 0;
    private final boolean elementwise;

    public ProductStatistic(String name, boolean elementwise) {
        super(name);
        this.elementwise = elementwise;
    }

    public void addStatistic(Statistic statistic) {
        if (!elementwise) {
            if (dimension == 0) {
                dimension = statistic.getDimension();
            } else if (dimension != statistic.getDimension()) {
                throw new IllegalArgumentException();
            }
        } else {
            dimension = 1;
        }
        statistics.add(statistic);
    }

    public int getDimension() {
        return dimension;
    }

    /**
     * @return mean of contained statistics
     */
    public double getStatisticValue(int dim) {

        double product = 1.0;

        for (Statistic statistic : statistics) {
            if (elementwise) {
                for (int j = 0; j < statistic.getDimension(); j++) {
                    product *= statistic.getStatisticValue(j);
                }
            } else {
                product *= statistic.getStatisticValue(dim);
            }
        }

        return product;
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private final List<Statistic> statistics = new ArrayList<Statistic>();
}
