/*
 * AlleleFraction.java
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

package dr.evolution.datatype;
import java.util.ArrayList;

/**
 * @author Pablo Bousquets
 *
 * Allele frequency data type
 */
public class AlleleFraction extends DataType {

    public static final String DESCRIPTION = "allele_fraction";
    public static int UNKNOWN_STATE_LENGTH = -1;
    private String name;

    public static final AlleleFraction INSTANCE = new AlleleFraction();

    public AlleleFraction() {}

    @Override
    public char[] getValidChars() {
        return null;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    /**
     * @return the description of the data type
     */
    public String getDescription() {
        return DESCRIPTION;
    }

    public int getType(){
        return ALLELE_FRAC;
    }

}
