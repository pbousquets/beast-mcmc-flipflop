/*
 * Blosum62.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.AminoAcids;
import dr.util.Author;
import dr.util.Citation;

import java.util.*;

/**
 * BLOSUM62 model of amino acid evolution
 * Henikoff, S., and J. G. Henikoff. 1992. PNAS USA 89:10915-10919.
 *
 * @version $Id: Blosum62.java,v 1.3 2005/05/24 20:25:58 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Blosum62 extends EmpiricalRateMatrix.AbstractAminoAcid {
	
	public static final Blosum62 INSTANCE = new Blosum62();

	// The rates below are specified assuming that the amino acids are in this order:
	// ARNDCQEGHILKMFPSTWYV
	// but the AminoAcids dataType wants them in this order:
	// ACDEFGHIKLMNPQRSTVWY
	// This is solved by calling the setEmpiricalRates and setEmpiricalFrequencies methods
	private Blosum62() { super("blosum62");

		int n = AminoAcids.INSTANCE.getStateCount();
		
		double[][] rate = new double[n][n];
		
		// Q matrix
		rate[0][1]=7.3579038969751e-01; 	rate[0][2]=4.8539105546575e-01;
		rate[0][3]=5.4316182089867e-01; 	rate[0][4]=1.4599953104700e+00;
		rate[0][5]=1.1997057046020e+00; 	rate[0][6]=1.1709490427999e+00;
		rate[0][7]=1.9558835749595e+00; 	rate[0][8]=7.1624144499779e-01;
		rate[0][9]=6.0589900368677e-01; 	rate[0][10]=8.0001653051838e-01;
		rate[0][11]=1.2952012667833e+00; 	rate[0][12]=1.2537582666635e+00;
		rate[0][13]=4.9296467974759e-01; 	rate[0][14]=1.1732759009239e+00;
		rate[0][15]=4.3250926870566e+00; 	rate[0][16]=1.7291780194850e+00;
		rate[0][17]=4.6583936772479e-01; 	rate[0][18]=7.1820669758623e-01;
		rate[0][19]=2.1877745220045e+00; 
	
		rate[1][2]=1.2974467051337e+00; 	rate[1][3]=5.0096440855513e-01;
		rate[1][4]=2.2782657420895e-01; 	rate[1][5]=3.0208336100636e+00;
		rate[1][6]=1.3605741904203e+00; 	rate[1][7]=4.1876330851753e-01;
		rate[1][8]=1.4561411663360e+00; 	rate[1][9]=2.3203644514174e-01;
		rate[1][10]=6.2271166969249e-01; 	rate[1][11]=5.4111151414889e+00;
		rate[1][12]=9.8369298745695e-01; 	rate[1][13]=3.7164469320875e-01;
		rate[1][14]=4.4813366171831e-01; 	rate[1][15]=1.1227831042096e+00;
		rate[1][16]=9.1466595456337e-01; 	rate[1][17]=4.2638231012175e-01;
		rate[1][18]=7.2051744121611e-01; 	rate[1][19]=4.3838834377202e-01;

		rate[2][3]=3.1801000482161e+00; 	rate[2][4]=3.9735894989702e-01;
		rate[2][5]=1.8392161469920e+00; 	rate[2][6]=1.2404885086396e+00;
		rate[2][7]=1.3558723444845e+00; 	rate[2][8]=2.4145014342081e+00;
		rate[2][9]=2.8301732627800e-01; 	rate[2][10]=2.1188815961519e-01;
		rate[2][11]=1.5931370434574e+00; 	rate[2][12]=6.4844127878707e-01;
		rate[2][13]=3.5486124922252e-01; 	rate[2][14]=4.9488704370192e-01;
		rate[2][15]=2.9041016564560e+00; 	rate[2][16]=1.8981736345332e+00;
		rate[2][17]=1.9148204624678e-01; 	rate[2][18]=5.3822251903674e-01;
		rate[2][19]=3.1285879799342e-01; 
	
		rate[3][4]=2.4083661480204e-01; 	rate[3][5]=1.1909457033960e+00;
		rate[3][6]=3.7616252083685e+00; 	rate[3][7]=7.9847324896839e-01;
		rate[3][8]=7.7814266402188e-01; 	rate[3][9]=4.1855573246161e-01;
		rate[3][10]=2.1813157759360e-01; 	rate[3][11]=1.0324479249521e+00;
		rate[3][12]=2.2262189795786e-01; 	rate[3][13]=2.8173069420651e-01;
		rate[3][14]=7.3062827299842e-01; 	rate[3][15]=1.5827541420653e+00;
		rate[3][16]=9.3418750943056e-01; 	rate[3][17]=1.4534504627853e-01;
		rate[3][18]=2.6142220896504e-01; 	rate[3][19]=2.5812928941763e-01;

		rate[4][5]=3.2980150463028e-01; 	rate[4][6]=1.4074889181440e-01;
		rate[4][7]=4.1820319228376e-01; 	rate[4][8]=3.5405810983129e-01;
		rate[4][9]=7.7489402279418e-01; 	rate[4][10]=8.3184264014158e-01;
		rate[4][11]=2.8507880090648e-01; 	rate[4][12]=7.6768882347954e-01;
		rate[4][13]=4.4133747118660e-01; 	rate[4][14]=3.5600849876863e-01;
		rate[4][15]=1.1971884150942e+00; 	rate[4][16]=1.1198313585160e+00;
		rate[4][17]=5.2766441887169e-01; 	rate[4][18]=4.7023773369610e-01;
		rate[4][19]=1.1163524786062e+00; 
	
		rate[5][6]=5.5289191779282e+00; 	rate[5][7]=6.0984630538281e-01;
		rate[5][8]=2.4353411311401e+00; 	rate[5][9]=2.3620245120365e-01;
		rate[5][10]=5.8073709318144e-01; 	rate[5][11]=3.9452776745146e+00;
		rate[5][12]=2.4948960771127e+00; 	rate[5][13]=1.4435695975031e-01;
		rate[5][14]=8.5857057567418e-01; 	rate[5][15]=1.9348709245965e+00;
		rate[5][16]=1.2774802945956e+00; 	rate[5][17]=7.5865380864172e-01;
		rate[5][18]=9.5898974285014e-01; 	rate[5][19]=5.3078579012486e-01;

		rate[6][7]=4.2357999217628e-01; 	rate[6][8]=1.6268910569817e+00;
		rate[6][9]=1.8684804693170e-01; 	rate[6][10]=3.7262517508685e-01;
		rate[6][11]=2.8024271516787e+00; 	rate[6][12]=5.5541539747043e-01;
		rate[6][13]=2.9140908416530e-01; 	rate[6][14]=9.2656393484598e-01;
		rate[6][15]=1.7698932389373e+00; 	rate[6][16]=1.0710972360073e+00;
		rate[6][17]=4.0763564893830e-01; 	rate[6][18]=5.9671930034577e-01;
		rate[6][19]=5.2425384633796e-01; 
	
		rate[7][8]=5.3985912495418e-01; 	rate[7][9]=1.8929629237636e-01;
		rate[7][10]=2.1772115923623e-01; 	rate[7][11]=7.5204244030271e-01;
		rate[7][12]=4.5943617357855e-01; 	rate[7][13]=3.6816646445253e-01;
		rate[7][14]=5.0408659952683e-01; 	rate[7][15]=1.5093262532236e+00;
		rate[7][16]=6.4143601140497e-01; 	rate[7][17]=5.0835892463812e-01;
		rate[7][18]=3.0805573703500e-01; 	rate[7][19]=2.5334079019018e-01;

		rate[8][9]=2.5271844788492e-01; 	rate[8][10]=3.4807220979697e-01;
		rate[8][11]=1.0225070358890e+00; 	rate[8][12]=9.8431152535870e-01;
		rate[8][13]=7.1453370392764e-01; 	rate[8][14]=5.2700733915060e-01;
		rate[8][15]=1.1170297629105e+00; 	rate[8][16]=5.8540709022472e-01;
		rate[8][17]=3.0124860078016e-01; 	rate[8][18]=4.2189539693890e+00;
		rate[8][19]=2.0155597175031e-01; 
	
		rate[9][10]=3.8909637733035e+00; 	rate[9][11]=4.0619358664202e-01;
		rate[9][12]=3.3647977631042e+00; 	rate[9][13]=1.5173593259539e+00;
		rate[9][14]=3.8835540920564e-01; 	rate[9][15]=3.5754441245967e-01;
		rate[9][16]=1.1790911972601e+00; 	rate[9][17]=3.4198578754023e-01;
		rate[9][18]=6.7461709322842e-01; 	rate[9][19]=8.3118394054582e+00;

		rate[10][11]=4.4557027426059e-01; 	rate[10][12]=6.0305593795716e+00;
		rate[10][13]=2.0648397032375e+00; 	rate[10][14]=3.7455568747097e-01;
		rate[10][15]=3.5296918452729e-01; 	rate[10][16]=9.1525985769421e-01;
		rate[10][17]=6.9147463459998e-01; 	rate[10][18]=8.1124585632307e-01;
		rate[10][19]=2.2314056889131e+00; 
	
		rate[11][12]=1.0730611843319e+00; 	rate[11][13]=2.6692475051102e-01;
		rate[11][14]=1.0473834507215e+00; 	rate[11][15]=1.7521659178195e+00;
		rate[11][16]=1.3038752007987e+00; 	rate[11][17]=3.3224304063396e-01;
		rate[11][18]=7.1799348690032e-01; 	rate[11][19]=4.9813847530407e-01;

		rate[12][13]=1.7738551688305e+00; 	rate[12][14]=4.5412362510273e-01;
		rate[12][15]=9.1872341574605e-01; 	rate[12][16]=1.4885480537218e+00;
		rate[12][17]=8.8810109815193e-01; 	rate[12][18]=9.5168216224591e-01;
		rate[12][19]=2.5758507553153e+00; 
	
		rate[13][14]=2.3359790962888e-01; 	rate[13][15]=5.4002764482413e-01;
		rate[13][16]=4.8820611879305e-01; 	rate[13][17]=2.0743248934965e+00;
		rate[13][18]=6.7472604308008e+00; 	rate[13][19]=8.3811961017754e-01;

		rate[14][15]=1.1691295777157e+00; 	rate[14][16]=1.0054516831488e+00;
		rate[14][17]=2.5221483002727e-01; 	rate[14][18]=3.6940531935451e-01;
		rate[14][19]=4.9690841067567e-01; 
	
		rate[15][16]=5.1515562922704e+00; 	rate[15][17]=3.8792562209837e-01;
		rate[15][18]=7.9675152076106e-01; 	rate[15][19]=5.6192545744165e-01;

		rate[16][17]=5.1312812689059e-01; 	rate[16][18]=8.0101024319939e-01;
		rate[16][19]=2.2530740511763e+00; 
	
		rate[17][18]=4.0544190065580e+00; 	rate[17][19]=2.6650873142646e-01;

		rate[18][19]=1.0000000000000e+00; 
		
		setEmpiricalRates(rate, "ARNDCQEGHILKMFPSTWYV");

		double[] f = new double[n];
		f[0]=0.074; f[1]=0.052; f[2]=0.045; f[3]=0.054;
		f[4]=0.025; f[5]=0.034; f[6]=0.054; f[7]=0.074;
		f[8]=0.026; f[9]=0.068; f[10]=0.099; f[11]=0.058;
		f[12]=0.025; f[13]=0.047; f[14]=0.039; f[15]=0.057;
		f[16]=0.051; f[17]=0.013; f[18]=0.032; f[19]=0.073;
		
		setEmpiricalFrequencies(f, "ARNDCQEGHILKMFPSTWYV");
	}

	@Override
	public Citation.Category getCategory() {
		return Citation.Category.SUBSTITUTION_MODELS;
	}

	@Override
	public String getDescription() {
		return "Blosum62 amino acid substitution model";
	}

	@Override
	public List<Citation> getCitations() {
		return Collections.singletonList(CITATION);
	}

	public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("S", "Henikoff"),
                    new Author("JG", "Henikoff")
            },
            "Amino acid substitution matrices from protein blocks",
            1992,
            "Proc Natl Acad Sci, USA",
            89,
            10915, 10919,
            Citation.Status.PUBLISHED
    );
}
