package test.dr.evomodel.flipflop;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dr.app.treestat.TreeStatData;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.AFsequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.tree.StarTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.util.FlipFlopUtils;
import dr.evomodel.treelikelihood.FlipFlopErrorModel;

import junit.framework.TestCase;

/**
 * Test the flipflop error model
 *
 * @author Pablo Bousquets
 */
public class TestFlipFlopErrorModel extends TestCase {

    interface Instance {

        int getNtips();
        int getNcells();
        int getNstates();

        double getStemCellParameter();
        double getDeltaParameter();
        double getEtaParameter();
        double getKappaParameter();

        double[] getSequence();

        String getSequenceString (int iTip);

        double[] getExpectedResult();

        double[] getExpectedResult(int iTip);

    }

    Instance test0 = new Instance() {
        //DATA
        protected final int nTips=7;
        protected final double[] sequence={0.1, 0.4, 0.14, 0.12, 0.99, .05, 0.01};

        //MODEL
        protected final int nCells=2;
        protected final int nStates=2 * nCells + 1;
        protected final double delta=0.05;
        protected final double eta=0.95;
        protected final double kappa=0.9;

        //EXPECTED

        protected final double [] expectedIdealAlpha= {0.045 , 0.2475, 0.45  , 0.6525, 0.855};

        protected final double [] expectedIdealBeta= {0.855 , 0.6525, 0.45  , 0.2475, 0.045};
        protected final double [] expectedResults = {
                -0.8984707443537294, 0.21490279494496956, 0.03657211353834405, -0.6749731588762, -2.678222651996067,
                -2.163589418547541, -0.6873845867301609, -0.5028839756181054, -0.8515979555139683, -2.492016156115154,
                -1.2132096860858579, -0.02249438812081994, -0.12348331086110648, -0.7576868246093116, -2.6835945590628394,
                -1.0693292669934076, 0.08551514086610923, -0.051344671979698886, -0.7214190758334251, -2.6831977003924745,
                -2.4353690352250608, 0.05345425063647635, 1.2505733157926897, 1.9144777899409842, 1.2866780433839557,
                -0.24435493400317204, 0.7177076889249079, 0.3880660911476625, -0.47479009763750163, -2.629350507127988,
                1.2866780433839566, 1.9144777899409862, 1.2505733157926906, 0.05345425063647652, -2.4353690352250603
        };

        //GETTERS
        public int getNtips(){return this.nTips;};
        public int getNcells(){return this.nCells;};
        public int getNstates(){return this.nStates;};
        public double getStemCellParameter(){
            return this.nCells;
        };
        public double getDeltaParameter(){
            return this.delta;
        };
        public double getEtaParameter(){
            return this.eta;
        };
        public double getKappaParameter(){
            return this.kappa;
        };
        public double[] getSequence(){
            return this.sequence;
        };

        public String getSequenceString(int iTip){
            return(String.valueOf(this.sequence[iTip])); //Implementation for only one patter/locus
        }
        public double[] getExpectedIdealAlpha() {
            return this.expectedIdealAlpha;
        }
        public double[] getExpectedIdealBeta() {
            return this.expectedIdealBeta;
        }

        public double[] getExpectedResult() {
            return(this.expectedResults);
        }

        public double[] getExpectedResult(int iTip) {
            return(Arrays.copyOfRange(expectedResults,iTip*nStates,iTip*nStates+nStates));
        }

    };

    Instance[] all = {test0}; //Add more instances of tests here

    public void testErrorModel() {
        for (Instance test : all) {
            // Create the parameters
            double delta = test.getDeltaParameter();
            double eta = test.getEtaParameter();
            double kappa = test.getKappaParameter();
            double cells = test.getStemCellParameter();

            // Get the pattern
            double[] seq = test.getSequence();

            // Compute the partials!
            int total_states = 2 * (int) cells + 1;
            double[] ideal_beta = new double [total_states];
            double[] beta = new double [total_states];
            double[] transformed_alpha = new double [total_states];
            double[] transformed_beta = new double [total_states];

            for (int state = 0; state < total_states; state ++){
                ideal_beta[state] = ((double) state) / (total_states - 1); // We'll have 0/N methylated alleles, 1/N allele, 2/N alleles, ... N/N methylated alleles (considering here N = total alleles)
                beta[state] = (eta - delta) * ideal_beta[state] + delta; // Rescale the expected peaks to linear scale

                //Convert mean/dispersion parameterization of a beta distribution
                transformed_alpha[state] = beta[state] * kappa;
                transformed_beta[state] =  (1 - beta[state]) * kappa;
            }

            //Compute log(p(y|z))
            double[] partials = FlipFlopUtils.beta_lpdf(seq, transformed_alpha, transformed_beta); //Format: partials per site per state: [site1state1,site1state2,site1state3,...,site2state1,site2state2,site2state3,...]
            System.out.println(Arrays.toString(partials));

            double[] expected = test.getExpectedResult();
            System.out.println(Arrays.toString(expected));
            assertEquals(expected.length, partials.length, 1e-15);
            
            for (int kk = 0; kk < expected.length; kk+=1) {
                assertEquals(expected[kk], partials[kk], 1e-5);
            }

        }
    };

    public void testErrorModelClass() {
        for (Instance test : all) {
            // Get the parameters
            double delta = test.getDeltaParameter();
            double eta = test.getEtaParameter();
            double kappa = test.getKappaParameter();
            double cells = test.getStemCellParameter();
            int Ntips = test.getNtips();
            int Nstates = test.getNstates();

            // Create PatternList to be able to use the FlipFlopErrorModel
            //Taxa
            Taxa taxonlist = new Taxa();
            for (int iTip=0;iTip<Ntips;iTip++){
                taxonlist.addTaxon(new Taxon("C" + String.valueOf(iTip)));
            }

            //Pattern
            Patterns patterns = new Patterns(new AFsequence(Nstates), taxonlist);
            List<int[]> seqlist = new ArrayList<int[]>();
            for (int iTip=0;iTip<Ntips;iTip++) {
                seqlist.add(new AFsequence(test.getSequenceString(iTip)).getSequence());
            }
            for (int site = 0; site < (seqlist.get(0)).length; site++){
                int[] current_pattern = new int[seqlist.size()];
                for (int seq_index = 0; seq_index < seqlist.size(); seq_index++){
                    current_pattern[seq_index] = seqlist.get(seq_index)[site];
                }
                patterns.addPattern(current_pattern);
            }

            //We will simulate a tree here. We only need the tree for the super(FlipFlopErrorModel) to set the tip partials properly and we do not really care about anything else
            CoalescentSimulator simulator = new CoalescentSimulator();
            DemographicModel constantPop = new ConstantPopulationModel(new Parameter.Default(50.0), dr.evolution.util.Units.Type.YEARS);
            Tree theTree=simulator.simulateTree(taxonlist,constantPop);

            //Creating the error model
            FlipFlopErrorModel errorModel= new FlipFlopErrorModel(taxonlist, new Taxa() ,new Parameter.Default(cells), new Parameter.Default(delta),new Parameter.Default(eta),new Parameter.Default(kappa));
            errorModel.setTree(theTree);

            //Initializing the error model
            for (int iTip=0;iTip<Ntips;iTip++) {
                errorModel.setStates(patterns,iTip,iTip,patterns.getTaxonId(iTip));
            }

            for (int iTip=0;iTip<Ntips;iTip++){
                double [] expectedPartials = test.getExpectedResult(iTip);
                double [] calculatedPartials = new double [Nstates];
                errorModel.getTipPartials(iTip,calculatedPartials);
                for (int kk = 0; kk < Nstates; kk++){
                    assertEquals(expectedPartials[kk], calculatedPartials[kk], 1e-5);
                }
            }
        }
    };

}
