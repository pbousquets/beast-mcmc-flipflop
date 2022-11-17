package test.dr.evomodel.flipflop;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.AFsequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.sitemodel.GammaSiteModel;
import dr.evomodel.substmodel.FlipFlopModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.evomodelxml.sitemodel.GammaSiteModelParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.Parameter;
import dr.evomodel.treelikelihood.FlipFlopErrorModel;

import dr.inference.operators.*;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import junit.framework.TestCase;

/**
 * Test the flipflop pipeline
 *
 * @author Pablo Bousquets
 */
public class TestFlipFlops extends TestCase {

    interface Instance {

        int getNtips();
        int getNstates();

        double getStemCellParameter();
        double getDeltaParameter();
        double getEtaParameter();
        double getKappaParameter();

        double getGammaParameter();
        double getLambdaParameter();
        double getMuParameter();
        boolean getUseFrequencyModel();

        String[] getSequences();
        double[] getExpectedResult(int iTip);

    }

    Instance test0 = new Instance() {
        //DATA
        protected final int nTips=3;
        protected final String[] sequences = {
                "0.1, 0.4, 0.14, 0.12, 0.99, .95, 0.01",
                "0.2, 0.3, 0.13, 0.42, 0.29, .01, 0.91",
                "0.6, 0.36, 0.1, 0.62, 0.24, .15, 0.73"
        };

        //MODEL
        protected final int nCells=2;
        protected final int nStates = (int) (0.5 * (nCells + 1) * (nCells + 2));
        protected final double delta=0.05;
        protected final double eta=0.95;
        protected final double kappa=0.9;

        protected final double gamma=0.05;
        protected final double lambda=0.95;
        protected final double mu=0.05;
        protected final boolean useFreqModel = false;

        //EXPECTED
        protected final double [] expectedResults = { //TODO!
        };

        //GETTERS
        public int getNtips(){return this.nTips;};
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

        public double getGammaParameter(){
            return this.gamma;
        };
        public double getLambdaParameter(){
            return this.lambda;
        };
        public double getMuParameter(){
            return this.mu;
        };
        public boolean getUseFrequencyModel(){
            return this.useFreqModel;
        };

        public String[] getSequences(){
            return this.sequences; //Implementation for only one patter/locus
        }
        public double[] getExpectedResult(int iTip) {
            return(Arrays.copyOfRange(expectedResults,iTip*nStates,iTip*nStates+nStates));
        }
    };

    Instance[] all = {test0}; //Add more instances of tests here

    public void testFlipFlops() {
        for (Instance test : all) {
            // Get the error model parameters
            double delta = test.getDeltaParameter();
            double eta = test.getEtaParameter();
            double kappa = test.getKappaParameter();
            double cells = test.getStemCellParameter();
            int Ntips = test.getNtips();
            int nstates = test.getNstates();
            AFsequence afseq = new AFsequence(nstates); // Needed to initialize the models later

            //Get the flipflop model parameters
            Parameter stemCellParameter = new Parameter.Default(1, test.getStemCellParameter());
            Parameter gammaParam = new Parameter.Default(1, test.getGammaParameter());
            Parameter lambdaParam = new Parameter.Default(1, test.getLambdaParameter());
            Parameter muParam = new Parameter.Default(1, test.getMuParameter());
            boolean useFreqModel = test.getUseFrequencyModel();

            //Create pi freqs and modify them so not all of them have the same value
            double[] freqs = new double[nstates];
            Arrays.fill(freqs, (double) 1/nstates);
            freqs[1] = freqs[1]+0.2;
            freqs[2] = freqs[2]-0.15;
            freqs[3] = freqs[3]-0.05;

            // Create PatternList to be able to use the FlipFlopErrorModel. Just adding C to each tip's number (short for Crypt)
            //Taxa
            Taxa taxonlist = new Taxa();
            for (int iTip=0;iTip<Ntips;iTip++){
                taxonlist.addTaxon(new Taxon("C" + String.valueOf(iTip)));
            }

            //Patterns
            List<int[]> seqlist = new ArrayList<int[]>();
            taxonlist = new Taxa();
            for (int idx = 0; idx < test.getSequences().length; idx ++){
                AFsequence seq = new AFsequence(test.getSequences()[idx]);
                seq.setTaxon(new Taxon("S"+idx));
                seqlist.add(seq.getSequence());
                taxonlist.addTaxon(seq.getTaxon());

            }

            Patterns patterns = new Patterns(new AFsequence(nstates), taxonlist);
            for (int site = 0; site < (seqlist.get(0)).length; site++){
                int[] current_pattern = new int[seqlist.size()];
                for (int seq_index = 0; seq_index < seqlist.size(); seq_index++){
                    current_pattern[seq_index] = seqlist.get(seq_index)[site];
                }
                patterns.addPattern(current_pattern);
            }
            patterns.setId("test");

            System.out.println("  - Taxa count = " + patterns.getTaxonCount());
            System.out.println("  - Site count = " + patterns.getPatternCount());
            System.out.println("  - Expected states = " + patterns.getStateCount());

            //Tree simulation. We do not care about the simulation parameters, we just need a tree object with the proper number of external nodes (tips) for the super(FlipFlopErrorModel) to initalize the tip states properly.
            // We do not care about anything else
            CoalescentSimulator simulator = new CoalescentSimulator();
            DemographicModel constantPop = new ConstantPopulationModel(new Parameter.Default(50.0), dr.evolution.util.Units.Type.YEARS);
            Tree theTree=simulator.simulateTree(taxonlist,constantPop);
            TreeModel treeModel = new TreeModel(theTree);

            //Creating the error model
            FlipFlopErrorModel errorModel= new FlipFlopErrorModel(taxonlist, new Taxa(), new Parameter.Default(cells), new Parameter.Default(delta),new Parameter.Default(eta),new Parameter.Default(kappa));

            //Initializing the error model
            errorModel.setTree(theTree);
            for (int iTip=0;iTip<patterns.getTaxonCount();iTip++) {
                errorModel.setStates(patterns,iTip,iTip,patterns.getTaxonId(iTip));
            }

            FrequencyModel freqModel = new FrequencyModel(afseq, freqs);

            FlipFlopModel model = new FlipFlopModel("test", stemCellParameter, gammaParam, lambdaParam, muParam, useFreqModel, freqModel);

            //siteModel
            GammaSiteModel siteModel = new GammaSiteModel(model);

            TreeLikelihood treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, null,
                    false, false, true, false, false);
            treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);

            System.out.println("\n\nTHE INITIAL LIKELIHOOD IS: " + treeLikelihood.getLogLikelihood() + "\n\n");

            /*
            // Operators
            OperatorSchedule schedule = new SimpleOperatorSchedule();

            MCMCOperator operator = new ScaleOperator(kappaParam, 0.5);
            operator.setWeight(1.0);
            schedule.addOperator(operator);

            Parameter rootHeight = treeModel.getRootHeightParameter();
            String TREE_HEIGHT = TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT;
            rootHeight.setId(TREE_HEIGHT);
            operator = new ScaleOperator(rootHeight, 0.5);
            operator.setWeight(1.0);
            schedule.addOperator(operator);

            Parameter internalHeights = treeModel.createNodeHeightsParameter(false, true, false);
            operator = new UniformOperator(internalHeights, 10.0);
            schedule.addOperator(operator);

            operator = new SubtreeSlideOperator(treeModel, 1, 1, true, false, false, false, CoercionMode.COERCION_ON);
            schedule.addOperator(operator);

            operator = new ExchangeOperator(ExchangeOperator.NARROW, treeModel, 1.0);
//        operator.doOperation();
            schedule.addOperator(operator);

            operator = new ExchangeOperator(ExchangeOperator.WIDE, treeModel, 1.0);
//        operator.doOperation();
            schedule.addOperator(operator);

            operator = new WilsonBalding(treeModel, 1.0);
//        operator.doOperation();
            schedule.addOperator(operator);

            // Log
            ArrayLogFormatter formatter = new ArrayLogFormatter(false);

            MCLogger[] loggers = new MCLogger[2];
            loggers[0] = new MCLogger(formatter, 1, false);
            loggers[0].add(treeLikelihood);
            loggers[0].add(rootHeight);
            //loggers[0].add(kappaParam);

            loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 1, false);
            loggers[1].add(treeLikelihood);
            loggers[1].add(rootHeight);
            //loggers[1].add(kappaParam);

            // MCMC
            MCMC mcmc = new MCMC("mcmc1");
            MCMCOptions options = new MCMCOptions(2);

            mcmc.setShowOperatorAnalysis(true);
            mcmc.init(options, treeLikelihood, schedule, loggers);
            mcmc.run();

            // time
            System.out.println(mcmc.getTimer().toString());

            // Tracer
            List<Trace> traces = formatter.getTraces();
            ArrayTraceList traceList = new ArrayTraceList("MCMCTest", traces, 0);

            for (int i = 1; i < traces.size(); i++) {
                traceList.analyseTrace(i);
            }
        */

        }
    };

}
