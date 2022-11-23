package test.dr.evomodel.flipflop;
import java.util.ArrayList;
import java.util.List;
import java.util.Arrays;

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
import junit.framework.TestCase;

/**
 * Test the flipflop pipeline
 *
 * @author Pablo Bousquets
 */
public class TestFlipFlops extends TestCase {
    class Test {
        public Test(String name, String[] sequences, int age, int nCells, double delta, double eta, double kappa, double gamma, double lambda, double mu, boolean useFreqModel, double[] stationaryDistribution){
            this.name = name;
            this.sequences = sequences;
            this.age = age;
            this.nCells = nCells;
            this.delta = delta;
            this.kappa = kappa;
            this.eta = eta;
            this.gamma = gamma;
            this.lambda = lambda;
            this.mu = mu;
            this.useFreqModel = useFreqModel;
            this.stationaryDistribution = stationaryDistribution;

            this.nStates = (int) (0.5 * (nCells + 1) * (nCells + 2));

            int Ntips = sequences.length;
            AFsequence afseq = new AFsequence(this.nStates);

            // Create the parameters
            Parameter deltaParameter = new Parameter.Default(delta);
            Parameter etaParameter = new Parameter.Default(eta);
            Parameter kappaParameter = new Parameter.Default(kappa);

            Parameter stemCellParameter = new Parameter.Default(1, nCells);
            Parameter gammaParam = new Parameter.Default(1, gamma);
            Parameter lambdaParam = new Parameter.Default(1, lambda);
            Parameter muParam = new Parameter.Default(1, mu);

            //Initiate taxonlist
            Taxa taxonlist = new Taxa();
            for (int iTip=0;iTip<Ntips;iTip++){
                taxonlist.addTaxon(new Taxon("C" + String.valueOf(iTip)));
            }

            //Create the seqlist
            List<int[]> seqlist = new ArrayList<int[]>();

            for (int iTip=0;iTip<Ntips;iTip++) {
                seqlist.add(new AFsequence(sequences[iTip]).getSequence());
            }

            // Create the patterns
            this.patterns = new Patterns(afseq, taxonlist);
            for (int site = 0; site < (seqlist.get(0)).length; site++){
                int[] current_pattern = new int[seqlist.size()];
                for (int seq_index = 0; seq_index < seqlist.size(); seq_index++){
                    current_pattern[seq_index] = seqlist.get(seq_index)[site];
                }
                this.patterns.addPattern(current_pattern);
            }

            System.out.println("Input read. Working with: ");
            System.out.println("  - Taxa count = " + patterns.getTaxonCount());
            System.out.println("  - Site count = " + patterns.getPatternCount());
            System.out.println("  - Expected states = " + patterns.getStateCount());

            //Tree simulation. We do not care about the simulation parameters, we just need a tree object with the proper number of external nodes (tips) for the super(FlipFlopErrorModel) to initalize the tip states properly.
            // We do not care about anything else
            CoalescentSimulator simulator = new CoalescentSimulator();
            DemographicModel constantPop = new ConstantPopulationModel(new Parameter.Default(age), dr.evolution.util.Units.Type.YEARS);
            Tree theTree=simulator.simulateTree(taxonlist,constantPop);

            TreeModel treeModel = new TreeModel(theTree);
            FlipFlopErrorModel errorModel= new FlipFlopErrorModel(taxonlist, new Taxa(), stemCellParameter, deltaParameter, etaParameter, kappaParameter);

            //Initializing the error model
            errorModel.setTree(theTree);
            for (int iTip=0;iTip<Ntips;iTip++) {
                errorModel.setStates(this.patterns,iTip,iTip,this.patterns.getTaxonId(iTip));
            }

            FrequencyModel freqModel = new FrequencyModel(afseq, this.stationaryDistribution);

            FlipFlopModel model = new FlipFlopModel("test", stemCellParameter, gammaParam, lambdaParam, muParam, useFreqModel, freqModel, null);

            GammaSiteModel siteModel = new GammaSiteModel(model);

            this.treeLikelihood = new TreeLikelihood(patterns, treeModel, siteModel, null, errorModel,
                    false, false, true, false, false);
            treeLikelihood.setId(TreeLikelihoodParser.TREE_LIKELIHOOD);
        }

        public TreeLikelihood getTreeLikelihoodModel(){
            return this.treeLikelihood;
        };
        public String getName(){
            return this.name;
        }

        protected final String[] sequences;
        protected final int age;
        protected final int nCells;
        protected final int nStates;
        protected final double delta;
        protected final double eta;
        protected final double kappa;
        protected final double gamma;
        protected final double lambda;
        protected final double mu;
        protected final boolean useFreqModel;
        protected final double[] stationaryDistribution;
        protected final Patterns patterns;
        protected final TreeLikelihood treeLikelihood;
        private String name;
    };


    String[] sequences = new String[] {
            "0.1, 0.4, 0.14, 0.12, 0.99, .95, 0.01",
            "0.2, 0.3, 0.13, 0.42, 0.29, .01, 0.91",
            "0.6, 0.36, 0.1, 0.62, 0.24, .15, 0.73"
    };

    int nCells=2;
    int age = 2;
    double delta=0.05;
    double eta=0.95;
    double kappa=0.9;
    double gamma=0.05;
    double lambda=0.95;
    double mu=0.05;
    double[] stationaryDistributionTest = new double[]{0.16666666666666666, 0.3666666666666667, 0.016666666666666663, 0.11666666666666665, 0.16666666666666666, 0.16666666666666666};

    Test test1 = new Test("test1", sequences, age, nCells, delta, eta, kappa, gamma, lambda, mu, true, stationaryDistributionTest);
    Test test2 = new Test("test2", sequences, age, nCells, delta, eta, kappa, gamma, lambda, mu, false, stationaryDistributionTest);

    String[] sequences2 = new String[] {
            "0.1, 0.4, 0.14, 0.12, 0.99, .95, 0.01",
            "0.2, 0.3, 0.13, 0.42, 0.29, .01, 0.91",
            "0.6, 0.36, 0.1, 0.62, 0.24, .15, 0.73",
            "0.1, 0.05, 0.2, 0.52, 0.20, .08, 0.6",
            "0.12, 0.16, 0.23, 0.42, 0.18, .11, 0.33"
    };

    Test test3 = new Test("test3", sequences2, age, nCells, delta, eta, kappa, gamma, lambda, mu, false, stationaryDistributionTest);
    Test test4 = new Test("test4", sequences2, age, nCells, delta, eta, kappa, gamma, lambda, mu, false, stationaryDistributionTest);

    Test[] all = {test1, test2, test3, test4}; //Add more instances of tests here

    public void testFlipFlops() {
        for (Test test : all) {
            TreeLikelihood treeModel = test.getTreeLikelihoodModel();
            System.out.println("The initial likelihood of " + test.getName() + " is " +  treeModel.getLogLikelihood());
/*
            OperatorSchedule schedule = new SimpleOperatorSchedule();

            Parameter kParam = new Parameter.Default(1, 0.5);
            MCMCOperator operator = new ScaleOperator(kParam, 0.5);
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
            //loggers[0].add(kParam);

            loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), 1, false);
            loggers[1].add(treeLikelihood);
            loggers[1].add(rootHeight);
            //loggers[1].add(kParam);

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
            }*/

        }
    }
}

