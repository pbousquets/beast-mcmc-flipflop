import re
from typing import List

from yattag import Doc
from yattag import indent

from readInputMethylation import readMethylation


class createXML:
    """
    This class builds the actual XML file
    """

    def __init__(self, age:int, stemCells: int = 3, delta: float = 0.2, eta: float = 0.7, kappa: float = 50, mu: float = 0.1, gamma: float = 0.1, Lambda: float = 1, normalize = False) -> None:
        self.age = age
        self.stemCells = stemCells
        self.delta = delta
        self.eta = eta
        self.kappa = kappa
        self.mu = mu
        self.gamma = gamma
        self.Lambda = Lambda
        self.normalize = "true" if normalize else "false"
        self.doc, self.tag, self.text, self.line = Doc().ttl()
        self.doc.asis('<?xml version="1.0" standalone="yes"?>')

    def addSamples(self, readMet: readMethylation):
        self.sequences: List[str] = readMet.samples
        self.names: List[str] = readMet.sample_names

    def buildDoc(self, output: str, iterations: int = 750000, sampling: int = 75, screenSampling: int = 75):
        with self.tag("beast"):

            self.newSection("The list of taxa to be analysed (can also include dates/ages).")
            with self.tag("taxa"):
                self.doc.attr(id="taxa")
                for name in self.names:
                    with self.tag("taxon", id=f"{name}"):
                        self.doc.stag("date", value=f"{self.age}", direction="forwards", units="years")

            self.newSection("The list of sequences to be analysed")
            with self.tag("afalignment", id="alignment"):
                for i, name in enumerate(self.names):
                    with self.tag("afsequence"):
                        self.doc.stag("taxon", idref=f"{name}")
                        self.text(self.sequences[i])
                with self.tag("stemCells"):
                    self.doc.stag("parameter", id="alignment.stemCells", value=self.stemCells)

            self.newSection("Initialize the error model")
            with self.tag("AFsequenceErrorModel", id="errorModel"):
                with self.tag("stemCells"):
                    self.doc.stag("parameter", idref="alignment.stemCells")
                with self.tag("deltaOffset"):
                    self.doc.stag("parameter", id="errorModel.deltaOffset", value=self.delta, lower="0.0", upper="1.0")
                with self.tag("etaOffset"):
                    self.doc.stag("parameter", id="errorModel.etaOffset", value=self.eta, lower="0.0", upper="1.0")
                with self.tag("kappaScale"):
                    self.doc.stag("parameter", id="errorModel.kappaScale", value=self.kappa, lower="0.0")

            self.newSection("A prior assumption that the population size has remained constant")
            self.newSection("throughout the time spanned by the genealogy", addNewline=False)
            with self.tag("constantSize", id="constant", units="years"):
                with self.tag("populationSize"):
                    self.doc.stag("parameter", id="constant.popSize", value="1", lower="0.0")

            self.newSection("Generate a random starting tree under the coalescent process")
            with self.tag("coalescentSimulator", id="startingTree"):
                self.doc.stag("taxa", idref="taxa")
                self.doc.stag("constantSize", idref="constant")

            self.newSection("Generate a tree model")
            with self.tag("treeModel", id="treeModel"):
                self.doc.stag("coalescentTree", idref="startingTree")
                with self.tag("rootHeight"):
                    self.doc.stag("parameter", id="treeModel.rootHeight")
                with self.tag("nodeHeights", internalNodes="true"):
                    self.doc.stag("parameter", id="treeModel.internalNodeHeights")
                with self.tag("nodeHeights", internalNodes="true", rootNode="true"):
                    self.doc.stag("parameter", id="treeModel.allInternalNodeHeights")

            self.newSection("Generate a coalescent likelihood")
            with self.tag("coalescentLikelihood", id="coalescent"):
                with self.tag("model"):
                    self.doc.stag("constantSize", idref="constant")
                with self.tag("populationTree"):
                    self.doc.stag("treeModel ", idref="treeModel")

            self.newSection("The strict clock (uniform rates across branches)")
            with self.tag("strictClockCenancestorBranchRates", id="branchRates"):
                with self.tag("rate"):
                    self.doc.stag("parameter", id="clock.rate", value="1")

            self.newSection("The substitution model")
            with self.tag("flipflopModel", id="flipflopSubstitutionModel", normalize = self.normalize):
                self.doc.stag("alignment", idref="alignment")
                with self.tag("stemCells"):
                    self.doc.stag("parameter", idref="alignment.stemCells")
                with self.tag("mu"):
                    self.doc.stag("parameter", id="flipflop.mu", value=self.mu, lower="0.0")
                with self.tag("gamma"):
                    self.doc.stag("parameter", id="flipflop.gamma", value=self.gamma, lower="0.0")
                with self.tag("lambda"):
                    self.doc.stag("parameter", id="flipflop.lambda", value=self.Lambda, lower="0.0")

            self.newSection("The site model")
            with self.tag("siteModel", id="siteModel"):
                with self.tag("substitutionModel"):
                    self.doc.stag("flipflopModel", idref="flipflopSubstitutionModel")

            self.newSection("The cenancestor frequency")
            with self.tag("flipflopCenancestorFrequency", id="cenancestorFrequencyModel", methylatedProportion="0.5"):
                self.doc.stag("alignment", idref="alignment")
                with self.tag("frequencies"):
                    self.doc.stag("parameter", id="cenancestor.frequencies", value="1")

            self.newSection("The cenancestor treelikelihood")
            with self.tag("cenancestorTreeLikelihood", id="treeLikelihood", useAmbiguities="false", heightRules="true"):
                self.doc.stag("alignment", idref="alignment")
                self.doc.stag("treeModel", idref="treeModel")
                self.doc.stag("siteModel", idref="siteModel")
                self.doc.stag("cenancestorFrequency", idref="cenancestorFrequencyModel")
                self.doc.stag("tipStatesModel ", idref="errorModel")
                with self.tag("cenancestorHeight"):
                    self.doc.stag("parameter", id="luca_height", value=f"{self.age}")
                with self.tag("cenancestorBranch"):
                    self.doc.stag("parameter", id="luca_branch", value="1", upper=f"{self.age}", lower="0.0")
                self.doc.stag("strictClockCenancestorBranchRates", idref="branchRates")

            self.newSection("Set the operators")
            with self.tag("operators", id="operators", optimizationSchedule="default"):

                ## error model
                self.newSection("Error model", addNewline=False)
                with self.tag("scaleOperator", scaleFactor="0.9", weight="0.75"):
                    self.doc.stag("parameter", idref="errorModel.deltaOffset")
                with self.tag("scaleOperator", scaleFactor="0.9", weight="0.75"):
                    self.doc.stag("parameter", idref="errorModel.etaOffset")
                with self.tag("scaleOperator", scaleFactor="0.9", weight="0.75"):
                    self.doc.stag("parameter", idref="errorModel.kappaScale")

                ## flipflop model
                self.newSection("flipflop model", addNewline=False)
                with self.tag("scaleOperator", scaleFactor="0.9", weight="0.5"):
                    self.doc.stag("parameter", idref="flipflop.mu")
                with self.tag("scaleOperator", scaleFactor="0.9", weight="0.5"):
                    self.doc.stag("parameter", idref="flipflop.gamma")
                with self.tag("scaleOperator", scaleFactor="0.9", weight="0.5"):
                    self.doc.stag("parameter", idref="flipflop.lambda")

                ## Demography
                self.newSection("Demography", addNewline=False)
                with self.tag("scaleOperator", scaleFactor="0.2", weight="3.0"):
                    self.doc.stag("parameter", idref="constant.popSize")

                ## Tree
                self.newSection("Tree", addNewline=False)
                with self.tag("subtreeSlide", size="2.5", gaussian="true", weight="15.0"):
                    self.newSection("2.5 years. They will be automatically optimized by BEAST though", addNewline=False)
                    self.doc.stag("parameter", idref="treeModel")
                with self.tag("narrowExchange", weight="15.0"):
                    self.doc.stag("parameter", idref="treeModel")
                with self.tag("wideExchange", weight="3.0"):
                    self.doc.stag("parameter", idref="treeModel")
                with self.tag("wilsonBalding", weight="3.0"):
                    self.doc.stag("parameter", idref="treeModel")
                with self.tag("scaleOperator", scaleFactor="0.25", weight="3.0"):
                    self.doc.stag("parameter", idref="treeModel.rootHeight")
                with self.tag("uniformOperator", weight="30.0"):
                    self.doc.stag("parameter", idref="treeModel.internalNodeHeights")
                with self.tag("upDownOperator", scaleFactor="0.25", weight="3.0"):
                    with self.tag("up"):
                        self.doc.stag("parameter", idref="flipflop.mu")
                        self.doc.stag("parameter", idref="flipflop.gamma")
                        self.doc.stag("parameter", idref="flipflop.lambda")
                    with self.tag("down"):
                        self.doc.stag("parameter", idref="treeModel.allInternalNodeHeights")
                        #self.doc.stag("parameter", idref="luca_height")

            self.newSection("Define MCMC")
            with self.tag("mcmc", id="mcmc", chainLength=f"{iterations}", autoOptimize="true", operatorAnalysis=f"{output}.ops"):
                with self.tag("posterior", id="posterior"):
                    with self.tag("prior", id="prior"):
                        self.newSection("Error model", addNewline=False)
                        with self.tag("betaPrior", shape="95.0", shapeB="5.0"):
                            self.doc.stag("parameter", idref="errorModel.etaOffset")
                        with self.tag("betaPrior", shape="5.0", shapeB="95.0"):
                            self.doc.stag("parameter", idref="errorModel.deltaOffset")
                        with self.tag("logNormalPrior", mean="4.56", stdev="0.3"):
                            self.newSection("Rounded from calculating the logscale mean and stdev for mean=100 and stdev=30", addNewline=False)
                            self.doc.stag("parameter", idref="errorModel.kappaScale")

                        self.newSection("Substitution model", addNewline=False)
                        with self.tag("halfNormalPrior", mean="0.0", stdev="0.05"):
                            self.doc.stag("parameter", idref="flipflop.mu")
                        with self.tag("halfNormalPrior", mean="0.0", stdev="0.05"):
                            self.doc.stag("parameter", idref="flipflop.gamma")
                        with self.tag("halfNormalPrior", mean="0.0", stdev="1"):
                            self.doc.stag("parameter", idref="flipflop.lambda")

                        self.newSection("Demography", addNewline=False)
                        with self.tag("oneOnXPrior"):
                            self.doc.stag("parameter", idref="constant.popSize")

                        self.newSection("Tree", addNewline=False)
                        self.doc.stag("coalescentLikelihood", idref="coalescent")

                        self.newSection("Cenancestor Prior on the height, since it is easier to have a meaningful prior on it (time of the initial development of the BE fragment)", addNewline=False)
                        with self.tag("uniformPrior", lower="0.0", upper=f"{self.age}"):
                            self.doc.stag("parameter", idref="luca_branch")

                    with self.tag("likelihood", id="likelihood"):
                        self.doc.stag("cenancestorTreeLikelihood", idref="treeLikelihood")

                self.doc.stag("operators", idref="operators")

                self.newSection("write log to screen")
                with self.tag("log", id="screenLog", logEvery=f"{screenSampling}"):
                    with self.tag("column", label="Posterior", dp="4", width="12"):
                        self.doc.stag("posterior", idref="posterior")
                    with self.tag("column", label="Prior", dp="4", width="12"):
                        self.doc.stag("prior", idref="prior")
                    with self.tag("column", label="Likelihood", dp="4", width="12"):
                        self.doc.stag("likelihood", idref="likelihood")

                self.newSection("write log to file")
                with self.tag("log", id="filelog", logEvery=f"{sampling}", fileName=f"{output}.log", overwrite="false"):
                    self.doc.stag("posterior", idref="posterior")
                    self.doc.stag("prior", idref="prior")
                    self.doc.stag("likelihood", idref="likelihood")
                    self.doc.stag("parameter", idref="errorModel.deltaOffset")
                    self.doc.stag("parameter", idref="errorModel.etaOffset")
                    self.doc.stag("parameter", idref="errorModel.kappaScale")
                    self.doc.stag("parameter", idref="flipflop.mu")
                    self.doc.stag("parameter", idref="flipflop.gamma")
                    self.doc.stag("parameter", idref="flipflop.lambda")
                    self.doc.stag("parameter", idref="alignment.stemCells")
                    self.doc.stag("parameter", idref="treeModel.rootHeight")
                    self.doc.stag("parameter", idref="luca_height")
                    self.doc.stag("parameter", idref="luca_branch")
                    self.doc.stag("parameter", idref="constant.popSize")
                    self.doc.stag("parameter", idref="clock.rate")
                    self.doc.stag("coalescentLikelihood", idref="coalescent")

                self.newSection("write tree log to file")
                with self.tag("logTree", id="treeFileLog", logEvery=f"{sampling}", nexusFormat="true", fileName=f"{output}.trees", sortTranslationTable="true"):
                    self.doc.stag("treeModel", idref="treeModel")

            with self.tag("report"):
                with self.tag("property", name="timer"):
                    self.doc.stag("mcmc", idref="mcmc")

    def newSection(self, text: str, addNewline: str = True) -> None:
        init = "\n" if addNewline else ""
        self.doc.asis(f"{init}<!-- {text} -->")

    def printDocument(self, output: str) -> None:
        text = self.__formatDoc()
        with open(f"{output}.xml", "w") as outfile:
            outfile.write(text)

    def __formatDoc(self) -> str:
        text = indent(self.doc.getvalue(), indentation="    ", newline="\n", indent_text=True, blank_is_text=True)  # Pretty format
        text = re.sub("\n\n", "\n", text)  # Remove unwanted extra spaces that the package includes
        return re.sub(" />", "/>", text)  # Remove final space before closure in single tag

    def __str__(self) -> str:
        return self.__formatDoc()
