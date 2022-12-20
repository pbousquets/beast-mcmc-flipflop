import click
from createXML import createXML
from readInputMethylation import readMethylation

"""
This script is intended to be used for the FlipFlop model. It can create the XML file needed to run the model in BEAST.
It only needs the beta values of the fCpGs in a data frame in comma-separated format. Rows represent each fCpG, columns
represent each sample
"""


@click.command(context_settings={"show_default": True})
@click.option("--input", required=True, help="Methylation array beta values in CSV format", type=str)
@click.option("--precision", default=3, help="Number of significant digits to consider when rounding the values", type=int)
@click.option("--iterations", default=20_000, help="Number of MCMC iterations (chain length)", type=int)
@click.option("--sampling", default=200, help="Frequency of sampling in the log file", type=int)
@click.option("--screenSampling", "screenSampling", default=200, help="Frequency of sampling to print in the screen", type=int)
@click.option("--stripRownames", "stripRownames", default=False, help="Whether the input contains row names or not (they need to be removed)", is_flag=True)
@click.option("--output", default="test", help="Output prefix for the analysis files (XML, trees, logs...)", type=str)
def main(input: str, precision: int = 3, iterations: int = 20_000, sampling: int = 200, screenSampling: int = 200, output: str = "test", stripRownames: bool = True) -> None:
    myobj = readMethylation(input, precision, stripRownames)
    myobj.parseSamples()

    XMLfile = createXML()
    XMLfile.addSamples(myobj)
    XMLfile.buildDoc(output, iterations=iterations, sampling=sampling, screenSampling=screenSampling)
    XMLfile.printDocument(output)


if __name__ == "__main__":
    main()
