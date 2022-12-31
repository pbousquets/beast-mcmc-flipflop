import click
from createXML import createXML
from readInputMethylation import readMethylation
from collections import OrderedDict

"""
This script is intended to be used for the FlipFlop model. It can create the XML file needed to run the model in BEAST.
It only needs the beta values of the fCpGs in a data frame in comma-separated format. Rows represent each fCpG, columns
represent each sample
"""

class CustomOption(click.Option):
    def __init__(self, *args, **kwargs):
        kwargs.update({'show_default': True})
        self.help_group = kwargs.pop('help_group', None)
        super(CustomOption, self).__init__(*args, **kwargs)


class CustomCommand(click.Command):
    def format_options(self, ctx, formatter):
        """Writes all the options into the formatter if they exist."""
        opts = OrderedDict()
        for param in self.get_params(ctx):
            rv = param.get_help_record(ctx)
            if rv is not None:
                if hasattr(param, 'help_group') and param.help_group:
                    opts.setdefault(str(param.help_group), []).append(rv)
                else:
                    opts.setdefault('Options', []).append(rv)

        for name, opts_group in opts.items():
            with formatter.section(name):
                formatter.write_dl(opts_group)


@click.command(cls=CustomCommand)
@click.option("--age", cls=CustomOption, required=True, help="Samples age in years", type=int, help_group="Inputs")
@click.option("--input", cls=CustomOption, required=True, help="Methylation array beta values in CSV format", type=str, help_group="Inputs")
@click.option("--stemCells", "stemCells", cls=CustomOption, default = 3, help="Prior for the stemCells parameter", type=int, help_group="Priors")
@click.option("--delta", cls=CustomOption, default = 0.2, help="Prior for the delta parameter", type=float, help_group="Priors")
@click.option("--eta", cls=CustomOption, default = 0.7, help="Prior for the eta parameter", type=float, help_group="Priors")
@click.option("--kappa", cls=CustomOption, default = 50, help="Prior for the kappa parameter", type=float, help_group="Priors")
@click.option("--mu", cls=CustomOption, default = 0.1, help="Prior for the delta parameter", type=float, help_group="Priors")
@click.option("--gamma", cls=CustomOption, default = 0.1, help="Prior for the eta parameter", type=float, help_group="Priors")
@click.option("--lambda", "Lambda", cls=CustomOption, default = 1, help="Prior for the kappa parameter", type=float, help_group="Priors")
@click.option("--iterations", cls=CustomOption, default=750_000, help="Number of MCMC iterations (chain length)", type=int, help_group="Extra parameters")
@click.option("--precision", cls=CustomOption, default=3, help="Number of significant digits to consider when rounding the values", type=int, help_group="Extra parameters")
@click.option("--sampling", cls=CustomOption, default=75, help="Frequency of sampling in the log file", type=int, help_group="Extra parameters")
@click.option("--screenSampling", "screenSampling", cls=CustomOption, default=75, help="Frequency of sampling to print in the screen", type=int, help_group="Extra parameters")
@click.option("--stripRownames", "stripRownames", cls=CustomOption, default=False, help="Whether the input contains row names or not (they need to be removed)", is_flag=True, help_group="Extra parameters")
@click.option("--output", cls=CustomOption, default="test", help="Output prefix for the analysis files (XML, trees, logs...)", type=str, help_group="Extra parameters")
def main(input: str, age: int, stemCells: int, delta:float, eta: float, kappa: float, mu: float, gamma: float, Lambda: float, iterations: int = 20_000, precision: int = 3, sampling: int = 200, screenSampling: int = 200, output: str = "test", stripRownames: bool = True) -> None:
    myobj = readMethylation(input, precision, stripRownames)
    myobj.parseSamples()

    XMLfile = createXML(age=age, stemCells=stemCells, delta=delta, eta=eta, kappa=kappa, mu=mu, gamma=gamma, Lambda=Lambda)
    XMLfile.addSamples(myobj)
    XMLfile.buildDoc(output, iterations=iterations, sampling=sampling, screenSampling=screenSampling)
    XMLfile.printDocument(output)


if __name__ == "__main__":
    main()
