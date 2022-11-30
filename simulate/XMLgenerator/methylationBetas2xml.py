import click
import re
import pandas as pd
from typing import List

from yattag import Doc
from yattag import indent

"""
This script is intended to be used for the FlipFlop model. It can create the XML file needed to run the model in BEAST.
It only needs the beta values of the fCpGs in a data frame in comma-separated format. Rows represent each fCpG, columns
represent each sample
"""


class methylationSample:
    """
    Class for reading samples and casting them into a single comma-separated element
    """

    def __init__(self, methylationArray: pd.Series, precision: int):
        self.methylationArray = methylationArray.round(precision)  # Round floats, we don't need long decimals

    def toString(self) -> str:
        return self.methylationArray.astype("string").str.cat(sep=",")  # Join all sites into a single comma-separated string


class readMethylation:
    """
    Class representing the whole dataset
    """

    def __init__(self, inputStr: str, precision: int, stripRownames: bool):
        index_col = 0 if stripRownames else None
        self.input = pd.read_csv(inputStr, index_col=index_col)
        self.precision = precision
        self.samples: List[str] = []
        self.sample_names: List[str] = self.input.columns.tolist()

    def parseSamples(self) -> None:
        for column in self.sample_names:
            sample: methylationSample = methylationSample(self.input[column], self.precision)
            self.samples.append(sample.toString())


class createXML:
    """
    This class builds the actual XML file
    """

    def __init__(self) -> None:
        self.doc, self.tag, self.text, self.line = Doc().ttl()
        self.doc.asis('<?xml version="1.0" standalone="yes"?>')

    def addSamples(self, readMet: readMethylation):
        self.sequences: List[str] = readMet.samples
        self.names: List[str] = readMet.sample_names

    def buildDoc(self):
        with self.tag("beast"):
            with self.tag("taxa"):
                self.doc.attr(id="taxa")
                for name in self.names:
                    with self.tag("taxon", id=f"{name}"):
                        pass

            with self.tag("afalignment", id="alignment"):
                for i, name in enumerate(self.names):
                    with self.tag("afsequence"):
                        with self.tag("taxon", idref=f"{name}"):
                            pass
                        self.text(self.sequences[i])
                with self.tag("states"):
                    with self.tag("parameter", id="alignment.states", value=6):
                        pass

    def printDocument(self, output: str) -> None:
        text = self.__formatDoc()

        with open(output, "w") as outfile:
            outfile.write(text)

    def __formatDoc(self) -> str:
        text = indent(self.doc.getvalue(), indentation="    ", newline="\n", indent_text=True) # Pretty format
        return re.sub("></[a-zA-Z]+>", "/>", text) # transforms empty tags to single format: <tag ...></tag> --> <tag .../>

    def __str__(self) -> str:
        return self.__formatDoc()


@click.command(context_settings={"show_default": True})
@click.option("--input", required=True, help="Methylation array beta values in CSV format", type=str)
@click.option("--precision", default=3, help="Number of significant digits to consider when rounding the values", type=int)
@click.option("--stripRownames", "stripRownames", default=False, help="Whether the input contains row names or not (they need to be removed)", is_flag=True)
@click.option("--output", default="text.xml", help="Output file name", type=str)
def main(input: str, precision: int = 3, output: str = "test.xml", stripRownames: bool = True) -> None:
    myobj = readMethylation(input, precision, stripRownames)
    myobj.parseSamples()

    XMLfile = createXML()
    XMLfile.addSamples(myobj)
    XMLfile.buildDoc()
    XMLfile.printDocument(output)


if __name__ == "__main__":
    main()
