import pandas as pd
from typing import List


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
