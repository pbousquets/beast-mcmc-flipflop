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
        self.checkValues()

    def parseSamples(self) -> None:
        for column in self.sample_names:
            sample: methylationSample = methylationSample(self.input[column], self.precision)
            self.samples.append(sample.toString())

    def checkValues(self) -> None:
        assert len(self.input.select_dtypes(include="object").columns) == 0, "\nStrings found in the dataframe. Did you to forget using the --stripRownames flag?\n"
        assert all(self.input.max() <= 1) and all(self.input.min()) >= 0, "\nThe data frame contains values not in (0,1) range. \nDid you to forget using the --stripRownames flag?\n"
