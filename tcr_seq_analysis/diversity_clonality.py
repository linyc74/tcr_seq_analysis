import pandas as pd
from .template import Processor


class DiversityClonality(Processor):

    df: pd.DataFrame
    sample_sheet: str
    group_column: str

    def main(
            self,
            df: pd.DataFrame,
            sample_sheet: str,
            group_column: str):
        self.df = df.copy()
        self.sample_sheet = sample_sheet
        self.group_column = group_column
