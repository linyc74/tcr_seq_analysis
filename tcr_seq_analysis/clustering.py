import pandas as pd
from .template import Processor


class Clustering(Processor):

    df: pd.DataFrame
    rpm_cutoff: float

    def main(self, df: pd.DataFrame, rpm_cutoff: float):
        self.df = df.copy()
        self.rpm_cutoff = rpm_cutoff
