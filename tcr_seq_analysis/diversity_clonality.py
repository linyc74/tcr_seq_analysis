import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu
from typing import List, Tuple, Dict, Union
from .template import Processor
from .grouping import AddGroupColumn


DSTDIR_NAME = 'diversity-clonality'


class DiversityClonality(Processor):

    count_df: pd.DataFrame
    sample_sheet: str
    group_column: str
    colors: list

    df: pd.DataFrame

    def main(
            self,
            count_df: pd.DataFrame,
            sample_sheet: str,
            group_column: str,
            colors: list):

        self.count_df = count_df.copy()
        self.sample_sheet = sample_sheet
        self.group_column = group_column
        self.colors = colors

        os.makedirs(f'{self.outdir}/{DSTDIR_NAME}', exist_ok=True)

        self.df = pd.DataFrame()
        sample_ids = pd.read_csv(self.sample_sheet, index_col=0).index
        for sample_id in sample_ids:
            values = self.count_df[sample_id]
            self.df.loc[sample_id, 'Observed clones'] = observed_clones(values)
            self.df.loc[sample_id, 'Shannon entropy'] = shannon_entropy(values)
            self.df.loc[sample_id, 'Clonality'] = 1 - normalized_shannon(values)
            self.df.loc[sample_id, 'Simpson concentration'] = simpson_concentration(values)
            self.df.loc[sample_id, 'Gini coefficient'] = gini_coefficient(values)

        self.df = AddGroupColumn(self.settings).main(
            df=self.df,
            sample_sheet=self.sample_sheet,
            group_column=self.group_column)

        self.df.to_csv(f'{self.outdir}/{DSTDIR_NAME}/diversity-clonality.csv', index=True)

        Plot(self.settings).main(df=self.df, colors=self.colors, group_column=self.group_column)

        MannWhitneyU(self.settings).main(df=self.df, group_column=self.group_column)


def gini_coefficient(values: np.ndarray) -> float:
    values = values[values != 0]  # remove zeros
    values = np.sort(values)  # low to high
    cumulative_sum = np.cumsum(values)
    lorenz_curve = cumulative_sum / values.sum()  # normalize by total sum
    B = np.sum(lorenz_curve) / len(values)  # the "B" area under the Lorenz curve
    return 1 - 2 * B


def observed_clones(values: np.ndarray) -> int:
    return int(np.count_nonzero(values))


def shannon_entropy(values: np.ndarray) -> float:
    total = float(values.sum())
    if total == 0.0:
        return 0.0
    proportions = values[values > 0] / total
    return float(-np.sum(proportions * np.log(proportions)))


def normalized_shannon(values: np.ndarray) -> float:
    total = float(values.sum())
    if total == 0.0:
        return 0.0
    clones = observed_clones(values)
    if clones <= 1:
        return 0.0
    entropy = shannon_entropy(values)
    return float(entropy / np.log(clones))


def simpson_concentration(values: np.ndarray) -> float:
    total = float(values.sum())
    if total == 0.0:
        return 0.0
    proportions = values[values > 0] / total
    return float(np.sum(proportions ** 2))


class Plot(Processor):

    df: pd.DataFrame
    colors: list
    group_column: str

    figsize: Tuple[float, float]
    box_width: float
    dpi: int
    line_width: float
    fontsize: int

    metrics: List[str]

    def main(self, df: pd.DataFrame, colors: list, group_column: str):
        self.df = df
        self.colors = colors
        self.group_column = group_column

        self.set_parameters()
        self.set_metrics()
        for metric in self.metrics:
            self.plot_one(metric=metric)

    def set_parameters(self):
        n_groups = self.df[self.group_column].nunique()
        if self.settings.for_publication:
            self.figsize = get_figsize_for_publication(n_groups=n_groups)
            self.box_width = 0.35
            self.dpi = 600
            self.line_width = 0.5
            self.fontsize = 7
        else:
            self.figsize = get_figsize(n_groups=n_groups)
            self.box_width = 0.5
            self.dpi = 300
            self.line_width = 1.0
            self.fontsize = 12

    def set_metrics(self):
        self.metrics = [
            c for c in self.df.columns if c != self.group_column
        ]

    def plot_one(self, metric: str):

        plt.rcParams['font.size'] = self.fontsize
        plt.rcParams['axes.linewidth'] = self.line_width

        plt.figure(figsize=self.figsize, dpi=self.dpi)

        sns.boxplot(
            data=self.df,
            x=self.group_column,
            hue=self.group_column,
            y=metric,
            palette=self.colors,
            width=self.box_width,
            linewidth=self.line_width,
            legend=False
        )
        sns.stripplot(
            data=self.df,
            x=self.group_column,
            hue=self.group_column,
            y=metric,
            palette=self.colors,
            linewidth=self.line_width,
            legend=False
        )

        plt.ylabel(metric)

        plt.tight_layout()

        f = metric.lower().replace(' ', '-')
        plt.savefig(f'{self.outdir}/{DSTDIR_NAME}/{f}.png', dpi=self.dpi)
        plt.close()


def get_figsize(n_groups: int) -> Tuple[float, float]:
    w = (3 * n_groups + 3) / 2.54
    h = 10 / 2.54
    return w, h


def get_figsize_for_publication(n_groups: int) -> Tuple[float, float]:
    w = (1.2 * n_groups + 2) / 2.54
    h = 4 / 2.54
    return w, h


class MannWhitneyU(Processor):

    df: pd.DataFrame
    group_column: str

    stats_data: List[Dict[str, Union[str, float]]]

    def main(self, df: pd.DataFrame, group_column: str):
        self.df = df
        self.group_column = group_column

        metrics = [c for c in self.df.columns if c != self.group_column]
        groups = self.df[self.group_column].unique()

        self.stats_data = []

        for metric in metrics:
            for group1, group2 in combinations(groups, 2):
                self.mann_whitney_u(metric=metric, group1=group1, group2=group2)

        pd.DataFrame(self.stats_data).to_csv(
            f'{self.outdir}/{DSTDIR_NAME}/diversity-clonality-mann-whitney-u.csv',
            index=False)

    def mann_whitney_u(self, metric: str, group1: str, group2: str):
        result = mannwhitneyu(
            x=self.df.loc[self.df[self.group_column] == group1, metric],
            y=self.df.loc[self.df[self.group_column] == group2, metric],
        )
        self.stats_data.append({
            'Diversity': metric,
            'Group 1': group1,
            'Group 2': group2,
            'P value': result.pvalue
        })
