import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.axes
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import mannwhitneyu
from typing import List, Dict, Any, Tuple
from statsmodels.stats.multitest import multipletests
from .template import Processor
from .grouping import AddGroupColumn


DSTDIR_NAME = 'motif-differential-abundance'


class DifferentialAbundance(Processor):

    df: pd.DataFrame
    sample_sheet: str
    group_column: str
    colors: List[Tuple[float, float, float, float]]
    p_value: float

    def main(
            self,
            df: pd.DataFrame,
            sample_sheet: str,
            group_column: str,
            colors: List[Tuple[float, float, float, float]],
            p_value: float):

        self.df = df
        self.sample_sheet = sample_sheet
        self.group_column = group_column
        self.colors = colors
        self.p_value = p_value

        self.df = normalize_to_log_rpm_pseudocount(self.df)
        self.df = self.df.transpose()

        self.df = AddGroupColumn(self.settings).main(
            df=self.df,
            sample_sheet=self.sample_sheet,
            group_column=self.group_column)

        MannwhitneyuTestsAndBoxplots(self.settings).main(
            df=self.df,
            group_column=self.group_column,
            colors=self.colors,
            p_value=self.p_value)

        self.zip_dstdir()

    def zip_dstdir(self):
        self.call(f'tar -C "{self.outdir}" -czf "{self.outdir}/{DSTDIR_NAME}.tar.gz" {DSTDIR_NAME}')
        self.call(f'rm -r "{self.outdir}/{DSTDIR_NAME}"')


def normalize_to_log_rpm_pseudocount(df: pd.DataFrame) -> pd.DataFrame:
    sum_per_column = df.sum(axis=0)
    df = df.divide(sum_per_column, axis=1) * 1000000
    df = df + 1
    df = np.log10(df)
    return df


class MannwhitneyuTestsAndBoxplots(Processor):

    df: pd.DataFrame
    group_column: str
    colors: List[Tuple[float, float, float, float]]
    p_value: float

    groups: List[str]

    def main(
            self,
            df: pd.DataFrame,
            group_column: str,
            colors: List[Tuple[float, float, float, float]],
            p_value: float):

        self.df = df
        self.group_column = group_column
        self.colors = colors
        self.p_value = p_value

        self.plot_all()

        self.groups = self.df[self.group_column].unique().tolist()

        for group_1, group_2 in combinations(self.groups, 2):
            self.process_group_pair(group_1=group_1, group_2=group_2)

    def plot_all(self):
        for motif_id in self.df.columns:
            if motif_id == self.group_column:
                continue
            dstdir = f'{self.outdir}/{DSTDIR_NAME}/all'
            os.makedirs(dstdir, exist_ok=True)
            Boxplot(self.settings).main(
                data=self.df,
                x=self.group_column,
                y=motif_id,
                colors=self.colors,
                title=motif_id,
                png=f'{dstdir}/{motif_id}.png'
            )

    def process_group_pair(self, group_1: str, group_2: str):
        dstdir = f'{self.outdir}/{DSTDIR_NAME}/{group_1}-{group_2}'
        os.makedirs(dstdir, exist_ok=True)

        color_1 = self.colors[self.groups.index(group_1)]
        color_2 = self.colors[self.groups.index(group_2)]
        stats_data = []

        for motif_id in self.df.columns:

            if motif_id == self.group_column:
                continue

            is_group_1 = self.df[self.group_column] == group_1
            is_group_2 = self.df[self.group_column] == group_2

            statistic, pvalue = mannwhitneyu(
                x=self.df.loc[is_group_1, motif_id],
                y=self.df.loc[is_group_2, motif_id]
            )

            if pvalue <= self.p_value:
                Boxplot(self.settings).main(
                    data=self.df[is_group_1 | is_group_2],
                    x=self.group_column,
                    y=motif_id,
                    colors=[color_1, color_2],
                    title=f'{motif_id}\n$p = {pvalue:.4f}$',
                    png=f'{dstdir}/{pvalue:.4f}_{motif_id}.png'
                )

            stats_data.append({
                'Motif ID': motif_id,
                'Mean 1 (%)': self.df.loc[is_group_1, motif_id].mean(),
                'Mean 2 (%)': self.df.loc[is_group_2, motif_id].mean(),
                'Statistics': statistic,
                'P value': pvalue,
            })

        self.__save_stats_data(stats_data=stats_data, dstdir=dstdir)

    def __save_stats_data(self, stats_data: List[Dict[str, Any]], dstdir: str):
        stats_df = pd.DataFrame(stats_data).sort_values(
            by='P value',
            ascending=True
        )
        rejected, pvals_corrected, _, _ = multipletests(
            stats_df['P value'],
            alpha=0.1,
            method='fdr_bh',  # Benjamini-Hochberg
            is_sorted=False,
            returnsorted=False)
        stats_df['Benjamini-Hochberg adjusted P value'] = pvals_corrected
        stats_df.to_csv(
            f'{dstdir}/Mann-Whitney-U.csv',
            index=False
        )


class Boxplot(Processor):

    WIDTH_PADDING = 1.5 / 2.54
    WIDTH_PER_GROUP = 1.25 / 2.54
    HEIGHT = 5 / 2.54
    DPI = 600
    FONT_SIZE = 6
    BOX_WIDTH = 0.5
    XLABEL = None
    YLABEL = 'Log10(RPM + 1)'
    LINEWIDTH = 0.5
    BOX_lINEWIDTH = 0.5
    MARKER_LINEWIDTH = 0.25
    YLIM = None

    data: pd.DataFrame
    x: str
    y: str
    colors: List[Tuple[float, float, float, float]]
    title: str
    png: str

    ax: matplotlib.axes.Axes

    def main(
            self,
            data: pd.DataFrame,
            x: str,
            y: str,
            colors: List[Tuple[float, float, float, float]],
            title: str,
            png: str):

        self.data = data
        self.x = x
        self.y = y
        self.colors = colors
        self.title = title
        self.png = png

        self.init()
        self.plot()
        self.config()
        self.save()

    def init(self):
        plt.rcParams['font.size'] = self.FONT_SIZE
        plt.rcParams['axes.linewidth'] = self.LINEWIDTH

        groups = len(self.data[self.x].unique())
        figsize = (groups * self.WIDTH_PER_GROUP + self.WIDTH_PADDING, self.HEIGHT)

        plt.figure(figsize=figsize)

    def plot(self):
        self.ax = sns.boxplot(
            data=self.data,
            x=self.x,
            y=self.y,
            hue=self.x,
            palette=self.colors,
            width=self.BOX_WIDTH,
            linewidth=self.BOX_lINEWIDTH,
            dodge=False,  # to align the boxes on the x axis
        )
        self.ax = sns.stripplot(
            data=self.data,
            x=self.x,
            y=self.y,
            hue=self.x,
            palette=self.colors,
            linewidth=self.MARKER_LINEWIDTH,
        )

    def config(self):
        self.ax.set_title(self.title)
        self.ax.set(xlabel=self.XLABEL, ylabel=self.YLABEL)
        plt.gca().xaxis.set_tick_params(width=self.LINEWIDTH)
        plt.gca().yaxis.set_tick_params(width=self.LINEWIDTH)
        plt.ylim(self.YLIM)
        plt.legend().remove()

    def save(self):
        plt.tight_layout()
        plt.savefig(self.png, dpi=self.DPI)
        plt.close()
