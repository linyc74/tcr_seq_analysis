import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from .template import Processor


class AddGroupColumn(Processor):

    NA_VALUE: str = 'None'  # Don't use 'NA', which would make dtype = `float` but not `str`, tricky for testing

    df: pd.DataFrame
    sample_sheet: str
    group_column: str

    sample_df: pd.DataFrame

    def main(
            self,
            df: pd.DataFrame,
            sample_sheet: str,
            group_column: str) -> pd.DataFrame:

        self.df = df.copy(deep=True)
        self.sample_sheet = sample_sheet
        self.group_column = group_column

        self.read_sample_sheet()
        self.add_group_column()
        self.reorder_columns()

        return self.df

    def read_sample_sheet(self):
        self.sample_df = pd.read_csv(self.sample_sheet, index_col=0)
        assert self.group_column in self.sample_df.columns, \
            f'No "{self.group_column}" column in {self.sample_sheet}'

    def add_group_column(self):
        self.df = self.df.merge(
            right=self.sample_df[self.group_column],
            how='left',
            left_index=True,
            right_index=True)
        self.df[self.group_column] = self.df[self.group_column].fillna(self.NA_VALUE)

    def reorder_columns(self):
        cols = list(self.df.columns)
        cols = [cols[-1]] + cols[:-1]
        self.df = self.df[cols]


class GetColors(Processor):

    sample_sheet: str
    group_column: str
    colormap: str
    invert_colors: bool

    def main(
            self,
            sample_sheet: str,
            group_column: str,
            colormap: str,
            invert_colors: bool) -> list:

        self.sample_sheet = sample_sheet
        self.group_column = group_column
        self.colormap = colormap
        self.invert_colors = invert_colors

        df = pd.read_csv(self.sample_sheet, index_col=0)
        n_groups = len(df[self.group_column].unique())

        if ',' in self.colormap:
            names = self.colormap.split(',')
            if len(names) != n_groups:
                self.logger.info(f'WARNING! Number of colors "{self.colormap}" does not match number of groups ({n_groups})')
            colors = [to_rgba(n) for n in names]
        else:
            cmap = plt.colormaps[self.colormap]
            colors = [cmap(i) for i in range(n_groups)]

        if self.invert_colors:
            colors = colors[::-1]

        return colors
