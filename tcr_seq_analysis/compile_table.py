import os
import pandas as pd
from os.path import join
from typing import List
from .template import Processor


class CompileTable(Processor):

    csv_dir: str
    csv_suffix: str
    sample_sheet: str
    clonal_index_columns: List[str]
    count_column: str

    count_df: pd.DataFrame
    rpm_df: pd.DataFrame  # reads per million

    def main(
            self,
            csv_dir: str,
            csv_suffix: str,
            sample_sheet: str,
            clonal_index_columns: List[str],
            count_column: str) -> pd.DataFrame:

        self.csv_dir = csv_dir
        self.csv_suffix = csv_suffix
        self.sample_sheet = sample_sheet
        self.clonal_index_columns = clonal_index_columns
        self.count_column = count_column

        self.read_csvs_and_merge()
        self.count_df.fillna(value=0, inplace=True)
        self.sort_cdr3_by_sum()
        self.normalize_to_rpm()
        self.count_df.to_csv(f'{self.settings.outdir}/count-table.csv', index=True)
        self.rpm_df.to_csv(f'{self.settings.outdir}/rpm-table.csv', index=True)

        return self.rpm_df

    def read_csvs_and_merge(self):
        self.count_df = pd.DataFrame(columns=self.clonal_index_columns)
        sample_ids = pd.read_csv(self.sample_sheet, index_col=0).index
        for sample_id in sample_ids:
            csv = f'{self.csv_dir}/{sample_id}{self.csv_suffix}'
            df: pd.DataFrame = pd.read_csv(
                csv,
                usecols=self.clonal_index_columns + [self.count_column]
            ).dropna(  # index columns should not contain any na
                subset=self.clonal_index_columns,
                how='any'
            ).rename(columns={
                self.count_column: sample_id
            })

            # clonal index columns need to be unique: groupby and sum identical TCR sequences
            df = df.groupby(by=self.clonal_index_columns).sum().reset_index(drop=False)
            assert_unique(df, self.clonal_index_columns)

            self.count_df = self.count_df.merge(
                right=df,
                left_on=self.clonal_index_columns,
                right_on=self.clonal_index_columns,
                how='outer',
            )

        self.count_df = self.count_df.set_index(keys=self.clonal_index_columns)

    def sort_cdr3_by_sum(self):
        self.count_df['sum'] = self.count_df.sum(axis=1)
        self.count_df = self.count_df.sort_values(
            by='sum',
            ascending=False
        ).drop(
            columns=['sum']
        )

    def normalize_to_rpm(self):
        sum_per_column = self.count_df.sum(axis=0)
        self.rpm_df = self.count_df.divide(sum_per_column, axis=1) * 1000000


def get_files(
        source: str = '.',
        startswith: str = '',
        endswith: str = '',
        isfullpath: bool = False) -> List[str]:

    ret = []
    s, e = startswith, endswith
    for path, dirs, files in os.walk(source):
        if path == source:
            ret = [f for f in files if (f.startswith(s) and f.endswith(e))]

    if isfullpath:
        ret = [join(source, f) for f in ret]

    if ret:
        ret.sort()  # make the order consistent across OS platforms
    return ret


def assert_unique(df: pd.DataFrame, columns: List[str]):
    assert len(df) == len(df.drop_duplicates(subset=columns))
