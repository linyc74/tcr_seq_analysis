import pandas as pd
from .template import Processor
from .clustering import Clustering
from .profile_plot import ProfilePlot
from .compile_table import CompileTable
from .diversity_clonality import DiversityClonality


class TcrSeqAnalysis(Processor):

    count_df: pd.DataFrame
    rpm_df: pd.DataFrame
    motif_count_df: pd.DataFrame
    motif_rpm_df: pd.DataFrame

    def main(
            self,
            csv_dir: str,
            csv_suffix: str,
            sample_sheet: str,
            clonal_index_column: str,
            count_column: str,
            group_column: str,
            rpm_cutoff: float,
            clustering_identity: float):

        self.count_df, self.rpm_df = CompileTable(self.settings).main(
            csv_dir=csv_dir,
            csv_suffix=csv_suffix,
            sample_sheet=sample_sheet,
            clonal_index_column=clonal_index_column,
            count_column=count_column)

        ProfilePlot(self.settings).main(
            count_df=self.count_df)

        DiversityClonality(self.settings).main(
            count_df=self.count_df,
            sample_sheet=sample_sheet,
            group_column=group_column)

        self.count_df, self.rpm_df, self.motif_count_df, self.motif_rpm_df = Clustering(self.settings).main(
            count_df=self.count_df,
            rpm_df=self.rpm_df,
            rpm_cutoff=rpm_cutoff,
            clustering_identity=clustering_identity)
