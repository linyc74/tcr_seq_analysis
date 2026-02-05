import pandas as pd
from .template import Processor
from .clustering import Clustering
from .profile_plot import ProfilePlot
from .compile_table import CompileTable
from .diversity_clonality import DiversityClonality
from .differential_abundance import DifferentialAbundance


class TcrSeqAnalysis(Processor):

    count_df: pd.DataFrame
    motif_count_df: pd.DataFrame

    def main(
            self,
            csv_dir: str,
            csv_suffix: str,
            sample_sheet: str,
            clonal_index_column: str,
            count_column: str,
            group_column: str,
            rpm_cutoff: float,
            clustering_identity: float,
            p_value: float):

        self.count_df = CompileTable(self.settings).main(
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

        self.count_df, self.motif_count_df = Clustering(self.settings).main(
            count_df=self.count_df,
            rpm_cutoff=rpm_cutoff,
            clustering_identity=clustering_identity)

        DifferentialAbundance(self.settings).main(
            df=self.motif_count_df,
            sample_sheet=sample_sheet,
            group_column=group_column,
            colors=colors,
            p_value=p_value)
