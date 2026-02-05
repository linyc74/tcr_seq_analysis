import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from .template import Processor
from .clustering import Clustering
from .profile_plot import ProfilePlot
from .compile_table import CompileTable
from .diversity_clonality import DiversityClonality
from .differential_abundance import DifferentialAbundance


class TcrSeqAnalysis(Processor):

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
            p_value: float,
            colormap: str,
            invert_colors: bool):

        count_df = CompileTable(self.settings).main(
            csv_dir=csv_dir,
            csv_suffix=csv_suffix,
            sample_sheet=sample_sheet,
            clonal_index_column=clonal_index_column,
            count_column=count_column)

        colors = GetColors(self.settings).main(
            sample_sheet=sample_sheet,
            group_column=group_column,
            colormap=colormap,
            invert_colors=invert_colors)

        ProfilePlot(self.settings).main(
            count_df=count_df)

        DiversityClonality(self.settings).main(
            count_df=count_df,
            sample_sheet=sample_sheet,
            group_column=group_column,
            colors=colors)

        count_df, motif_count_df = Clustering(self.settings).main(
            count_df=count_df,
            rpm_cutoff=rpm_cutoff,
            clustering_identity=clustering_identity)

        DifferentialAbundance(self.settings).main(
            df=motif_count_df,
            sample_sheet=sample_sheet,
            group_column=group_column,
            colors=colors,
            p_value=p_value)

        count_df.to_csv(f'{self.outdir}/count-table.csv', index=True)
        motif_count_df.to_csv(f'{self.outdir}/motif-count-table.csv', index=True)

        self.call(f'mkdir -p {self.outdir}/log')
        self.call(f'mv {self.outdir}/*.log {self.outdir}/log/')


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
