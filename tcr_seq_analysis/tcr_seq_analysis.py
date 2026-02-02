import pandas as pd
from yltk import Processor


class TcrSeqAnalysis(Processor):

    df: pd.DataFrame

    def main(self):
        self.df = CompileTable(self.settings).main(
            csvdir=CSVDIR,
            sample_sheet=SAMPLE_SHEET,
            csv_suffix=CSV_SUFFIX,
            clonal_index_columns=CLONAL_INDEX_COLUMNS,
            count_column=COUNT_COLUMN)

        DiversityClonality(self.settings).main(df=self.df, sample_sheet=SAMPLE_SHEET, group_column=GROUP_COLUMN)

        Tree(self.settings).main(df=self.df, rpm_cutoff=RPM_CUTOFF)

        Clustering(self.settings).main(df=self.df, rpm_cutoff=RPM_CUTOFF)
