import os
from .template import Settings
from .tcr_seq_analysis import TcrSeqAnalysis


def main(
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
        invert_colors: bool,
        for_publication: bool,
        threads: int,
        debug: bool,
        outdir: str):

    prefix = os.path.basename(outdir)
    for c in [' ', ',', '(', ')']:
        prefix = prefix.replace(c, '_')
    workdir = get_temp_path(prefix=f'./{prefix}_')
    
    settings = Settings(
        workdir=workdir,
        outdir=outdir,
        threads=threads,
        debug=debug,
        mock=False,
        for_publication=for_publication)

    for d in [settings.workdir, settings.outdir]:
        os.makedirs(d, exist_ok=True)

    TcrSeqAnalysis(settings).main(
        csv_dir=csv_dir,
        csv_suffix=csv_suffix,
        sample_sheet=sample_sheet,
        clonal_index_column=clonal_index_column,
        count_column=count_column,
        group_column=group_column,
        rpm_cutoff=rpm_cutoff,
        clustering_identity=clustering_identity,
        p_value=p_value,
        colormap=colormap,
        invert_colors=invert_colors)


def get_temp_path(
        prefix: str = 'temp',
        suffix: str = '') -> str:

    i = 1
    while True:
        fpath = f'{prefix}{i:03}{suffix}'
        if not os.path.exists(fpath):
            return fpath
        i += 1
