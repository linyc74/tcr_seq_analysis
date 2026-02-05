import os
import numpy as np
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt
from typing import List, Tuple
from .template import Processor


SEQUENCE_LOGO_DSTDIR = 'motif-sequence-logo'


class Clustering(Processor):

    count_df: pd.DataFrame
    rpm_cutoff: float
    clustering_identity: float

    rpm_df: pd.DataFrame
    mmseqs_df: pd.DataFrame
    motif_count_df: pd.DataFrame

    def main(
            self,
            count_df: pd.DataFrame,
            rpm_cutoff: float,
            clustering_identity: float) -> Tuple[pd.DataFrame, pd.DataFrame]:

        self.count_df = count_df.copy()
        self.rpm_cutoff = rpm_cutoff
        self.clustering_identity = clustering_identity

        self.rpm_df = normalize_to_rpm(self.count_df)
        self.filter_by_rpm_and_write_fasta()
        self.run_mmseqs()
        self.merge_motif_id_to_count_df()
        self.motif_count_df = self.count_df.groupby('motif_id').sum()

        self.generate_sequence_logo()
        self.zip_sequence_logo()

        return self.count_df, self.motif_count_df

    def filter_by_rpm_and_write_fasta(self):
        avg_rpm = self.rpm_df.mean(axis=1)
        tcr_sequences = self.rpm_df[avg_rpm >= self.rpm_cutoff].index
        tcr_sequences = sorted(tcr_sequences)  # to make mmseqs deterministic
        self.logger.info(f'''\
Filtering out TCRs with average RPM less than {self.rpm_cutoff}
Unique TCRs before filtering: {len(self.rpm_df)}, after filtering: {len(tcr_sequences)}''')

        with FastaWriter(file=f'{self.workdir}/tcr.fasta', mode='w') as fasta:
            for tcr in tcr_sequences:
                fasta.write(header=tcr, sequence=tcr)
        
    def run_mmseqs(self):
        lines = [
            'mmseqs easy-cluster',
            f'{self.workdir}/tcr.fasta',
            f'{self.workdir}/mmseqs',  # output prefix
            f'{self.workdir}/mmseqs-tmp',  # temporary directory
            f'--min-seq-id {self.clustering_identity}',
            f'-c {self.clustering_identity}',
            f'--cov-mode 0',
            f'--alignment-mode 3',
            f'--threads {self.threads}',
            f'1> {self.outdir}/mmseqs.log',
            f'2> {self.outdir}/mmseqs.log',
        ]
        cmd = self.CMD_LINEBREAK.join(lines)
        self.call(cmd)

        self.mmseqs_df = pd.read_csv(
            f'{self.workdir}/mmseqs_cluster.tsv',
            sep='\t',
            header=None,
            names=['representative_cdr3_amino_acid_sequence', 'cdr3_amino_acid_sequence']
        )
        self.mmseqs_df.sort_values(
            ['representative_cdr3_amino_acid_sequence', 'cdr3_amino_acid_sequence'],
            inplace=True,
            ignore_index=True
        )

        # representative_sequence -> motif_id
        unique_rep_seqs = self.mmseqs_df['representative_cdr3_amino_acid_sequence'].unique()
        d = {seq: f'motif_{count+1:06d}' for count, seq in enumerate(unique_rep_seqs)}
        self.mmseqs_df['motif_id'] = self.mmseqs_df['representative_cdr3_amino_acid_sequence'].map(d)

        self.mmseqs_df = self.mmseqs_df.drop(
            columns=['representative_cdr3_amino_acid_sequence']
        ).set_index(
            keys='cdr3_amino_acid_sequence'
        )

        self.logger.info(f'Number of motifs: {len(self.mmseqs_df["motif_id"].unique())}')

    def merge_motif_id_to_count_df(self):
        self.count_df = self.count_df.merge(
            right=self.mmseqs_df,
            left_index=True,
            right_index=True,
            how='left'
        )
        self.count_df['motif_id'] = self.count_df['motif_id'].fillna(value='low_abundance')
        columns = self.count_df.columns.tolist()
        reorder = [columns[-1]] + columns[:-1]  # move motif_id to the first column
        self.count_df = self.count_df[reorder]

    def generate_sequence_logo(self):
        dstdir = f'{self.outdir}/{SEQUENCE_LOGO_DSTDIR}'
        mafft_dir = f'{self.workdir}/mafft-tmp'
        for d in [dstdir, mafft_dir]:
            os.makedirs(d, exist_ok=True)

        for motif_id, df in self.mmseqs_df.groupby('motif_id', sort=True):
            sequences = df.index.tolist()

            faa = f'{mafft_dir}/{motif_id}.faa'
            aln = f'{mafft_dir}/{motif_id}.aln.faa'
            with FastaWriter(faa, mode='w') as fasta:
                for s in sequences:
                    fasta.write(header=s, sequence=s)

            self.run_mafft(src=faa, dst=aln)
            alignment = read_fasta(file=aln)

            logo_df = lm.alignment_to_matrix(
                alignment,
                characters_to_ignore='',  # include gaps '-', do not leave out any character
                pseudocount=0.0,
                to_type='probability'
            )

            logo_df.index = range(1, len(logo_df) + 1)  # position starts from 1
            logo_df.to_csv(f'{dstdir}/{motif_id}.csv', index=True)

            plt.rcParams.update({'font.size': 7})
            plt.rcParams.update({'lines.linewidth': 0.5})

            seq_length = len(logo_df)
            w = (1.3 + seq_length * 0.7) / 2.54
            h = 3.5 / 2.54

            color_scheme = lm.src.colors.get_color_dict(color_scheme='chemistry', chars='ACDEFGHIKLMNPQRSTVWY')
            color_scheme['-'] = np.array([0.9, 0.9, 0.9])  # set the gap character '-' to light gray

            lm.Logo(logo_df, color_scheme=color_scheme, figsize=(w, h))
            
            plt.title(motif_id, fontsize=7)
            plt.xlabel('Position', fontsize=7)
            plt.xticks(range(1, seq_length + 1))
            plt.ylabel('Frequency', fontsize=7)
            plt.tight_layout()
            plt.savefig(f'{dstdir}/{motif_id}.png', dpi=600)
            plt.close()

    def run_mafft(self, src: str, dst: str):
        """
        Args:
            src: input fasta file

            dst: output aligned fasta file
        """
        cmd = self.CMD_LINEBREAK.join([
            'mafft',
            f'--auto',
            f'--thread {self.threads}',
            f'{src}',
            f'1> {dst}',
            f'2>> {self.outdir}/mafft.log',
        ])
        self.call(cmd)

    def zip_sequence_logo(self):
        self.call(f'tar -C "{self.outdir}" -czf "{self.outdir}/{SEQUENCE_LOGO_DSTDIR}.tar.gz" {SEQUENCE_LOGO_DSTDIR}')
        self.call(f'rm -r "{self.outdir}/{SEQUENCE_LOGO_DSTDIR}"')


class FastaWriter:

    def __init__(self, file: str, mode: str = 'w'):
        """
        Args:
            file: path-like

            mode: 'w' for write or 'a' for append
        """
        assert mode in ['w', 'a']

        self.__fasta = open(file, mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __wrap(self, seq: str, length: int) -> str:
        """
        Wraps a single line of string by a specified length

        Args:
            seq: a single line of DNA or protein sequence without '\n'

            length: length of each line
        """
        if len(seq) <= length:
            return seq  # no need to wrap

        w = length
        list_ = []
        for i in range(int(len(seq)/w) + 1):
            list_.append(seq[i*w:(i+1)*w])

        return '\n'.join(list_)

    def write(self, header: str, sequence: str, wrap: int = 80):
        """
        Args:
            header: Fasta header

            sequence: DNA or protein sequence without '\n'

            wrap: length of each wrapped line for the sequence
        """
        seq = self.__wrap(sequence, wrap)
        self.__fasta.write(f'>{header}\n{seq}\n')

    def close(self):
        self.__fasta.close()


def read_fasta(file: str) -> List[str]:
    """
    Args:
        file: fasta file path
    """
    sequences: List[str] = []
    buffer: List[str] = []
    with open(file) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if buffer:
                    sequences.append(''.join(buffer))
                    buffer = []
                continue
            buffer.append(line)
    if buffer:
        sequences.append(''.join(buffer))
    return sequences


def normalize_to_rpm(df: pd.DataFrame) -> pd.DataFrame:
    sum_per_column = df.sum(axis=0)
    return df.divide(sum_per_column, axis=1) * 1000000
