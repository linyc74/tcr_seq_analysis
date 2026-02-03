import os
import numpy as np
import pandas as pd
import logomaker as lm
import matplotlib.pyplot as plt
from typing import List
from .template import Processor


class Clustering(Processor):

    df: pd.DataFrame
    abundance_cutoff: float
    clustering_identity: float
    mmseqs_df: pd.DataFrame  # motif_id, cdr3_amino_acid_sequence

    def main(self, df: pd.DataFrame, abundance_cutoff: float, clustering_identity: float):
        self.df = df.copy()
        self.abundance_cutoff = abundance_cutoff
        self.clustering_identity = clustering_identity

        self.filter_by_abundance()
        self.write_fasta()
        self.run_mmseqs()
        self.generate_sequence_logo()

    def filter_by_abundance(self):
        n_before = len(self.df)
        avg_abundance = self.df.mean(axis=1)
        self.df = self.df[avg_abundance >= self.abundance_cutoff]
        n_after = len(self.df)
        self.logger.info(f'''\
Filtering out TCRs with average abundance less than {self.abundance_cutoff}
Unique TCRs before filtering: {n_before}, after filtering: {n_after}''')

    def write_fasta(self):
        with FastaWriter(file=f'{self.workdir}/cdr3.faa', mode='w') as fasta:
            for tcr in self.df.index:
                fasta.write(header=tcr, sequence=tcr)
        
    def run_mmseqs(self):
        lines = [
            'mmseqs easy-cluster',
            f'{self.workdir}/cdr3.faa',
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

        unique_seqs = self.mmseqs_df['representative_cdr3_amino_acid_sequence'].unique()
        d = {seq: f'motif_{count+1:06d}' for count, seq in enumerate(unique_seqs)}
        self.mmseqs_df['motif_id'] = self.mmseqs_df['representative_cdr3_amino_acid_sequence'].map(d)
        self.mmseqs_df.drop(columns=['representative_cdr3_amino_acid_sequence'], inplace=True)
        self.mmseqs_df.to_csv(f'{self.outdir}/mmseqs.tsv', index=False, sep='\t')

    def generate_sequence_logo(self):
        dstdir = f'{self.outdir}/sequence-logo'
        os.makedirs(dstdir, exist_ok=True)
        mafft_dir = f'{self.workdir}/mafft-tmp'
        os.makedirs(mafft_dir, exist_ok=True)

        for motif_id, df in self.mmseqs_df.groupby('motif_id', sort=True):
            sequences = df['cdr3_amino_acid_sequence'].tolist()

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
            color_scheme['-'] = np.array([0.9, 0.9, 0.9])  # light gray

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
