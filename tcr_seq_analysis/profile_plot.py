import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from .template import Processor
from .diversity_clonality import gini_coefficient


class ProfilePlot(Processor):

    count_df: pd.DataFrame

    def main(self, count_df: pd.DataFrame):
        self.count_df = count_df

        self.rank_count_plots()
        self.lorenz_curves()
        self.tcr_length_histogram()
    
    def rank_count_plots(self):
        dstdir = f'{self.outdir}/rank-count-plot'
        os.makedirs(dstdir, exist_ok=True)
        
        global_max_count = self.count_df.max().max()
        global_max_ranks = 0
        for sample_id in self.count_df.columns:
            count_values = self.count_df[sample_id]
            n_ranks = len(count_values[count_values != 0])
            if n_ranks > global_max_ranks:
                global_max_ranks = n_ranks

        for sample_id in self.count_df.columns:
            count_values = self.count_df[sample_id]
            nonzero_counts = count_values[count_values != 0]
            sorted_counts = nonzero_counts.sort_values(ascending=False)
            ranks = range(1, len(sorted_counts) + 1)

            plt.figure(figsize=(5/2.54, 5/2.54))
            plt.rcParams.update({'font.size': 7})
            plt.rcParams.update({'lines.linewidth': 0.5})
            plt.plot(ranks, sorted_counts.values, linewidth=1.0, clip_on=False, color='black')
            plt.xlabel('Rank')
            plt.ylabel('Count')
            plt.yscale('log')
            plt.ylim(0.5, global_max_count*1.5)
            plt.xlim(-global_max_ranks*0.05, global_max_ranks*1.05)
            plt.title(sample_id)
            plt.tight_layout()
            plt.savefig(f'{dstdir}/{sample_id}.png', dpi=600)
            plt.close()
    
    def lorenz_curves(self):
        dstdir = f'{self.outdir}/lorenz-curve'
        os.makedirs(dstdir, exist_ok=True)
        
        for sample_id in self.count_df.columns:
            count_values = self.count_df[sample_id]
            nonzero_counts = count_values[count_values != 0]
            sorted_counts = nonzero_counts.sort_values(ascending=True)  # low to high
            cumulative_sum = sorted_counts.cumsum().values
            lorenz_curve = cumulative_sum / max(cumulative_sum)

            ranks = np.arange(1, len(lorenz_curve) + 1)
            normalized_ranks = ranks / max(ranks)
            
            gini = gini_coefficient(count_values)

            plt.figure(figsize=(5.5/2.54, 5.5/2.54))
            plt.rcParams.update({'font.size': 7})
            plt.rcParams.update({'lines.linewidth': 0.5})

            ax = plt.gca()
            for spine in ax.spines.values():
                spine.set_zorder(0)  # push axis lines behind
                ax.patch.set_zorder(-1)  # push axes background behind (safe)
            
            plt.plot(normalized_ranks, lorenz_curve, linewidth=1.0, clip_on=False, color='red', zorder=10)
            plt.plot(normalized_ranks, normalized_ranks, linewidth=1.0, clip_on=False, color='black', linestyle='--', zorder=10)
            plt.xlabel('Rank')
            plt.ylabel('Cumulative Count')
            plt.ylim(0, 1)
            plt.xlim(0, 1)
            ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
            ax.set_xticklabels(['0%', '25%', '50%', '75%', '100%'])
            ax.set_yticklabels(['0%', '25%', '50%', '75%', '100%'])
            plt.title(f'{sample_id}\n(Gini: {gini:.4f})')
            plt.tight_layout()
            plt.savefig(f'{dstdir}/{sample_id}.png', dpi=600)
            plt.close()
    
    def tcr_length_histogram(self):
        dstdir = f'{self.outdir}/tcr-length-histogram'
        os.makedirs(dstdir, exist_ok=True)

        global_tcr_length_min = min([len(sequence) for sequence in self.count_df.index])
        global_tcr_length_max = max([len(sequence) for sequence in self.count_df.index])
        
        for sample_id in self.count_df.columns:
            count_values = self.count_df[sample_id]
            nonzero_counts = count_values[count_values != 0]
            tcr_sequences = nonzero_counts.index
            tcr_lengths = [len(sequence) for sequence in tcr_sequences]
            plt.figure(figsize=(5/2.54, 5/2.54))
            plt.rcParams.update({'font.size': 7})
            plt.rcParams.update({'lines.linewidth': 0.5})
            plt.hist(tcr_lengths, bins=range(global_tcr_length_min, global_tcr_length_max + 1), color='#595959')  
            plt.xlabel('TCR Length')
            plt.ylabel('Count')
            plt.title(sample_id)
            plt.tight_layout()
            plt.savefig(f'{dstdir}/{sample_id}.png', dpi=600)
            plt.close()
