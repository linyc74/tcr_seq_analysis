import argparse
import tcr_seq_analysis


__VERSION__ = '1.0.0-beta'


PROG = 'python tcr_seq_analysis'
DESCRIPTION = f'T Cell Receptor Sequencing Analysis Suite (version {__VERSION__}) by Yu-Cheng Lin (ylin@nycu.edu.tw)'
REQUIRED = [
    {
        'keys': ['-c', '--csv-dir'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the directory containing all input CSV files',
        }
    },
    {
        'keys': ['-s', '--csv-suffix'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'suffix of input CSV files',
        }
    },
    {
        'keys': ['-m', '--sample-sheet'],
        'properties': {
            'type': str,
            'required': True,
            'help': 'path to the sample sheet CSV file; the first column is the sample ID',
        }
    },
]
OPTIONAL = [
    {
        'keys': ['-o', '--outdir'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'tcr_seq_analysis_outdir',
            'help': 'path to the output directory (default: %(default)s)',
        }
    },
    {
        'keys': ['--clonal-index-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'cdr3_amino_acid_sequence',
            'help': 'clonal index column in the input CSV files (default: %(default)s)',
        }
    },
    {
        'keys': ['--count-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'read_count',
            'help': 'count column in the input CSV files (default: %(default)s)',
        }
    },
    {
        'keys': ['--group-column'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'Group',
            'help': 'group column in the sample sheet CSV file (default: %(default)s)',
        }
    },
    {
        'keys': ['--rpm-cutoff'],
        'properties': {
            'type': float,
            'required': False,
            'default': 1000,
            'help': 'RPM cutoff before motif clustering (default: %(default)s)',
        }
    },
    {
        'keys': ['--clustering-identity'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.80,
            'help': 'sequence identity for motif clustering (default: %(default)s)',
        }
    },
    {
        'keys': ['--p-value'],
        'properties': {
            'type': float,
            'required': False,
            'default': 0.05,
            'help': 'p value cutoff for motif differential abundance (default: %(default)s)',
        }
    },
    {
        'keys': ['--colormap'],
        'properties': {
            'type': str,
            'required': False,
            'default': 'Set1',
            'help': 'matplotlib colormap for plotting, or comma-separated color names, e.g. "darkred,lightgreen,skyblue" (default: %(default)s)',
        }
    },
    {
        'keys': ['--invert-colors'],
        'properties': {
            'action': 'store_true',
            'help': 'invert the order of colors',
        }
    },
    {
        'keys': ['--publication-figure'],
        'properties': {
            'action': 'store_true',
            'help': 'plot figures in the form and quality for paper publication',
        }
    },
    {
        'keys': ['-t', '--threads'],
        'properties': {
            'type': int,
            'required': False,
            'default': 4,
            'help': 'number of CPU threads (default: %(default)s)',
        }
    },
    {
        'keys': ['-d', '--debug'],
        'properties': {
            'action': 'store_true',
            'help': 'debug mode',
        }
    },
    {
        'keys': ['-h', '--help'],
        'properties': {
            'action': 'help',
            'help': 'show this help message',
        }
    },
    {
        'keys': ['-v', '--version'],
        'properties': {
            'action': 'version',
            'version': __VERSION__,
            'help': 'show version',
        }
    },
]


class EntryPoint:

    parser: argparse.ArgumentParser

    def main(self):
        self.set_parser()
        self.add_required_arguments()
        self.add_optional_arguments()
        self.run()

    def set_parser(self):
        self.parser = argparse.ArgumentParser(
            prog=PROG,
            description=DESCRIPTION,
            add_help=False,
            formatter_class=argparse.RawTextHelpFormatter)

    def add_required_arguments(self):
        group = self.parser.add_argument_group('required arguments')
        for item in REQUIRED:
            group.add_argument(*item['keys'], **item['properties'])

    def add_optional_arguments(self):
        group = self.parser.add_argument_group('optional arguments')
        for item in OPTIONAL:
            group.add_argument(*item['keys'], **item['properties'])

    def run(self):
        args = self.parser.parse_args()
        print(f'Start running TCR Sequencing Analysis Suite version {__VERSION__}\n', flush=True)
        tcr_seq_analysis.main(
            csv_dir=args.csv_dir,
            csv_suffix=args.csv_suffix,
            sample_sheet=args.sample_sheet,
            clonal_index_column=args.clonal_index_column,
            count_column=args.count_column,
            group_column=args.group_column,
            rpm_cutoff=args.rpm_cutoff,
            clustering_identity=args.clustering_identity,
            p_value=args.p_value,
            colormap=args.colormap,
            invert_colors=args.invert_colors,
            publication_figure=args.publication_figure,
            threads=args.threads,
            debug=args.debug,
            outdir=args.outdir)


if __name__ == '__main__':
    EntryPoint().main()
