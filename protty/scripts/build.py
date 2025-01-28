import argparse
import logging
import os
import random
import sys

from rich.progress import track

from Bio import SeqIO

from pyhmmer.easel import Alphabet, MSAFile
from pyhmmer.plan7 import Background, Builder

from protty.core import merops
from protty.core.wrappers import ClustalOmega

logger = logging.getLogger(__name__)


def filter_fasta(
    output_dir: str, family: str, min_records: int, max_records: int
) -> bool:
    with open(
        os.path.join(output_dir, 'raw', f'{family}.lib'), errors='ignore'
    ) as file:
        records = list(SeqIO.parse(file, 'fasta'))

    if len(records) < min_records:
        logger.info(
            f'File for family {family} was ignored (fewer than {min_records} records).'
        )
        return False
    else:
        if len(records) > max_records:
            records = random.sample(records, max_records)
            logger.info(
                f'File for family {family} was trimmed (exceeded {max_records} records).'
            )

        SeqIO.write(
            records, os.path.join(output_dir, 'filtered', f'{family}.fasta'), 'fasta'
        )
        return True


def build_profile_hmm(output_dir: str, family: str) -> None:
    with MSAFile(
        os.path.join(output_dir, 'msa', f'{family}.fasta'), digital=True
    ) as msafile:
        msa = msafile.read()
        msa.name = family.encode()

    alphabet = Alphabet.amino()
    builder = Builder(alphabet)
    background = Background(alphabet)

    hmm, *_ = builder.build_msa(msa, background)

    with open(
        os.path.join(output_dir, 'profiles', f'{family}.hmm'), mode='wb'
    ) as hmmfile:
        hmm.write(hmmfile)

    logger.info(f'Successfully built profile HMM for {family}.')


def parse_args() -> argparse.Namespace:
    """Parses command-line arguments.

    Returns:
        argparse.Namespace: An object containing the parsed command-line arguments.
    """
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--min-records',
        default=5,
        type=int,
        help='Files with fewer than the minimum number of records will be ignored (default: 5).',
        metavar='INT',
    )
    parser.add_argument(
        '--max-records',
        default=1000,
        type=int,
        help='Files with more than the maximum number of records will be randomly trimmed (default: 1000).',
        metavar='INT',
    )
    parser.add_argument(
        '--seed',
        default=42,
        type=int,
        help='Random seed for subsampling (default: 42).',
        metavar='INT',
    )
    parser.add_argument(
        '--clustalo',
        default='clustalo',
        type=str,
        help='Path to the Clustal Omega executable (default: clustalo).',
        metavar='EXECUTABLE',
    )
    parser.add_argument(
        '--threads',
        default=os.cpu_count(),
        type=int,
        help='Number of threads to use for Clustal Omega (default: number of CPU cores).',
        metavar='INT',
    )
    parser.add_argument(
        '--output-dir',
        default='.',
        type=str,
        help='Path to the output directory (default: current working directory).',
        metavar='PATH',
    )

    return parser.parse_args()


def main() -> None:
    # Parse command line arguments
    args = parse_args()

    # Set random seed
    random.seed(args.seed)

    # Create the output directory if necessary
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # Configure the logger
    logging.basicConfig(
        filename=os.path.join(args.output_dir, 'protty_build.log'),
        filemode='w',
        format='%(asctime)s %(levelname)s %(message)s',
        level=logging.INFO,
    )

    # Initialize Clustal Omega
    clustalo = ClustalOmega(args.clustalo)

    # Create subdirectories for output
    for sub_dir in ('raw', 'filtered', 'msa', 'profiles'):
        os.mkdir(os.path.join(args.output_dir, sub_dir))

    # Download MEROPS database
    if not merops.download(os.path.join(args.output_dir, 'raw')):
        sys.exit(
            'An error occurred while downloading MEROPS data. See log for details.'
        )

    # Filter FASTA files, align protein sequences and build profile HMMs
    for filename in track(
        os.listdir(os.path.join(args.output_dir, 'raw')), description='Building HMMs...'
    ):
        family, _ = os.path.splitext(filename)

        if filter_fasta(args.output_dir, family, args.min_records, args.max_records):
            clustalo.run(
                os.path.join(args.output_dir, 'filtered', f'{family}.fasta'),
                os.path.join(args.output_dir, 'msa', f'{family}.fasta'),
                args.threads,
            )
            build_profile_hmm(args.output_dir, family)

    # Merge files
    with open(os.path.join(args.output_dir, 'merops.hmm'), mode='w') as outfile:
        for filename in os.listdir(os.path.join(args.output_dir, 'profiles')):
            with open(os.path.join(args.output_dir, 'profiles', filename)) as infile:
                for line in infile:
                    outfile.write(line)
