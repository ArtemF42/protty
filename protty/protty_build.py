import argparse
import logging
import os
import random
import sys

from rich.progress import track

from Bio import SeqIO

from pyhmmer.easel import Alphabet, MSAFile
from pyhmmer.plan7 import Background, Builder

from .merops import download
from .wrappers import ClustalOmega, ProgramNotFoundError

logger = logging.getLogger(__name__)


def filter_fasta(outdir: str, filename: str, min_records: int, max_records: int) -> None:
    '''Filters FASTA files
    
    Args:
        outdir (str): output directory
        filename (str): name of the FASTA file
        min_records (int): minimum number of records in FASTA file
        max_records (int): maximum number of records in FASTA file
    '''
    #TODO apply user-defined filter

    with open(f'{outdir}/raw/{filename}', errors='ignore') as file:
        records = list(SeqIO.parse(file, 'fasta'))

    if len(records) < min_records:
        logger.info(f'File {filename} was ignored (< {min_records} records)')
    else:
        if len(records) > max_records:
            records = random.sample(records, max_records)
            logger.info(f'File {filename} was trimmed to {max_records} records')
        
        SeqIO.write(records, f'{outdir}/filtered/{filename}', 'fasta')


def build_profile_hmm(outdir: str, family: str) -> None:
    '''Builds profile HMM for the given family

    Args:
        outdir (str): output directory
        family (str): name of the family
    '''
    with MSAFile(f'{outdir}/msa/{family}.fasta', digital=True) as msafile:
        msa = msafile.read()
        msa.name = family.encode()
    
    alphabet = Alphabet.amino()
    builder = Builder(alphabet)
    background = Background(alphabet)

    hmm, *_ = builder.build_msa(msa, background)

    with open(f'{outdir}/profiles/{family}.hmm', 'wb') as hmmfile:
        hmm.write(hmmfile)
    
    logger.info(f'Successfully built profile HMM for {family}')


def parse_args() -> argparse.Namespace:
    '''Parses command line arguments

    Returns:
        argparse.Namespace: object containing parsed command line arguments
    '''
    parser = argparse.ArgumentParser()

    parser.add_argument('--min', default=5, type=int,
                        help='minimum number of records in FASTA file. ' \
                             'Files with fewer records will be ignored ' \
                             'during filtering (default: 5)')
    parser.add_argument('--max', default=1000, type=int,
                        help='maximum number of records in FASTA file. ' \
                             'Files with larger number of records will ' \
                             'be randomly trimmed (default: 1000)')
    parser.add_argument('--clustalo', default='clustalo',
                        help='path to Clustal Omega executable')
    parser.add_argument('--threads', default=os.cpu_count(), type=int,
                        help='number of threads to use')
    parser.add_argument('outdir')

    return parser.parse_args()


def main() -> None:
    # Parse command line arguments
    args = parse_args()

    # Configure logger
    logging.basicConfig(filename=f'{args.outdir}/protty_build.log',
                        filemode='w',
                        format='%(asctime)s %(levelname)s %(message)s',
                        level=logging.INFO)

    # Create output directories
    for subdir in ('raw', 'filtered', 'msa', 'profiles'):
        path = f'{args.outdir}/{subdir}'

        if not os.path.exists(path):
            os.makedirs(path)

    # STEP 1. Download MEROPS data
    if not download(f'{args.outdir}/raw'):
        sys.exit('An error occurred while downloading MEROPS data. '
                 'See log file for details.')
        
    # STEP 2. Filter FASTA files
    for filename in track(os.listdir(f'{args.outdir}/raw'),
                          description='Filtering...'):
        filter_fasta(args.outdir, filename, args.min, args.max)
    
    # STEP 3. Align protein sequences & build profile HMMs
    try:
        clustalo = ClustalOmega(args.clustalo)
    except ProgramNotFoundError as error:
        sys.exit('Could not find Clustal Omega')
    
    for filename in track(os.listdir(f'{args.outdir}/filtered'),
                          description='Building HMMs...'):
        family = os.path.splitext(filename)[0]
        clustalo.run(f'{args.outdir}/filtered/{filename}',
                     f'{args.outdir}/msa/{family}.fasta',
                     threads=args.threads)
        build_profile_hmm(family)
