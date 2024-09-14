import argparse
import logging
import os
import random
import sys

from Bio import SeqIO

from pyhmmer.easel import Alphabet, MSAFile
from pyhmmer.plan7 import Background, Builder

from .wrappers import ClustalOmega, ProgramNotFoundError


def filter_fasta(filename: str, min_records: int, max_records: int) -> bool:
    with open(f'{OUTDIR}/raw/{filename}', errors='ignore') as file:
        records = list(SeqIO.parse(file, 'fasta'))

    if len(records) < min_records:
        logging.info(f'File {filename} was ignored (< {min_records} records)')
        return False
    else:
        if len(records) > max_records:
            records = random.sample(records, max_records)
            logging.info(f'File {filename} was trimmed to {max_records} records')
        
        SeqIO.write(records, f'{OUTDIR}/filtered/{filename}', 'fasta')
        return True


def build_profile_hmm(family: str) -> None:
    with MSAFile(f'{OUTDIR}/msa/{family}.fasta', digital=True) as msafile:
        msa = msafile.read()
        msa.name = family.encode()
    
    alphabet = Alphabet.amino()
    builder = Builder(alphabet)
    background = Background(alphabet)

    hmm, *_ = builder.build_msa(msa, background)

    with open(f'{OUTDIR}/profiles/{family}.hmm', 'wb') as hmmfile:
        hmm.write(hmmfile)
    
    logging.info(f'Successfully built profile HMM for {family}')


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
    # parser.add_argument('--skip')
    parser.add_argument('outdir')

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    global OUTDIR
    OUTDIR = args.outdir

    
    for subdir in ('raw', 'filtered', 'msa', 'profiles'):
        os.makedirs(f'{OUTDIR}/{subdir}')
    
    try:
        clustalo = ClustalOmega(args.clustalo)
    except ProgramNotFoundError as error:
        sys.exit('Could not find Clustal Omega')

    download_merops_data()

    for filename in os.listdir(f'{OUTDIR}/raw/'):
        family = os.path.splitext(filename)[0]

        if filter_fasta(filename, args.min, args.max):
            clustalo.run(f'{OUTDIR}/filtered/{filename}',
                         f'{OUTDIR}/msa/{family}.fasta',
                         threads=args.threads)
            
            build_profile_hmm(family)
