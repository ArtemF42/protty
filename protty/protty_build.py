import argparse
import ftplib
import logging
import os
import random
import re
import sys

from Bio import SeqIO

from .wrappers import ClustalOmega, ProgramNotFoundError


def download_merops_data() -> None:
    HOST = 'ftp.ebi.ac.uk'
    PATH = 'pub/databases/merops/current_release/seqlib/'

    logging.info('Started to download MEROPS data')

    try:
        with ftplib.FTP(HOST) as ftp:
            ftp.login()
            ftp.cwd(PATH)

            for filename in ftp.nlst():
                if re.match(r'[acgimnpstu]\d+\.lib', filename):
                    with open(f'{OUTDIR}/raw/{filename}', 'wb') as file:
                        ftp.retrbinary(f'RETR {filename}', file.write)
            
            logging.info('Successfully downloaded MEROPS data')
    except ftplib.all_errors as error:
        sys.exit('Failed to download MEROPS data')


def filter_fasta(filename: str, min_records: int, max_records: int) -> bool:
    records = list(SeqIO.parse(f'{OUTDIR}/raw/{filename}', 'fasta'))

    if len(records) < min_records:
        logging.info(f'File {filename} was ignored (< {min_records} records)')
        return False
    else:
        if len(records) > max_records:
            records = random.sample(records, max_records)
            logging.info(f'File {filename} was trimmed to {max_records} records')
        
        SeqIO.write(records, f'{OUTDIR}/filtered/{filename}', 'fasta')
        return True


def build_profile() -> None:
    ...



def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument('--min', default=5, type=int)
    parser.add_argument('--max', default=1000, type=int)
    parser.add_argument('--clustalo', default='clustalo')
    parser.add_argument('--threads', default=os.cpu_count(), type=int)
    parser.add_argument('outdir')

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    global OUTDIR
    OUTDIR = args.outdir

    logging.basicConfig(filename=f'{OUTDIR}/protty_build.log',
                        filemode='w',
                        format='%(asctime)s %(levelname)s %(message)s',
                        level=logging.INFO)
    
    for subdir in ('raw', 'filtered', 'msa', 'profiles'):
        os.makedirs(f'{OUTDIR}/{subdir}')
    
    try:
        clustalo = ClustalOmega(args.clustalo)
    except ProgramNotFoundError as error:
        sys.exit('Could not find Clustal Omega')

    download_merops_data()

    for filename in os.listdir(f'{OUTDIR}/raw/'):
        if filter_fasta(filename, args.min, args.max):
            pass

    
