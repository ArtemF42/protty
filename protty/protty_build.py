import argparse
import ftplib
import logging
import os
import random
import re
import sys

from Bio import SeqIO

from .wrappers import ClustalOmega, HMMER3


def download_merops_data() -> None:
    HOST = 'ftp.ebi.ac.uk'
    PATH = 'pub/databases/merops/current_release/seqlib/'

    try:
        with ftplib.FTP(HOST) as ftp:
            ftp.login()
            ftp.cwd(PATH)

            for filename in ftp.nlst():
                if re.fullmatch(r'[acgimnpstu]\d+\.lib', filename):
                    with open(f'{OUTDIR}/raw/{filename}', 'wb') as file:
                        ftp.retrbinary(f'RETR {filename}', file.write)
                    
                    logging.info(f'File {filename} was successfully downloaded')
    except ftplib.all_errors as error:
        # logging.error()
        sys.exit('Failed to download MEROPS data')


def filter_fasta(filename: str, min_sequences: int, max_sequences: int) -> bool:
    records = list(SeqIO.parse(f'{OUTDIR}/raw/{filename}', 'fasta'))

    if len(records) > min_sequences:
        if len(records) > max_sequences:
            records = random.sample(records, max_sequences)
            logging.info(f'File {filename} was trimmed')
        
        SeqIO.write(records, f'{OUTDIR}/filtered/{filename}', 'fasta')
        return True
    else:
        logging.info(f'File {filename} was ignored (< {min_sequences} sequences)')
        return False


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument('--min', default=5)
    parser.add_argument('--max', default=1000)
    parser.add_argument('--threads', '-t', default=os.cpu_count())
    parser.add_argument('outdir')

    return parser.parse_args()


def main() -> None:
    args = parse_arguments()

    clustalo = ClustalOmega()
    hmmer3 = HMMER3()

    global OUTDIR
    OUTDIR = args.outdir

    logging.basicConfig(filename=f'{OUTDIR}/protty.log', filemode='w',
                        format='%(levelname)s: %(message)s', level=logging.INFO)

    for subdir in ('raw', 'filtered', 'msa', 'profiles'):
        os.makedirs(f'{OUTDIR}/{subdir}')

    download_merops_data()

    for filename in os.listdir(f'{OUTDIR}/raw'):
        if filter_fasta(filename, args.min, args.max):
            family, _ = os.path.splitext(filename)
            
            clustalo.align(f'{OUTDIR}/filtered/{filename}',
                           f'{OUTDIR}/msa/{family}.alignment.fasta',
                           threads=args.threads)
            hmmer3.build(f'{OUTDIR}/profiles/{family}.hmm',
                         f'{OUTDIR}/msa/{family}.alignment.fasta', family)
            
    with open(f'{OUTDIR}/merops.hmm', 'w') as infile:
        for filename in os.listdir(f'{OUTDIR}/profiles'):
            with open(f'{OUTDIR}/profiles/{filename}') as outfile:
                for line in infile:
                    outfile.write(line)
