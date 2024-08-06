from argparse import ArgumentParser
from ftplib import FTP
import logging
import os
import re

from pyhmmer.easel import Alphabet, MSAFile
from pyhmmer.plan7 import Background, Builder

from protty.clustalo import ClustalOmega

logging.basicConfig(format='%(levelname)s: %(message)s', level=logging.INFO)


def download_merops_database(outdir: str) -> None:
    HOST = 'ftp.ebi.ac.uk'
    PATH = 'pub/databases/merops/current_release/seqlib/'

    with FTP(HOST) as ftp:
        ftp.login()
        ftp.cwd(PATH)

        for filename in ftp.nlst():
            if re.fullmatch(r'[acgimnpstu]\d+\.lib', filename):
                with open(f'{outdir}/{filename}', 'wb') as file:
                    ftp.retrbinary(f'RETR {filename}', file.write)
                
                logging.info(f'File {filename} was successfully downloaded')
    

def main() -> None:
    args = parser.parse_args()

    for subdir in ('sequences', 'alignments', 'hmms'):
        os.makedirs(f'{args.outdir}/{subdir}')
    
    download_merops_database(f'{args.outdir}/sequences')
    
    clustalo = ClustalOmega(args.clustalo)

    alphabet = Alphabet.amino()
    background = Background(alphabet)
    builder = Builder()

    for family in map(lambda filename: os.path.splitext(filename)[0],
                      os.listdir(f'{args.outdir}/sequences')):
        clustalo(f'{args.outdir}/sequences/{family}.lib',
                 f'{args.outdir}/alignments/{family}.fasta',
                 threads=args.threads)
        
        with MSAFile(f'{args.outdir}/alignments/{family}.fasta',
                     digital=True,
                     alphabet=alphabet) as file:
            msa = file.read()

        msa.name = family        
        hmm, _, _ = builder.build_msa(msa, background)

        with open(f'{args.outdir}/hmms/{family}.hmm') as file:
            hmm.write(file)

        logging.INFO(f'Successfully built profile HMM for {family.upper()}')


parser = ArgumentParser()

parser.add_argument('--clustalo', default='clustalo',
                    help='Path to Clustal Omega executable')
parser.add_argument('--threads', default=None,
                    help='Number of processors to use')
parser.add_argument('outdir', help='Directory to store all result files')
