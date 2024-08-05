from argparse import ArgumentParser
from ftplib import FTP
import logging
import re

logging.basicConfig(format='%(level)s: %(message)s', level=logging.INFO)

HOST = 'ftp.ebi.ac.uk'
PATH = 'pub/databases/merops/current_release/seqlib/'


def download(outdir: str) -> None:
    with FTP(HOST) as ftp:
        ftp.login()
        ftp.cwd(PATH)

        for filename in ftp.nlst():
            if re.fullmatch(r'[acgimnpstu]\d+\.lib', filename):
                with open(f'{outdir}/{filename}', 'wb') as file:
                    ftp.retrbinary(f'RETR {filename}', file.write)
                
                logging.info(f'File {filename} was successfully downloaded')


def main() -> None:
    ...
