from ftplib import FTP, all_errors
from logging import getLogger
from rich.progress import track
import re

logger = getLogger(__name__)

HOST = 'ftp.ebi.ac.uk'
PATH = 'pub/databases/merops/current_release/seqlib/'


def download(dst: str) -> bool:
    '''Downloads MEROPS data to the specified directory

    Args:
        dst (str): output directory

    Returns:
        bool: True if the database was downloaded, False otherwise
    '''
    logger.info('Started to download MEROPS data')

    try:
        with FTP(HOST) as ftp:
            ftp.login()
            ftp.cwd(PATH)

            for filename in track(ftp.nlst(), description='Downloading...'):
                if re.match(r'[acgimnpstu]\d+\.lib', filename):
                    with open(f'{dst}/{filename}', mode='wb') as file:
                        ftp.retrbinary(f'RETR {filename}', file.write)
                
            logger.info('Successfully downloaded MEROPS data')
            return True
    except all_errors as error:
        logger.error('Failed to download MEROPS data', exc_info=error)
        return False
