import ftplib
import logging
import os
import re

from rich.progress import track

logger = logging.getLogger(__name__)

HOST = 'ftp.ebi.ac.uk'
PATH = 'pub/databases/merops/current_release/seqlib/'


def download(dst: str) -> bool:
    """Downloads MEROPS data to the specified directory.

    Args:
        dst (str): Directory for storing downloaded files.

    Returns:
        bool: True if the download was successful, False otherwise.
    """
    logger.info('Started downloading the MEROPS database.')

    try:
        with ftplib.FTP(HOST) as ftp:
            ftp.login()
            ftp.cwd(PATH)

            for filename in track(ftp.nlst(), description='Downloading...'):
                if re.match(r'[a-z]\d+\.lib', filename):
                    with open(os.path.join(dst, filename), mode='wb') as file:
                        ftp.retrbinary(f'RETR {filename}', file.write)

        logger.info('Successfully downloaded the MEROPS database.')
        return True
    except ftplib.all_errors as error:
        logger.error('Failed to download the MEROPS database.', exc_info=error)
        return False
