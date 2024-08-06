import os
import subprocess
import shutil


class ClustalOmega:
    '''A simple wrapper around Clustal Omega (http://www.clustal.org/omega/)'''

    def __init__(self, executable: str = 'clustalo') -> None:
        if shutil.which(executable):
            self.executable = executable
        else:
            raise FileNotFoundError(f'failed to find Clustal Omega')
    
    def __call__(self, infile: str, outfile: str, threads: int = None) -> None:
        threads = str(threads if threads else os.cpu_count())
        subprocess.run((self.executable,
                        '--infile', infile,
                        '--outfile', outfile,
                        '--threads', threads))
