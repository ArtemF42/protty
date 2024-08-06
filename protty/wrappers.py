import logging
import os
import subprocess
import shutil
from typing import Tuple


class ProgramNotFoundError(Exception):
    pass


class _BaseWrapper:
    def run(self, executable: str, *args: Tuple[str]) -> None:
        if not shutil.which(executable):
            raise ProgramNotFoundError(self.name)

        try:
            subprocess.run((executable, *map(str, args)), capture_output=True, check=True)
        except subprocess.CalledProcessError as error:
            pass
    
    @property
    def name(self) -> str:
        raise NotImplementedError


class ClustalOmega(_BaseWrapper):
    '''A simple wrapper around Clustal Omega (http://www.clustal.org/omega/)'''

    def align(self, infile: str, outfile: str, threads: int = os.cpu_count()) -> None:
        '''Align sequences'''

        super().run('clustalo', '--threads', threads, '-i', infile, '-o', outfile)
    
    @property
    def name(self) -> str:
        return 'Clustal Omega'


class HMMER3(_BaseWrapper):
    '''A simple wrapper around HMMER3 (http://hmmer.org/)'''

    def build(self, hmmfile: str, msafile: str) -> None:
        '''Construct profiles from multiple sequence alignments'''

        super().run('hmmbuild', hmmfile, msafile)
    
    def press(self, hmmfile: str) -> None:
        '''Prepare a profile database for hmmscan'''

        super().run('hmmpress', hmmfile)
    
    def scan(self, hmmfile: str, seqfile: str, tblout: str, evalue: float = 1e-3) -> None:
        '''Search sequence(s) against a profile database'''
        
        super().run('hmmscan', '--tblout', tblout, '-E', evalue, hmmfile, seqfile)

    @property
    def name(self) -> str:
        return 'HMMER3'
