import logging
import shutil
import subprocess


class ProgramNotFoundError(Exception):
    pass


class _BaseWrapper:
    def __init__(self, executable: str) -> None:
        if shutil.which(executable):
            self.executable = executable
        else:
            raise ProgramNotFoundError(self.name)

    def run(self, *args: tuple) -> bool:
        try:
            subprocess.run((self.executable, *map(str, args)),
                           capture_output=True, check=True, text=True)
        except subprocess.CalledProcessError as error:
            #TODO
            return False
        else:
            return True

    @property
    def name(self) -> str:
        raise NotImplementedError


class ClustalOmega(_BaseWrapper):
    '''A simple wrapper around Clustal Omega (http://www.clustal.org/omega/)'''
    
    def run(self, infile: str, outfile: str, threads: int) -> bool:
        super().run('-i', infile, '-o', outfile, '--threads', threads)
    
    @property
    def name(self) -> str:
        return 'Clustal Omega'
