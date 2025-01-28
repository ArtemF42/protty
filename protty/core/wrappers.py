import shutil
import subprocess
import logging

logger = logging.getLogger(__name__)


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
            subprocess.run(
                (self.executable, *map(str, args)),
                capture_output=True,
                check=True,
                text=True,
            )
        except subprocess.CalledProcessError as error:
            logger.error(f'Failed to run {self.name}.', exc_info=error)

    @property
    def name(self) -> str:
        raise NotImplementedError


class ClustalOmega(_BaseWrapper):
    """A simple wrapper for Clustal Omega (http://www.clustal.org/omega/)."""

    def run(self, infile: str, outfile: str, threads: int) -> bool:
        """Executes Clustal Omega with the specified input file, output file, and number
        of threads.

        Args:
            infile (str): Path to the input file.
            outfile (str): Path to the output file.
            threads (int): Number of threads to use.

        Returns:
            bool: True if the program completed successfully, False otherwise.
        """
        super().run('-i', infile, '-o', outfile, '--threads', threads)

    @property
    def name(self) -> str:
        return 'Clustal Omega'
