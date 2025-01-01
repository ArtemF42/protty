import argparse
import os
import sys
from typing import Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from pyhmmer.plan7 import HMMFile
from pyhmmer.easel import SequenceFile
from pyhmmer.hmmer import hmmscan


def compute_physical_properties(seq: Seq) -> Tuple[float, float]:
    molecular_weight, isoelectric_point = None, None

    if not set(seq).difference('ACDEFGHIKLMNPQRSTVWY'):
        analysed_seq = ProteinAnalysis(seq)
        molecular_weight = round(analysed_seq.molecular_weight() / 1000, 1)
        isoelectric_point = round(analysed_seq.isoelectric_point(), 2)

    return molecular_weight, isoelectric_point


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '--tsv',
        default='./predicted_proteases.tsv',
        help='Path to the output file in TSV format '
        '(default is ./predicted_proteases.tsv).',
    )
    parser.add_argument(
        '--faa',
        default='./predicted_proteases.faa',
        help='Path to the output file in FASTA format '
        '(default is ./predicted_proteases.faa).',
    )
    parser.add_argument(
        '--evalue',
        '-e',
        default=1e-3,
        type=float,
        help='E-value threshold for filtering results (default is 1e-3).',
    )
    parser.add_argument('hmmfile', help='Path to the HMM database.')
    parser.add_argument('seqfile', help='Path to the query sequences.')

    args = parser.parse_args()

    if not os.path.exists(args.hmmfile):
        sys.exit(f'File {args.hmmfile} not found.')

    if not os.path.exists(args.seqfile):
        sys.exit(f'File {args.seqfile} not found.')

    return args


def main() -> None:
    args = parse_args()

    with (
        open(args.tsv, mode='w') as tsv_output_file,
        open(args.faa, mode='w') as faa_output_file,
        HMMFile(args.hmmfile) as profiles,
        SequenceFile(args.seqfile, digital=True) as queries,
    ):
        print(
            'accession',
            'family',
            'molecular_weight',
            'isoelectric_point',
            sep='\t',
            file=tsv_output_file,
        )

        for hits in hmmscan(queries, profiles, E=args.evalue):
            if hits:
                accession = hits.query.name.decode()
                family = hits[0].name.decode().upper()

                seq = Seq(hits.query.textize().sequence)
                molecular_weight, isoelectric_point = compute_physical_properties(seq)

                print(
                    accession,
                    family,
                    molecular_weight,
                    isoelectric_point,
                    sep='\t',
                    file=tsv_output_file,
                )
                SeqIO.write(
                    SeqRecord(
                        seq, id=accession, description=f'Putative {family} peptidase'
                    ),
                    faa_output_file,
                    'fasta',
                )
