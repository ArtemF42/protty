import argparse

from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile
from pyhmmer.hmmer import hmmscan


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument('--tsv', default='./predicted_proteases.tsv')
    parser.add_argument('--faa', default='./predicted_proteases.faa')
    parser.add_argument('--evalue', '-e', default=1e-3, type=int)
    parser.add_argument('hmmfile')
    parser.add_argument('seqfile')

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    sequences = SeqIO.index(args.seqfile, 'fasta')

    with (SequenceFile(args.seqfile, digital=True) as queries,
          HMMFile(args.hmmfile) as profiles,
          open(args.tsv, mode='w') as tsv_outfile,
          open(args.faa, mode='w') as faa_outfile):
        print('seq_id', 'family', 'molecular_weight', 'isoelectric_point',
              sep='\t', file=tsv_outfile)

        for hits in hmmscan(queries, profiles, E=args.evalue):
            if hits:
                seq_id = hits.query_name.decode()
                family = hits[0].name.decode()

                mw, ip = None, None

                if not 'X' in (seq := sequences[seq_id]):
                    mw = round(molecular_weight(seq, seq_type='protein'))
                    ip = round(IsoelectricPoint(seq).pi(), 2)
                
                print(seq_id, family, mw, ip, sep='\t', file=tsv_outfile)
                SeqIO.write(seq, faa_outfile, 'fasta')
