from argparse import ArgumentParser

from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile
from pyhmmer.hmmer import hmmscan

from Bio import SeqIO
from Bio.SeqUtils import molecular_weight
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint

FIELDS = ('seqid', 'family', 'evalue', 'mw', 'pi')


def main() -> None:
    args = parser.parse_args()

    sequences = SeqIO.index(args.seqfile, 'fasta')

    with (open(args.tsv, 'w') as tsv_outfile,
          open(args.faa, 'w') as faa_outfile,
          SequenceFile(args.seqfile, digital=True) as queries,
          HMMFile(args.hmmfile) as profiles):
        print(*FIELDS, sep='\t', file=tsv_outfile)

        for hits in hmmscan(queries, profiles):

            if hits:
                seqid = hits.query_name.decode()
                family = hits[0].name.decode()
                evalue = hits[0].evalue

                if 'X' in (seq := sequences[seqid]):
                    mw, pi = None, None
                else:
                    mw = molecular_weight(seq, seq_type='protein')
                    pi = IsoelectricPoint(seq).pi()

                print(seqid, family, evalue, mw, pi, sep='\t', file=tsv_outfile)
                SeqIO.write(seq, faa_outfile, 'fasta')



parser = ArgumentParser()

parser.add_argument('--tsv', default='./predicted_proteases.tsv')
parser.add_argument('--faa', default='./predicted_proteases.faa')
parser.add_argument('hmmfile')
parser.add_argument('seqfile')
