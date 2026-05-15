#!/usr/bin/env python3
import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True, help="Input FASTA alignment")
parser.add_argument("-o", "--output", required=True, help="Output FASTA alignment with gaps")
parser.add_argument("-p", "--gap-prob", type=float, default=0.05,
                    help="Probability of replacing each character with a gap")
parser.add_argument("--seed", type=int, default=123)
args = parser.parse_args()

random.seed(args.seed)

records_out = []

for rec in SeqIO.parse(args.input, "fasta"):
    seq = str(rec.seq)
    gapped = "".join(
        "-" if random.random() < args.gap_prob else aa
        for aa in seq
    )

    records_out.append(
        SeqRecord(
            Seq(gapped),
            id=rec.id,
            description=""
        )
    )

SeqIO.write(records_out, args.output, "fasta")