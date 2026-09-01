#!/usr/bin/env python3

"""
Reads a FASTA file and if >1 sequence has the same description line,
it only keeps the longest sequence. It outputs all the sequences to stdout
when complete.
"""
from Bio import SeqIO
import sys

if len(sys.argv) == 1:
    print("ERROR: Please enter the filename to read as the first argument after the program name")
    sys.exit()
else:
    file = sys.argv[1]

seqs = {}
new = 0
existing = 0
for seq_record in SeqIO.parse(file, "fasta"):
    if seq_record.name not in seqs:
        seqs[seq_record.name] = seq_record.seq
        new += 1
    else:
        existing += 1
        if len(seqs[seq_record.name]) <= len(seq_record.seq):
            seqs[seq_record.name] = seq_record.seq

for name, seq in seqs.items():
    print(">" + name)
    print(seq)