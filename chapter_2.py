#!/usr/bin/env python

from io import StringIO
import requests

from Bio import SeqIO

def read_orchids():
    ls_orchid_fasta = StringIO(
        requests.get("https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.fasta").text
    )

    ls_orchid_gb = StringIO(
        requests.get("https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/ls_orchid.gbk").text
    )

    print("Using fastas...")
    for record in SeqIO.parse(ls_orchid_fasta, "fasta"):
        print(record.id)
        print(record.seq)
        print(len(record))

    print("Using genbanks...")
    for record in SeqIO.parse(ls_orchid_gb, "genbank"):
        print(record.id)
        print(record.seq)
        print(len(record))


def main():
    read_orchids()


if __name__ == "__main__":
    main()
