#!/usr/bin/env python

"""
We're going to download a fasta, do an alignment on it, convert the
alignment into a tree, then draw the tree with matplotlib.
"""

from io import StringIO, BytesIO
import subprocess as sp
import tempfile
import os

import shutil

import numpy as np
import requests

from Bio import AlignIO, Phylo, SeqIO

def cleanup(tmpdir):
    shutil.rmtree(tmpdir)


def get_the_fastas(tmpdir):
    """
    Get fastas from github and write to a file in tmpdir.
    """
    res = requests.get(
        "https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples/opuntia.fasta"
    )

    if res.status_code != 200:
        raise ValueError("Not 200.")

    with open(f"{tmpdir}/dna.fasta", "w") as f:
        f.write(res.text)


    return f"{tmpdir}/dna.fasta"


def make_aln(tmpdir, fasta):
    """
    Align the fasta with clustal omega, convert to phylip in tmpdir.
    """
    out = sp.run(
        f"clustalo -i {fasta} --outfmt=clu".split(),
        capture_output=True
    )

    if out.returncode != 0:
        print("Fail alignment.")
        print(out.stderr)
        cleanup(tmpdir)
        raise ValueError

    with open(f"{tmpdir}/dna.aln", 'w') as f:
        f.write(out.stdout.decode("utf-8"))

    AlignIO.convert(f"{tmpdir}/dna.aln", "clustal", f"{tmpdir}/dna.phy", "phylip")
    
    return f"{tmpdir}/dna.phy"


def make_tree(tmpdir, aln):
    """
    Use phyml to make a newick tree from a phylip alignment.
    """
    cur = os.getcwd()

    os.chdir(tmpdir)

    out = sp.run(
        f"phyml -i {aln}".split(),
        capture_output=True
    )

    if out.returncode != 0:
        print("Fail tree.")
        print(out.stderr)
        os.chdir(cur)
        cleanup(tmpdir)
        raise ValueError
    

    os.chdir(cur)

    return f"{tmpdir}/dna.phy_phyml_tree.txt"


def plot_tree(phyml):
    """
    Plot that tree.
    """
    tree = Phylo.read(phyml, "newick")
    Phylo.draw(tree)


def main():
    tmpdir = tempfile.mkdtemp()
    fasta = get_the_fastas(tmpdir)
    aln = make_aln(tmpdir, fasta)
    phyml = make_tree(tmpdir, aln)
    plot_tree(phyml)
    cleanup(tmpdir)

if __name__ == "__main__":
    main()

