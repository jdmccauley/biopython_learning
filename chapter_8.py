#!/usr/bin/env python

from io import StringIO, BytesIO
import subprocess as sp

import numpy as np


from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio import AlignIO

"""
This is on MultipleSequenceAlignments, using AlignIO.
AlignIO is the sister to SeqIO, and should be used for reading/parsing/writing
multiple sequence alignments.
"""

def reading_alns():
    """
    Like SeqIO: read() for one, parse() for multiple. parse() is an iterator.
    """
    alns = AlignIO.parse("opuntia.aln", "clustal")
    for aln in alns:
        print(aln.alignment)


def ambig_fasta_alns():
    """
    Sometimes people write alignments to fasta, which are hard to tell their
    type: multiple alignment, multiple multiple, or multiple pair. See examples
    of each and how to handle them.

    Know your fasta to know how what the alignments are!
    """
    multiple_aln = AlignIO.parse("multiple.fa", "fasta")
    multiple_multiple_aln = AlignIO.parse(
        "multiple_multiple.fa", "fasta", seq_count=3
    )
    multiple_pair_aln = AlignIO.parse(
        "multiple_pair.fa", "fasta", seq_count=2
    )

    print("Multiple:")
    for aln in multiple_aln:
        print(aln)
    print()

    print("Multiple multiple:")
    for aln in multiple_multiple_aln:
        print(aln)
    print()

    print("Multiple pair:")
    for aln in multiple_pair_aln:
        print(aln)
    print()


def aln_writing():
    """
    How to write MultipleSeqAlignments.
    """
    aln1 = MultipleSeqAlignment([
        SeqRecord(Seq("ATGC"), "alpha"),
        SeqRecord(Seq("ATG-"), "beta"),
        SeqRecord(Seq("-TGC"), "gamma")
    ])

    aln2 = MultipleSeqAlignment([
        SeqRecord(Seq("AAAAGGGG"), "alpha"),
        SeqRecord(Seq("----GGGG"), "beta"),
        SeqRecord(Seq("AAAA----", "gamma"))
    ])

    alns = [aln1, aln2]
    s = StringIO()
    
    AlignIO.write(alns, s, "phylip")

    s.seek(0)

    print("Alns as phylip:")
    print(s.read())


def aln_conversion():
    """
    This shows how to read/write alignments and convert between them, which
    is a big use case for reading/writing alignments in bulk.
    """
    alns = AlignIO.parse("PF05371_seed.sth", "stockholm")
    s = StringIO()
    AlignIO.write(alns, s, "clustal")

    s.seek(0)

    print(s.read())

    # If you really want to convert between files, just do convert:
    # AlignIO.convert("PF05371_seed.sth", "stockholm", "seed.aln", "clustal")


def aln_slicing():
    """
    You can slice alignments just like you do sequences.
    """
    aln = AlignIO.read("protein.aln", "clustal")

    print("aln[2, 6]: ", aln[2, 6])
    print("aln[2][6]: ", aln[2][6])
    print("aln[:, 6]: ", aln[:, 6])
    print("aln[3:6, :6]: ", aln[1:4, :8])

    print("You can add alignments, they records just need to have")
    print("the same ids. Because otherwise it doesn't make sense!")
    # You can also only add if same number of rows.

    print("aln:")
    print(aln)

    new_aln = aln[:, 8:11] + aln[:, 14:18]
    print("new aln of columns 8-10 + 14:17")
    print(new_aln)


def misc_aln_ops():
    """
    Miscellaneous alignment operations: getting the alignment as an np.array,
    checking substitutions.
    """
    aln = AlignIO.read("protein.aln", "clustal")

    print("As arr:")
    print(np.array(aln))

    print("Substitutions for DYKA:")
    print(aln.substitutions.select("DYKA"))


def aln_summary():
    """
    The AlignInfo.SummaryInfo class lets us see a consensus sequence,
    position specific score of alignment, other info, and get the sub
    matrix.
    """
    # Currently in 1.82 there's Deprecation warnings, so just ignore them
    # for now.
    import warnings
    warnings.filterwarnings("ignore")

    # Just until 10th residue so we can print the whole thing easily.
    aln = AlignIO.read("protein.aln", "clustal")[:, :10]
    summary = AlignInfo.SummaryInfo(aln)

    # Can set theshold to be % needed to declare consensus,
    # and what the ambig char is.
    print("Consensus:")
    print(summary.dumb_consensus())

    print("Position specific score matrix:")
    print(summary.pos_specific_score_matrix())


def get_real_sub_matrix():
    """
    Like example from last chapter: get the real observed substitution matrix.
    """
    aln = AlignIO.read("protein.aln", "clustal")

    subs = aln.substitutions.select("DEHKR")
    total = subs.sum()

    norm: np.array  = subs / total
    sym_mat: np.array = (norm + norm.transpose()) / 2

    priors = sym_mat.sum(0)

    expected = priors[:, None].dot(priors[None, :])

    print("Expected:")
    print(expected)

    diff = np.log2(sym_mat / expected)

    print(format(diff, "%4f"))  

    # print(format(priors, "%4f"))


def using_outside_commands():
    """
    To use outside command line tools, use subprocess.
    We'll do an alignment with clustal omega here.
    Clustal omega can use files, but we'll use a pipe to make it more
    fun.
    """
    # cat the file, but don't do it until we start the clustal
    # command with sp.run. Popen prevents this.
    cat = sp.Popen("cat opuntia.fasta".split(), stdout=sp.PIPE)

    # Keep the output, we'll need it. It's in bytes though.
    clustal_out = sp.run(
        "clustalo -i - --outfmt=clu".split(), 
        stdin=cat.stdout,
        capture_output=True
    )

    # We can use StringIO to use the stdout as a 'file' for AlignIO.
    s = StringIO(clustal_out.stdout.decode("utf-8"))

    aln = AlignIO.read(s, "clustal")

    print(aln)



def main():
    reading_alns()
    ambig_fasta_alns()
    aln_writing()
    aln_conversion()
    aln_slicing()
    misc_aln_ops()
    aln_summary()
    get_real_sub_matrix()
    using_outside_commands()

if __name__ == "__main__":
    main()

