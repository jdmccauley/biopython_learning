#!/usr/bin/env python

import numpy as np

from Bio import SeqIO
from Bio.Align import PairwiseAligner, Alignment, Alignments, reverse_complement, read
from Bio.Align.substitution_matrices import Array


"""
This is on pairwise alignments.
"""

def pairwise_align():
    """
    Pairwise alignment can be done on two sequences to return an Alignment
    object.
    """
    seq1 = "GAACT"
    seq2 = "GAT"

    # Default match_score is 1.
    aligner = PairwiseAligner()
    alns = aligner.align(seq1, seq2)
    print("Alignments:")
    for aln in alns:
        print("Subject: ", aln.target, "Query: ", aln.query)
        print("Score: ", aln.score)
        print(aln)

    seq3 = "AAACAAA"
    seq4 = "AAAGAAA"
    alns2 = aligner.align(seq3, seq4)
    for aln in alns2:
        print("Subject: ", aln.target, "Query: ", aln.query)
        print("Score: ", aln.score)
        print(aln)

    print("Now local alignment, see that subject/query may be trimmed:")
    aligner.mode = "local"
    seq5 = "AGAACTC"
    seq6 = "GAACT"
    alns3 = aligner.align(seq5, seq6)
    for aln in alns3:
        print("Subject: ", aln.target, "Query: ", aln.query)
        print("Score: ", aln.score)
        print(aln)

    print("And here's the aligner and all its params:")
    print(aligner)

    print("The alg:")
    print(aligner.algorithm)


def scoring():
    """
    A PairwiseAligner can have EITHER a substitution matrix, OR at least one of
    a match score and mismatch score.
    """
    from Bio.Align import substitution_matrices
    aligner = PairwiseAligner()
    print(
        "matrix: ", 
        aligner.substitution_matrix, 
        "match score: ", 
        aligner.match_score,
        "mismatch score: ",
        aligner.mismatch_score
    )

    print("Now set sub matrix to BLASTN:")
    aligner.substitution_matrix = substitution_matrices.load("BLASTN")
    print(
        "matrix: ", 
        aligner.substitution_matrix, 
        "match score: ", 
        aligner.match_score,
        "mismatch score: ",
        aligner.mismatch_score
    )


    print("Now set match score to 5:")
    aligner.match_score = 5.0
    print(
        "matrix: ", 
        aligner.substitution_matrix, 
        "match score: ", 
        aligner.match_score,
        "mismatch score: ",
        aligner.mismatch_score
    )

    print("You can also use a custom function for gaps, taking start and len as args.")
    def gap_function(start, length):
        if start == 2:
            return -1000
        return -1 * length
    
    print("Before -1000 penalty for gap at 2:")
    alns = aligner.align("AACTT", "AATT")
    print("Alns score: ", alns.score)
    for aln in alns:
        print(aln)
        print("Score: ", aln.score)
    
    print("After penalty for -1000 at 2:")
    aligner.query_gap_score = gap_function
    alns = aligner.align("AACTT", "AATT")
    print("Alns score: ", alns.score)
    for aln in alns:
        print(aln)
        print("Score: ", aln.score)
    for aln in alns:
        print(aln)
        print("Score: ", aln.score)


def reverse_alignments():
    """
    You can do alignments of the rcomp of a target and query sequence,
    you just need to rcomp the query, NOT THE TARGET, and set arg 
    "strand='-'". The strand is what says to use the rcomp of the target,
    but you still need to rcomp the query manually.

    NOTE that rcomp alignments may have different scores than original
    because of differring right/left gap penalties. Make them equal to
    keep alignments the same.
    """
    target = "AAAACCC"
    query = "AACC"
    aligner = PairwiseAligner()
    # I'm setting this so there are no gaps if possible.
    aligner.query_internal_open_gap_score = -2.0

    print("Normal alignment score: ")
    print(next(aligner.align(target, query)))

    print("Using rcomp query but no neg strand: ")
    print(next(aligner.align(target, reverse_complement(query))))

    print("Using rcomp query and neg strand: ")
    print(next(aligner.align(target, reverse_complement(query), strand="-")))
    print("Notice when 'strand='-'' that bases are like the original.")
    print("So this is likely for when rcomps are available for either target or query.")


def sub_mats():
    """
    substitution_matrices are subclasses of np arrays that can use indices of
    bases. As such, you can do np things with them.
    Just be sure to set the alphabet correctly.
    """
    # Note that dims=1 is just for a frequency count array, not a sub matrix.
    arr = Array(
        "ATGC",
        dims=2,
        data = np.arange(16).reshape(4, 4)
    )

    print("Using range of 0-16:")
    print(arr)

    print("arr['A', 'G']:")
    print(arr['A','G'])

    print("arr + arr:")
    print(arr + arr)


# THIS IS A REALLY GOOD ONE
def get_sub_mat_from_aln():
    """
    With some quick math, we can find out what the substitution matrix
    is for a given alignment.
    """
    ecoli = SeqIO.read("ecoli.fa", "fasta")
    bsubtilis = SeqIO.read("bsubtilis.fa", "fasta")

    # Let's use the blastn substitution matrix first for alignment.
    aligner = PairwiseAligner(scoring="blastn")
    # Just take the first alignment.
    aln = next(aligner.align(ecoli, bsubtilis))

    subs = aln.substitutions
    print("Raw match counts:")
    print(subs)

    # We want a symmetric normalized matrix, so take avg of
    # norm and transposed norm.

    total = subs.sum()
    norm = subs / total
    sym_subs = (norm + norm.transpose()) / 2.0

    print("Symmetric substitution matrix:")
    print(format(sym_subs, "%4f"))

    # Now let's get the full prob of each base on their own, no priors.
    background: Array = sym_subs.sum(0)
    print("Background probs:")
    print(format(background, "%4f"))

    # And now lets take the dot product to get the expected sub matrix.
    expected = background[:, None].dot(background[None, :])
    print("Expected sub matrix:")
    print(format(expected, "%4f"))

    print("Ratios on log2 scale:")
    print(np.log2(sym_subs / expected))


def file_sub_mat_from_aln():
    """
    The same process as above, but from a file.
    """
    aln: Alignment = read("protein.aln", "clustal")

    # Just a subset of the amino acids for this example.
    subs: Array = aln.substitutions.select("DEHKR")

    total = subs.sum()

    # Make the symmetric matrix.
    norm = subs / total
    symmetric = (
        norm + norm.transpose()
    ) / 2.0
    print("Norm symmetric:")
    print(format(symmetric, "%4f"))

    # Make the priors matrix.
    priors = symmetric.sum(0)
    print("Priors:")
    print(format(priors, "%4f"))

    # Get the expeected frequencies.
    expected = priors[:, None].dot(priors[None, :])

    print("Frequency change on log scale:")
    print(format(
        np.log2(symmetric / expected), "%4f"
    ))


    

    


def main():
    pairwise_align()
    scoring()
    reverse_alignments()
    sub_mats()
    get_sub_mat_from_aln()
    file_sub_mat_from_aln()


if __name__ == "__main__":
    main()

