#!/usr/bin/env python

from io import StringIO, BytesIO

# Numpy used for holding numerical data in alignments.
import numpy as np

from Bio import Align 
from Bio.Align import Alignment, Alignments
from Bio.SeqUtils import gc_fraction

"""
The main package for alignment is Bio.Align, and the main class is Alignment.
"""

def manual_alignment():
    """
    Creating a manual alignment with sequences and their coordinates of
    alignment.
    """
    seqA = "CCGGTTTTT"
    seqB = "AGTTTAA"
    seqC = "AGGTTT"
    seqs = [seqA, seqB, seqC]

    # Now let's manually define the alignment.
    coords = np.array([
        [1, 3, 4, 7, 9],
        [0, 2, 2, 5, 5],
        [0, 2, 3, 6, 6]
    ])
    # This shows seqA[1:3], seqB[0:2], and seqC[0:2] are aligned
    # seqA[3:4], seqC[2:3] are aligned, wuth gap in seqB.
    # seqA[4:7], seqB[2:5], adn seqC[3:6] are aligned.
    # seqA[7:9] has no aligned bases in seqB or seqC.
    aln = Alignment(seqs, coords)
    print(aln)

    return aln

def infer_coordinates():
    """
    Using sequences to infer the coordinates of alignment.
    """
    # Assuming a given alignment from an aligner using dashes for missing
    # bases in alignment, like from manual_alignment:
    #   1 CGGTTTTT 9
    #   0 AG-TTT-- 5
    #   0 AGGTTT-- 6
    aligned_seqs = [
        "CGGTTTTT", "AG-TTT--", "AGGTTT--"
    ]
    seqs = [seq.replace("-", "") for seq in aligned_seqs]
    coords = Alignment.infer_coordinates(aligned_seqs)
    print(coords)

    # But these are not the original coordinates because the first base in
    # seqA was removed in the alignment, as was the last 2 in seqB.
    aligned_seqs[0] = "C" + aligned_seqs[0]
    aligned_seqs[1] = aligned_seqs[1] + "AA"

    # And to update the coords.
    coords[0, :] += 1

    print(coords)

    aln = Alignment(seqs, coords)
    print(aln)


def np_arr_as_sequences():
    """
    This is to show that np.array(i32) can be used as sequences in alignments.
    """
    letters = ['A', 'T', 'G', 'C']

    seqs_char = np.array([
        letters, letters
    ])

    seqs_int = np.array([
        [ord(x) for x in letters],
        [ord(x) for x in letters]
    ])

    coords = np.array([
        [0, 4],
        [0, 4]
    ])

    print(f"Using np {seqs_char.dtype} array for letters:")
    print(Alignment(seqs_char, coords))

    print(f"Using np {seqs_int.dtype} array of ascii codes for 'ATGC':")
    print(Alignment(seqs_int, coords))

def example_aln():
    seqs = ["CCGGTTTTT", "AGTTTAA", "AGGTTT"]
    coords = np.array([
        [1, 3, 4, 7, 9],
        [0, 2, 2, 5, 5],
        [0, 2, 3, 6, 6]
    ])
    return Alignment(seqs, coords)


def slicing():
    """
    Shows slicing of alignments.
    """
    aln = example_aln()
    print(aln)

    print("First row:")
    print(aln[0])

    print("First row, columns 1, 2, 4:")
    # Use an iterable to slice columns!
    print(aln[0, [1, 2, 4]])

    print("First column:")
    print(aln[:, 0])

    print("All rows, columns 0 to 2:")
    print(aln[:, 0:2])


def aln_shapes():
    """
    Shows shape information of the alignment.
    """
    aln = example_aln()
    print("Number of sequences in alignment:")
    print(len(aln))

    print("Number of matches+mismatches+len_gaps (so bases total):")
    print(aln.length)

    print("Shape of alignment:")
    print(aln.shape)


def comparing_alns():
    """
    Note alignments are equal to each other if and only if they have the same
    sequences and coordinates. Otherwise this is false.
    """
    aln = example_aln()

    print("Same aln?")
    print(aln == aln)

    # Now for slightly longer aln but still correct.
    seqs = ["CCGGTTTTTG", "AGTTTAAC", "AGGTTTA"]
    coords = np.array([
        [1, 3, 4, 7, 9],
        [0, 2, 2, 5, 5],
        [0, 2, 3, 6, 6]
    ])
    aln2 = Alignment(seqs, coords)

    print("Same aln with just one extra base on each seq?")
    print(aln == aln2)


def aligned():
    """
    Shows the start and end indicies in target/query sequence as an np array.
    This only works for TWO sequences.
    Format is [[target/subject pairs], [query pairs]]
    """
    # Make pairwise alignment (just 2 seqs.)
    pair_aln = example_aln()[:2, :]
    print("Pairwise alignment:")
    print(pair_aln)

    print("Pairwise_alignment.aligned")
    print(pair_aln.aligned)

    print("Annotated:")
    print("Subject:")
    print(pair_aln.aligned[0])
    print("Query:")
    print(pair_aln.aligned[1])


def other_alignment_attrs():
    """
    Shows alignment attributes such as counts of identities, mismatches, gaps,
    letter frequencies, and substitutions.
    """
    aln = example_aln()
    aln.map()

    print("Alignment:")
    print(aln)

    print("Alignment.counts:")
    print(aln.counts())

    print("Alignment.frequencies:")
    print(aln.frequencies)

    print("Alignment.substitutions:")
    print("Format is 'row' that are aligned to 'col'")
    print(aln.substitutions)


def aln_as_arr():
    """
    You can turn an alignment into an np array of letters. Default type is |S1,
    being bytes, but you can do 'U' for unicode. This lets you do calculations
    faster than with the sequence objects.
    """
    aln = np.array(example_aln(), dtype="U")

    print("Aln as np.arr:")
    print(aln)


def aln_ops():
    """
    Alignment operations. You can sort how you please, and make rcomps.
    You can also add alignments given that they are of the same number of seqs.
    Why you would add them I have no idea, because they don't have to be related
    to add them.
    """
    aln = example_aln()

    print("Alignment:")
    print(aln)

    print("Alignment sorted by gc:")
    aln_sorted = example_aln()
    aln_sorted.sort(key=gc_fraction)
    print(aln_sorted)

    # Annotations are removed, column numbers are kept (in reverse).
    print("Alignment rcomps:")
    print(aln.reverse_complement())


    # Common annotations are kept, others are lost.
    print("Added Alignments:")
    print(aln + aln)


def aln_mapping():
    """
    Alignment.map() allows you to align the query sequence of one PAIR alignment
    to the subject of another PAIR alignment. This is useful for evaluating
    multiple assemblies.
    """
    # See the following example of chromosome, transcript, and rnaseq seq
    # where the rnaseq is the final (spliced) version of the transcript.
    chr = "AAAAAAAACCCCCCCAAAAAAAAAAAGGGGGGAAAAAAAA"
    transcript = "CCCCCCCGGGGGG"
    transcript_coords = np.array([
        [8, 15, 26, 32],
        [0, 7, 7, 13]
    ])

    transcript_aln = Alignment([chr, transcript], transcript_coords)

    print("Transcript alignment:")
    print(transcript_aln)

    rnaseq = "CCCCGGGG"
    rna_coords = np.array([
        [3, 11],
        [0, 8]
    ])
    rna_aln = Alignment([transcript, rnaseq], rna_coords)
    print("rnaseq alignment:")
    print(rna_aln)

    print("rnaseq aligned to chr with .map():")
    print(transcript_aln.map(rna_aln))

def alignment_collection():
    """
    Alignments can be stored in the Alignments class, which is an iterator.
    It's better than a list because it can hold arbitrary properties.
    Use Alignments.rewind() to reset the iterator.
    """
    alns = Alignments([example_aln(), example_aln(), example_aln()])

    print("Alignments:")
    print(alns)
    for aln in alns:
        print(aln)

    alns.rewind()

    print("Some arbitrary attrs I'm setting:")
    alns.name = "wow"
    alns.origin = "my brain"
    alns.endianness = "little"
    print(".name: ", alns.name)
    print(".origin: ", alns.origin)
    print(".endianness: ", alns.endianness)


def reading_and_writing_alns():
    """
    Use Bio.Align.parse for reading multiple, .read for one. This is like SeqIO.
    Use Bio.Align.write for any number of alignments.
    """
    alns = Align.parse("opuntia.aln", "clustal")
    print("Alignments found: ", len(alns))
    
    aln = next(alns)

    print(aln)

    print("Note that was a clustal formatted alignment.")
    # format() is a method implemented by a class with the format_spec method.
    print("We can use format(aln, <fmt>) to make it maf.")
    print(format(aln, "maf"))

    print("Writing to a file (stringio):")
    s = StringIO()
    Align.write(aln, s, "maf")
    s.seek(0)
    print(s.read())
    print("Worked!")


def alignment_convert():
    """
    We can use Align.read and Align.write to easily convert between types.
    YOU MUST put them in an Alignments first
    """
    alns = Align.parse("opuntia.aln", "clustal")
    s = StringIO()
    print("Watch this fails.")
    try:
        Align.write(alns, s, "maf")
        return
    except ValueError:
        print("Fail on just Align.write(Align.parse())")
    finally:
        alns.rewind()


    input("Now enter to see Alignments work:")
    alns2 = Alignments(alns)
    Align.write(alns2, s, "maf")
    s.seek(0)

    # NOTE AlignIO fixes this complication easily to work as expected.

    print(s.read())



def main():
    manual_alignment()
    infer_coordinates()
    np_arr_as_sequences()
    slicing()
    aln_shapes()
    comparing_alns()
    aligned()
    other_alignment_attrs()
    aln_as_arr()
    aln_ops()
    aln_mapping()
    alignment_collection()
    reading_and_writing_alns()
    alignment_convert()

    # And stopped reading on alignment file formats because that can be
    # done later.

if __name__ == "__main__":
    main()



