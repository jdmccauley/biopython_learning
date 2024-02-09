#!/usr/bin/env python

from io import StringIO, BytesIO
import requests

from Bio.Seq import Seq, reverse_complement, transcribe, translate
from Bio.SeqUtils import gc_fraction

# Seq object is the focus of this chapter.
# It handles sequences only, while SeqIO handles annotations+ids+etc.

# Seq objects have methods for the central dogma,
# complement, revcomp, transcript, translate.

def seq_indexing():
    """
    Showing some Seq indexing.
    """
    seq = Seq("ATGCATGCATGC")

    print("Printing all bases...\n")
    for i, base in enumerate(seq):
        print(f"{i}\t{base}")

    print("Printing specific bases...\n")
    print(seq[0])
    print(seq[2])


    print("Counting patterns: ATGC")
    # Note that strings allow this as well: nonoverlapping counts.
    print(seq.count("ATGC"))

    # Get gc fraction
    print(f"GC fraction: {gc_fraction(seq)}")

def slicing():
    """
    Showing some slicing.
    """
    seq = Seq("ATGCATGCATGC")

    print("every other 3rd base...")
    print(seq[::3])

    print("Reversed seq...")
    # Note this also works for regular strings.
    print(seq[::-1])


def concats():
    """
    Showing some concatenation.
    """
    dna = Seq("ATGC")
    protein = Seq("DYKDDDDK")

    # Note that Seq does not prevent DNA + protein sequences.
    print(dna + protein)

    dnas = [dna, dna, dna]
    # Join = iter[0], subject, iter[1], subject...
    print(protein.join(dnas))


def cases():
    """
    Upper and lower like strings.
    """
    seq = Seq("ATGCatgcAA")

    print("Upper:")
    print(seq.upper())

    print("Lower:")
    print(seq.lower())

    print("Is aa in seq?")
    print("aa" in seq)
    print("Is aa in seq.lower()?")
    print("aa" in seq.lower())


def complements():
    """
    Showing complements and reverse complements for DNA and proteins.
    """
    dna = Seq("ATGCATGC")
    protein = Seq("DYE")

    print("Seq:")
    print(dna)
    print("Complement:")
    print(dna.complement())
    print("rcomp:")
    print(dna.reverse_complement())

    print("Protein:")
    print(protein)
    # Note this works because some of the bases I put for the protein
    # are IUPAC nucleotides. Ex: D = A or G or T, and H = A or C or T.
    # So complement/rcomp makes the NUCLEOTIDE complement if the base is any nucleotide,
    # otherwise it does not change the base.
    print("Complement:")
    print(protein.complement())
    print("rcomp:")
    print(protein.reverse_complement())

def transcription():
    """
    Remember transcription happens on the 3' strand to make 5'->3'
    So 5'->3' is the coding strand, 3'->5' is the template strand!
    Note "sense strand" == "coding stand"
         "antisense strand" == "template strand"

    Ex:
      Seq:
     5' ATGCATGC 3'
     3' TACGTACG 5'

    Making the mRNA:
     5' AUGCAUGC 3'
     3' TACGTACG 5'

    So, given a coding region 5' to 3', just switch the T->U
    to get the mRNA.
    """
    print("Note all strands printed 5'->3'")

    coding_strand = Seq("ATGCATGC")
    template_stand: Seq = coding_strand.reverse_complement()

    print("Coding strand:")
    print(coding_strand)
    print("Template strand:")
    print(template_stand)

    print("Quick transcription using coding strand:")
    print(coding_strand.transcribe())

    print("Actual transcription with template strand:")
    print(template_stand.reverse_complement().transcribe())

    rna = coding_strand.transcribe()
    print("Back transcribing:")
    print(rna.back_transcribe())


def translation():
    """
    Translating from the coding stand and mRNA.
    """
    coding_strand = Seq("ATGATGGAATAACTG")

    # Can translate from RNA, and coding strand.
    print("From rna:")
    print(coding_strand.transcribe().translate())

    print("From coding strand:")
    print(coding_strand.translate())

    # You can change the stop symbol with stop_symbol=""
    # And you can say to stop at a stop codon with to_stop=True

    # You can use trandlation tables. Note that in bacteria, GTG
    # is a valid start codon. To use this, use the 'bacterial' table,
    # and set 'cds=True' to let the GTG be a M instead of a V like normal.

    b_strand = Seq("GTGGTGTAA")
    print("Default bacterial translation:")
    print(b_strand.translate(table="Bacterial"))
    print("With cds=True:")
    print(b_strand.translate(table="Bacterial", cds=True))

    # Note that cds=True throws an exception if the Seq is not one.

def translation_tables():
    """
    Here we talk about translation tables.
    """

    from Bio.Data import CodonTable

    standard = CodonTable.ambiguous_dna_by_name["Standard"]
    mito = CodonTable.ambiguous_dna_by_name["Vertebrate Mitochondrial"]

    print(standard)
    print(mito)

    print(standard.stop_codons)
    print(mito.start_codons)
    print(mito.forward_table["ACG"])


def comparing_seqs():
    """
    Seqs are compared as raw strings.
    """
    s = Seq("ATG")

    print(s == "ATG")


def partial_and_empty_seqs():
    """
    For sequences with no data or partial data.
    """
    from Bio.Seq import UndefinedSequenceError
    none_seq = Seq(None, 10)

    print("Noneseq len:")
    print(len(none_seq))

    print("Noneseq data:")
    try:
        print(none_seq)
    except UndefinedSequenceError:
        print("None")

    region = Seq({
        100: "ATG"
    }, length=10000)

    print("region at 0:99")
    try:
        print(region[:99])
    except UndefinedSequenceError:
        print("None")

    print("region 100:103")
    print(region[100:103])

    new_seq = none_seq + region

    print("noneseq concat with region, at 110:113")
    print(new_seq[110:113])

def mutable_seq():
    """
    Seq is immutable by default, and can thus be used as a dictionary key.
    To modify Seq, make it a MutableSeq. Note it can then no longer be a dictionary key.
    Also MutableSeq.reverse() reverses in place!
    """
    from Bio.Seq import MutableSeq
    seq = Seq("ATGC")
    mut = MutableSeq(seq)

    print("seq, mut:")
    print(seq, mut)

    mut[0] = 'T'

    print("Modified mut:")
    print(seq, mut)

    print("reversed mut:")
    mut.reverse()
    print(mut)


def just_strings():
    """
    You can do any Seq method on a string.
    """

    seq = "ATGATG"
    print("str:")
    print(seq)

    print("rcomp:")
    print(reverse_complement(seq))


    print("transcribed:")
    print(transcribe(seq))
    print("translated:")
    print(translate(seq))

def main():
    seq_indexing()
    slicing()
    concats()
    cases()
    complements()
    transcription()
    translation()
    translation_tables()
    comparing_seqs()
    partial_and_empty_seqs()
    mutable_seq()
    just_strings()

if __name__ == "__main__":
    main()

    
