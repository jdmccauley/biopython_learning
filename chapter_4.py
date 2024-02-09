#!/usr/bin/env python

from io import StringIO
import requests

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import *

from Bio import SeqIO

def records():
    """
    A SeqRecord is a Seq with extra info like id, name, description, and features.
    Note the id is what's the name in a fasta file.
    """

    seq = SeqRecord(
        seq=Seq("ATGATGTAA"),
        id="1",
        name="My seq",
        description="An example seq.",
        features=[
            SeqFeature(
                SimpleLocation(0, 9, strand=1),
                type="CDS"
            )
        ],
        annotations={"source": "I made it up."}
    )
    seq.annotations["test_numbers"] = [1, 2, 3, 4]

    print("seq:")
    print(seq)
    for feature in seq.features:
        print(feature)


def fastas():
    """
    Using fastas for a SeqRecord.
    Note it has no annotations.
    A fasta is formatted like this:
        >(id) (description)
        (sequence)
    NOTE THE SPACE BETWEEN THE ID AND DESCRIPTION
    """
    dna = StringIO(
        requests.get(
            "https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.fna"
        ).text
    )

    rec = SeqIO.read(dna, "fasta")

    print("Read fasta:")
    print(rec)

def gbs():
    dna = StringIO(requests.get(
        "https://raw.githubusercontent.com/biopython/biopython/master/Tests/GenBank/NC_005816.gb"
    ).text)

    rec = SeqIO.read(dna, "genbank")

    print("Read gb:")
    print(rec)
    print("Annotations:")
    for a in rec.annotations:
        print(a)
    print("Features:")
    print(len(rec.features))

def fuzzy_features():
    rec = SeqRecord(
        Seq("ATGGCGTAAATGTAA"),
        id="My seq",
        name="Sample text"
    )

    start_codon = SeqFeature(SimpleLocation(
        BeforePosition(2), AfterPosition(2)
    ), type="CDS", qualifiers={
        "note": ["Start codon"]
    })

    rec.features.append(start_codon)

    stop_codon = SeqFeature(SimpleLocation(
        WithinPosition(6, 6, 9), OneOfPosition(8, [8, 9, 10])
    ), type="CDS", qualifiers={
        "note": ["Stop codon"]
    })

    rec.features.append(stop_codon)

    print("Added features:")
    for feature in rec.features:
        print(feature)


def location_testing():
    """
    You can probe for a location in all features.
    """
    rec = SeqRecord(
        Seq("ATGATGGCCACCTAA"),
        features=[
            SeqFeature(SimpleLocation(0, 3, 1),type="CDS"),
            SeqFeature(SimpleLocation(3, 6, 1),type="CDS"),
            SeqFeature(SimpleLocation(6, 9, 1),type="CDS")
        ]
    )

    snp_1 = 10
    snp_2 = 3

    for feature in rec.features:
        if snp_1 in feature:
            print("1")
        if snp_2 in feature:
            print("2")


def example_seqrec():
    return SeqRecord(
        Seq("ATGATGGCCACCTAA"),
        features=[
            SeqFeature(SimpleLocation(0, 3, 1),type="CDS"),
            SeqFeature(SimpleLocation(3, 6, 1),type="CDS"),
            SeqFeature(SimpleLocation(6, 9, 1),type="CDS")
        ],
        # This is required for genbank for format correctly.
        # Can be genomic, cDNA...
        annotations={
            "molecule_type": ["dna"],
            "note": ["This is a test"],
            "topology": ["circular"]
        },
        id="samplerec",
        description="A sample record.",
    )


def feature_seq():
    """
    Features themselves don't have sequences, but they can be used to extract
    slices from any Seq or SeqRecord object.
    """
    rec = example_seqrec()

    feat_1: SeqFeature = rec.features[0]

    print("Feature 1 extracted with extract()")
    print(feat_1.extract(rec.seq))

    print("Feature 1 with manual slicing")
    print(rec.seq[
        rec.features[0].location.start:rec.features[0].location.end
    ])

    print("Feature 1 extraced from another seq:")
    seq2 = Seq("DYKK")
    print(feat_1.extract(seq2))


def comparing():
    """
    Comparing SeqRecords is intentionally not simple. You can't compare SeqRecords
    themselves, only the attributes.
    """
    rec = example_seqrec()

    rec2 = rec

    rec3 = example_seqrec()

    print("rec == rec2?")
    try:
        print(rec == rec2)
    except NotImplementedError:
        print("Not implemented!")

    print("rec == rec3")
    try:
        print(rec == rec3)
    except NotImplementedError:
        print("Not implemented!")

    print("rec.seq == rec2.seq?")
    print(rec.seq == rec2.seq)

    print("rec.seq == rec3.seq?")
    print(rec.seq == rec3.seq)

def formatting():
    """
    You can print strings formatting the SeqRecord as you please.
    """
    rec = example_seqrec()
    print("fasta:")
    print(rec.format("fasta"))

    print("genbank:")
    print(rec.format("genbank"))


def slicing():
    """
    Note that slicing does the following:
     * keeps features, but adjusts their locations
     * removes many annotations that may not apply, like toplology and notes
     * id and description are kept

    Combining slices:
     * similarly keeps features but removes annototations
     * id and description kept

    To keep annotations: new_rec = rec.annotations.copy()
    """
    rec = example_seqrec()[6:9]
    print("Slice of rec 6:9")
    print(rec)
    print(rec.features)
    print(rec.format("genbank"))


    print("Rearranged rec by slicing:")
    left = example_seqrec()[:5]
    right = example_seqrec()[6:]
    new_rec = example_seqrec()[:5] + example_seqrec()[6:]

    print(new_rec)
    print(new_rec.features)
    print(new_rec.format("genbank"))


def rcomp():
    """
    SeqRecords can be reverse complemented, but the id+description+annotations
    are dropped because they're no longer relevant and could mess up tracking. 
    Features and per letter annotations (quality) are kept for obvious reasons.
    """
    rec = example_seqrec()
    rev = rec.reverse_complement()

    print("Reverse seq:")
    print(rev)


def main():
    records()
    fastas()
    gbs()
    fuzzy_features()
    location_testing()
    feature_seq()
    comparing()
    formatting()
    slicing()
    rcomp()

if __name__ == "__main__":
    main()
