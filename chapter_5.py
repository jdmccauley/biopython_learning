#!/usr/bin/env python

from io import StringIO, BytesIO
import gzip
import bz2
import tempfile
import shutil

import requests

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

"""
SeqIO provides simple interfaces for working with sequence files.
Note it's importing them as SeqRecords!
"""

EXAMPLES_SLUG="https://raw.githubusercontent.com/biopython/biopython/master/Doc/examples"


def print_rec_info(rec: SeqRecord):
    """
    Writing this to just summarize a record: id, description, length, 
    and len features.
    """
    print(f"\nID:\t{rec.id}")
    print(f"Description:\t{rec.description}")
    print(f"Length:\t{len(rec)}")
    print(f"Number features:\t{len(rec.features)}\n")


def reading():
    """
    Can read one file with SeqIO.read().
    Read multiple with SeqIO.parse(), which returns an iterator.
    """
    multiple = StringIO(requests.get(f"{EXAMPLES_SLUG}/ls_orchid.fasta").text)


    many_rec = [rec for rec in SeqIO.parse(multiple, "fasta")]

    print(f"Records in multiple: {len(many_rec)}")
    print("First record in the mutliple:")
    print_rec_info(many_rec[0])

    print("Reading one and only one record by just next(SeqIO.parse()):")
    # Remember to seek(0) for already-read StringIO/BytesIO
    multiple.seek(0)
    print_rec_info(next(SeqIO.parse(multiple, "fasta")))


def reading_compressed():
    """
    You can read compressed files with SeqIO too, assuming you use a compression
    library to turn the compressed file into a handle.
    """
    # NOTE requests.get().text assumes plaintext content and guesses encoding.
    # When bytes, use requests.get().content.
    gz = BytesIO(
        requests.get(f"{EXAMPLES_SLUG}/ls_orchid.gbk.gz").content
    )

    bz = BytesIO(
        requests.get(f"{EXAMPLES_SLUG}/ls_orchid.gbk.bz2").content
    )

    # genbanks must be opened in text mode
    # use rb if file is in binary (not the encoding, the actual content)
    with gzip.open(gz, 'rt') as f:
        print_rec_info(next(SeqIO.parse(f, "genbank")))

    # genbanks must be opened in text mode
    with bz2.open(bz, 'rt') as f:
        print_rec_info(next(SeqIO.parse(f, "genbank")))


# SKIPPING 5.3 for now because I don't anticipate needing to get sequences
# from Entrez or SwissProt online databases.
        

def indexing_dict():
    """
    For creating in-memory or on-disk databases of sequences!
    """
    multiple_fastas = StringIO(
        requests.get(f"{EXAMPLES_SLUG}/ls_orchid.fasta").text
    )

    # As dictionary: the 'in-memory' database.
    # This method does check for duplicate keys, so it's nice to use over
    # constucting the db ourselves.
    # The key is the id.
    dict_db = SeqIO.to_dict(
        SeqIO.parse(multiple_fastas, "fasta")
    )
    multiple_fastas.seek(0)

    # Note dict keys/values are not iterables, but they can be made into them.
    some_key = next(iter(dict_db.keys()))

    print("As dict, first key and value ('next'):")
    print(f"Key:\t{some_key}")
    print(f"Value:\t{dict_db[some_key]}")


    print("\nNow to use another way to create an id: with our own function!\n")

    def get_accension_key(rec: SeqRecord):
        return rec.id.split("|")[3]
    
    dict_db_acc = SeqIO.to_dict(
        SeqIO.parse(multiple_fastas, "fasta"),
        key_function=get_accension_key
    )
    multiple_fastas.seek(0)

    some_key = next(iter(dict_db_acc.keys()))

    print("As dict, first key and value ('next'):")
    print(f"Key:\t{some_key}")
    print(f"Value:\t{dict_db_acc[some_key]}")


    print("\nFinally, let's use someone else's function and a lambda.\n")

    from Bio.SeqUtils.CheckSum import seguid

    multiple_gbs = StringIO(
        requests.get(f"{EXAMPLES_SLUG}/ls_orchid.gbk").text
    )

    seguid_dict = SeqIO.to_dict(
        SeqIO.parse(multiple_gbs, "genbank"),
        lambda rec: seguid(rec.seq)
    )

    some_key = next(iter(seguid_dict.keys()))
    print(f"Length of the seguid_dict: {len(seguid_dict)}")
    print("As dict, first key and value ('next'):")
    print(f"Key:\t{some_key}")
    print(f"Value:\t{seguid_dict[some_key]}")


def indexing_file():
    """
    A dict db is large in memory, while indexing a file and reading as needed
    saves a lot of space. Note this must be done with a file, not a handle.

    The file can be fasta, gb, or bgz. Note bgz is handled automatically with
    no bgz module needed, but it does need to know if it's still a fasta or gb.
    """
    with open("ls_orchid.fasta", 'r') as f:
        id_dict = SeqIO.to_dict(SeqIO.parse(f, "fasta"))

    id_file = SeqIO.index("ls_orchid.fasta", "fasta")

    # # This also works:
    # id_file = SeqIO.index("ls_orchid.fasta.bgz", "fasta")

    some_key = next(iter(id_dict.keys()))

    from sys import getsizeof

    print(f"Size of dict: {getsizeof(id_dict)}")
    print(f"Size of indexed file: {getsizeof(id_file)}")

    # Note 3328 vs 48, that's 69 times the memory footprint,
    # at the cost of now dealing with file I/O.

    print(f"Value of {some_key} in dict:")
    print(id_dict[some_key])

    print(f"Value of {some_key} in indexed file:")
    print(id_file[some_key])


def indexing_db():
    """
    You can store records in a database too. A .idx file is made (a SQLite3 db)
    on disk. You can store records from multiple files here, given that each
    record identifier is unique.
    """
    temp_dir = tempfile.mkdtemp()

    print("Putting fastas in db...")
    db = SeqIO.index_db(
        f"{temp_dir}/ls_orchid.idx",
        ["ls_orchid.fasta"],
        "fasta", 
        # Note key function is for the raw record data, 
        # not SeqRecord since those aren't made. This also holds true for 
        # SeqIO.index. Only to_dict uses SeqRecords.
        key_function= lambda rec: rec.split("|")[3]
    )
    print("Done making db")

    print("Getting record Z78533.1:")
    print(db["Z78533.1"])

    shutil.rmtree(temp_dir)

def raw_data():
    """
    You can get the raw information from a record instead of the SeqRecord
    representation from get_raw(). NOTE THIS ONLY WORKS FOR SeqIO.index()
    and SeqIO.index_db()!
    """
    dict_db = SeqIO.index(
        "ls_orchid.fasta",
        "fasta",
        key_function= lambda rec: rec.split("|")[3]
    )

    print("The raw data for Z78533.1:")
    # Note that it's bytes, so we need to decode it.
    raw_bytes: bytes = dict_db.get_raw("Z78533.1")
    print(raw_bytes.decode("utf-8"))


def writing_files():
    """
    Bio.SeqIO.write is very flexible, it lets you take a list or iterator
    of SeqRecords. Note that any file written to is overwritten with no
    warning!
    """
    temp_dir = tempfile.TemporaryDirectory()
    print(temp_dir.name)
    newfile = f"{temp_dir.name}/new_ls_orchid.fasta"

    recs = [seq for seq in SeqIO.parse("ls_orchid.fasta", "fasta")]

    print("Original first record:")
    print(recs[0])

    SeqIO.write(
        recs, # or just SeqIO.parse(...)
        newfile,
        "fasta"
    )

    first_record = next(SeqIO.parse(newfile, "fasta"))
    print("New first record:")
    print(first_record)


def converting_files():
    """
    You can also just easily convert between filetypes.
    """

    tempdir = tempfile.TemporaryDirectory()
    newfile = f"{tempdir.name}/new_ls_orchid.fasta"

    SeqIO.convert("ls_orchid.gbk", "genbank", newfile, "fasta")

    print("Converted fasta record from gb:")
    print(next(SeqIO.parse(newfile, "fasta")))


    newrcomp = f"{tempdir.name}/rcomp_ls_orchid.fasta"

    # Note this makes a generator, not a list!
    # It saves a lot of memory!
    rcomps = (
        rec.reverse_complement()
        for rec in SeqIO.parse("ls_orchid.fasta", "fasta")
    )

    SeqIO.write(rcomps, newrcomp, "fasta")

    print("Converted fasta rcomp record from fasta:")
    print(next(SeqIO.parse(newrcomp, "fasta")))


def low_level_parsers():
    """
    You can use low level parsers as desired since those are called by
    SeqIO.parse(), and you can reduce the overhead if you call them
    directly.
    """
    from Bio.SeqIO.FastaIO import SimpleFastaParser

    count = 0
    total_len = 0
    with open("ls_orchid.fasta") as f:
        for title, seq in SimpleFastaParser(f):
            count += 1
            total_len += len(seq)

    print("Total records:")
    print(count)
    print("Total bases:")
    print(total_len)
        




def main():
    reading()
    reading_compressed()
    indexing_dict()
    indexing_file()
    indexing_db()
    raw_data()
    writing_files()
    converting_files()
    low_level_parsers()


if __name__ == "__main__":
    main()

