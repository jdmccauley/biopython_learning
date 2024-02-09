#!/usr/bin/env python

import os

from io import StringIO
from tempfile import TemporaryFile
from copy import deepcopy

import subprocess as sp

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree

from Bio import AlignIO


def tree_intro():
    """
    See a tree.
    """
    tree = Phylo.read("simple.dnd", "newick")

    Phylo.draw(tree)


def color_and_width():
    """
    You can add color and width to clade branches, but it's only saved in
    PhyloXML files, not simple dnd files.
    """
    tree: Tree = Phylo.read("simple.dnd", "newick")

    # Set root
    tree.rooted = True
    tree.root.color = "gray"

    ancestor = tree.common_ancestor(
        {"name": "E"},
        {"name": "F"}
    )
    ancestor.color = "salmon"

    # Or second child of root's first child:
    grandchild = tree.clade[0, 1]
    grandchild.color = "blue"

    Phylo.draw(tree)

    # Now write to file.
    s = StringIO()
    Phylo.write(tree, s, "phyloxml")
    s.seek(0)
    print(s.read())


def io_phylo():
    """"
    Phylo IO works like SeqIO and AlignIO, but in the base package alone!
    """
    one_tree = Phylo.read("int_node_labels.nwk", "newick")
    print(Phylo.draw_ascii(one_tree))

    print("Writing file:")
    with TemporaryFile() as f:
        x = Phylo.write(one_tree, f, "newick")
        if x:
            print("Done")


    many_tree = [
        tree for tree in 
        Phylo.parse("phyloxml_examples.xml", "phyloxml")
    ]

    print("Many trees number: ", len(many_tree))

    print("Converting:")

    with TemporaryFile() as f:
        v = Phylo.convert("int_node_labels.nwk", "newick", "tmp.xml", "phyloxml")
        if v:
            print("Done")

def more_annotations():
    """
    You can annotate branches with labels! Use branch_length to get branch
    length.
    """
    tree = Phylo.read("example.xml", "phyloxml")
    Phylo.draw(tree, branch_labels=lambda clade: clade.branch_length)


def tree_traversal():
    """
    Shows how to traverse a tree.
    """
    tree = Phylo.read("simple.dnd", "newick")
    Phylo.draw_ascii(tree)
    
    print("\nDFS traversal:")
    # Can use regex for attribute values.
    nodes = tree.find_clades(
       {"name": ".*"},
       order="preorder"
    )
    for node in nodes:
        print(node.name)

    print("\nBFS traversal:")
    nodes = tree.find_clades(
       {"name": ".*"},
       order="level"
    )
    for node in nodes:
        print(node.name)

    print("Levels:")
    print(tree.depths())



def modifications():
    """
    You can modify trees, do copy them first.
    """

    tree = Phylo.read("my_simple.dnd", "newick")
    print("tree:")
    Phylo.draw_ascii(tree)

    # Now to do some work.
    newtree = deepcopy(tree)

    print("B branch length: ", newtree.find_any("B").branch_length)

    print("Removing the leaf A with prune:")
    newtree.prune("A")
    Phylo.draw_ascii(newtree)

    print("B branch length: ", newtree.find_any("B").branch_length)
    print("See that B had its length changed since ancestor AB is no longer relevant.")

    print("Deleting CD with collapse:")
    newtree.collapse("CD")
    Phylo.draw_ascii(newtree)

    print("Deleting F with collapse:")
    newtree.collapse("F")
    Phylo.draw_ascii(newtree)

    print("Deleting EFG with collapse:")
    newtree.collapse("EFG")
    Phylo.draw_ascii(newtree)


    # Use collapse generally.

    print("collapse_all:")
    newtree.collapse_all()
    Phylo.draw_ascii(newtree)


def changing_rooting():
    """
    You can change how trees are rooted!
    """
    tree = Phylo.read("my_simple.dnd", "newick")

    mid = deepcopy(tree)
    mid.root_at_midpoint()

    print("original:")
    Phylo.draw_ascii(tree)

    print("midpoint rooted tree:")
    Phylo.draw_ascii(mid)

    # Now root at A. Just use identifier in the function.
    a_root = deepcopy(tree)
    a_root.root_with_outgroup({"name": "A"})
    
    print("Rooted at A ancestor")
    Phylo.draw_ascii(a_root)

    print("original:")
    print(tree)

    print("rooted a:")
    print(a_root)


def using_alignment():
    """
    You can use an alignment as an input to an external application to do
    alignments.
    """
    # phyml uses phylip input.
    success = AlignIO.convert("opuntia.aln", "clustal", "opuntia.phy", "phylip")

    if not success:
        raise ValueError
    res = sp.run(
        "phyml -i opuntia.phy".split(),
        capture_output=True
    )

    if res.returncode != 0:
        print(res.stderr)
        raise ValueError
    
    tree = Phylo.read("opuntia.phy_phyml_tree.txt", "newick")
    Phylo.draw_ascii(tree)


def main():
    tree_intro()
    color_and_width()
    io_phylo()
    more_annotations()
    tree_traversal()
    modifications()
    changing_rooting()
    using_alignment()


if __name__ == "__main__":
    main()

