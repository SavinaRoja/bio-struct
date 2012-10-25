#! /usr/bin/python

"""
This section provides code for adjusting the hyperedges to the multiply aligned
sequence of the protein family. Conduct a multiple alignment, then extract a
copy of the appropriate parent sequence. This sequence should only differ from
your structure sequence by the inclusion of gaps.
"""

from Bio.PDB.PDBParser import PDBParser


parser = PDBParser()


def get_alignment_map(pdb_sequence, parent_sequence, initial_number):
    map = [0] * (initial_number + len(pdb_sequence))
    i = 0
    j = 0
    while i < len(pdb_sequence):
        pdb = pdb_sequence[i]
        while pdb != parent_sequence[j]:
            j += 1
        map[i + initial_number] = j + 1
        j += 1
        i += 1
    return map

def get_aligned_hyperedges(hyperedges, pdb_sequence, parent_sequence,
                             initial_number):
    alignment_map = get_alignment_map(pdb_sequence, parent_sequence,
                                      initial_number)
    aligned_hyperedges = set()
    for he in hyperedges:
        hyper_list = []
        for node in he:
            new_index = alignment_map[node]
            print(node, new_index)
            hyper_list.append(new_index)
        aligned_hyperedges.add(frozenset(hyper_list))
    return aligned_hyperedges

