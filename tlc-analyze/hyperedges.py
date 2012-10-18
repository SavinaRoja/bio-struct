#! /usr/bin/python

from Bio.PDB.PDBParser import PDBParser
import math
import time


parser = PDBParser()

PROXIMITY = 4.6  # Threshold for residues contacting, must be at least this close
CONSTRICT = 35  # Threshold for residue contact sanity check

#Default backbone atom definitions for the different biopolymers
#Most crystal structures do not include hydrogen atoms unlike NMR structures,
#it is wise in this case to either adjust the cutoff parameter accordingly or
#to exclude hydrogens.
PROTEIN_BACKBONE = ['N', 'CA', 'C', 'O']
DNA_BACKBONE = ['P', 'OP1', 'OP2', 'O5\'', 'C5\'', 'C4\'', 'O4\'', 'C3\'',
                'O3\'', 'C2\'','C1\'', 'H5\'', 'H5\'\'', 'H4\'', 'H3\'',
                'H2\'', 'H2\'\'', 'H1\'', 'HO5\'']  # HO5' occurs at the ends
RNA_BACKBONE = []

def find_mutuals(residue, other_residues=[], proximity=PROXIMITY,
                 constrict=CONSTRICT, exclude_backbone=True,
                 exclude_hydrogens=True):
    '''
    This provides the recursive core for finding mutual residue contacts. It
    will return a list of all contact paths it found. These represent all
    possible hyperedges, including many which are subsets of others. In
    typical use, the highest order unique hyperedges are of interest, and the
    results of this function will be slimmed down to those.
    
    This function detects residue contacts between the residue argument and all
    members of the other_residues argument. exclude_backbone is an optional
    argument, which may use defaults provided in this module when set to True,
    or may employ a custom atom backbone definition when supplied with such a
    list. If set to False, no backbone atoms will be excluded.
    
    exclude_hydrogens will cause the script to overlook all hydrogen atoms,
    these are most often found in NMR structures, 
    '''

    contacts = []

    for other in other_residues:  # Check all other residues
        contact = False
        for resi_atom in residue:  # Iterate over all residue atoms
            if resi_atom.get_name() in backbone:  # Skip backbone atoms
                continue
            for other_atom in other:  # Iterate over all other residue atoms
                if other_atom.get_name() in backbone:  # Skip backbone atoms
                    continue
                dist = atomic_distance(resi_atom, other_atom)
                #A sanity check for residue contacts, no pair of atoms should
                #be this far apart in touching protein residues.
                if dist >= CONSTRICT:
                    break
                elif dist <= CUTOFF:
                    contact = True
                    break
            if contact:
                break
        if contact:
            contacts.append(other)
            #contacts.append(other.get_id()[1])
    if not contacts:
        return [[residue]]
    elif len(contacts) == 1:
        return [[residue, contacts[0]]]
    else:
        new_contacts = []
        i = 0
        while i < len(contacts):
            new_contacts += find_mutuals(contacts[i],contacts[i + 1:])
            i += 1
        contacts = []
        for item in new_contacts:
            contacts.append([residue] + item)
        
        return contacts

def atomic_distance(atom_one, atom_two):
    '''
    Calculates the distance between two PDB Atom objects.
    '''
    #Retrieve atom coordinates
    point_one = atom_one.get_coord()
    point_two = atom_two.get_coord()
    #Calculate differences
    point_diffs = []
    for i in xrange(3):
        point_diffs.append(point_two[i] - point_one[i])
    #Sum the square of differences
    diff_sum = 0
    for num in point_diffs:
        diff_sum += math.pow(num, 2)
    #Return the square root of the sum of squared differences
    return math.sqrt(diff_sum)
    

def get_hyperedges(pdbfile, modelnumber=0, chain_keys=[], no_hetero=True,
                   proximity=PROXIMITY, constrict=CONSTRICT,
                   exclude_backbone=True, exclude_hydrogens=True):
    '''
    This function operates on a structure to compile a list of hyperedges. A
    hyperedge is defined as the mutual contact between two or more residues.
    This code should function for all biopolymer structures. Selection of
    residues within a structure may be prescribed with a single model number
    and a list of chains.
    '''
    
    start = time.time()
    structure = parser.get_structure('test', pdbfile)
    
    if len(structure) == 1:
        print('Structure file contains 1 model.')
    else:
        print('Structure file contains {0} models.'.format(len(structure)))
    
    model = structure[modelnumber]
    chains=[]
    
    if not chain_keys:  # Assume all chains if none are given
        chains = model.get_list()
    else:  # Only use the chains as listed in the chain_keys argument
        for c in chain_keys:
            chains.append(model[c])
    
    #Compile all residues into a single list
    residues = []
    for c in chains:
        for resi in c:
            if no_hetero:
                if resi.get_id()[0].strip():
                    continue
            residues.append(resi)
    
    hyperedges = set()
    #Find all the mutually contacting residue side groups (hyperedge)
    while residues:
        resi = residues.pop(0)
        for edge in find_mutuals(resi, residues, proximity=proximity,
                                 constrict=constrict,
                                 exclude_backbone=exclude_backbone,
                                 exclude_hydrogens=exclude_hydrogens):
            if len(edge) > 1:
                res_numbers = [i.get_id()[1] for i in edge]
                hyperedge = frozenset(res_numbers)
                is_subset = False
                subbed = False
                for he in hyperedges:
                    if hyperedge.issubset(he):
                        is_subset = True
                        break
                    if hyperedge.issuperset(he):
                        subbed = he
                if subbed:
                    hyperedges.remove(subbed)
                if not is_subset:
                    hyperedges.add(hyperedge)
    print('Hyperedges calculated in {0} seconds'.format(time.time() - start))
    return hyperedges
