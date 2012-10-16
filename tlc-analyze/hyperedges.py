#! /usr/bin/python

from Bio.PDB.PDBParser import PDBParser
import math


parser = PDBParser()

CUTOFF = 4.6  # Threshold for residues contacting, must be at least this close
CONSTRICT = 35  # Threshold for residue contact sanity check

backbone = ['N', 'CA', 'C', 'O']  # Protein backbone atom names

def find_mutuals(residue, other_residues):
    '''
    
    '''
    resseq = residue.get_id()[1]
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
            contacts.append(other.get_id()[1])
    
    if not contacts:
        return [resseq]
    else:
        new_contacts = []
        
        
        return new_contacts

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
    

def get_hyperedges(pdbfile, modelnumber=0, chain_keys=[], no_hetero=True):
    '''
    This function operates on a structure to compile a list of hyperedges. A
    hyperedge is defined as the mutual contact between two or more residues.
    This code should function for all biopolymer structures. Selection of
    residues within a structure may be prescribed with a single model number
    and a list of chains.
    '''
    
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
    
    hyperedges = []
    #Find all the mutually contacting residue side groups (hyperedge)
    while residues:
        resi = residues.pop(0)
        for hyperedge in find_mutuals(resi, residues):
            print(hyperedge)
    return hyperedges

###Development testing section###

get_hyperedges('./data/3EYC.pdb', chain_keys=['A'])
structure = parser.get_structure('3EYC', './data/3EYC.pdb')
model = structure[0]
chain_a = model['A']
residue = chain_a[13]
residue.get_id()[1]
ca = residue['CA']
n = residue['N']









