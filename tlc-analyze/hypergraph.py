'''
A class-based approach to my protein hypergraph code.
'''

from Bio.PDB.PDBParser import PDBParser
import Bio.SeqIO
import Bio.AlignIO

import math

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

class Hypergraph(object):
    '''
    A class for my protein hypergraph, presents a coherent united framework for
    the hypergraph analysis.
    '''
    def __init__(self, structure_file_path):
        self.pdb_structure = parser.get_structure('test', structure_file_path)
        self.hyperedges = set()
        self.hyperresidues = set()
    
    def hyper_analyze(self):
        self.get_hyperedges()
        self.unaligned_parent = str(Bio.SeqIO.read('./data/3eyc_a_from_pdb.fa', 'fasta').seq)
        self.aligned_parent = str(Bio.SeqIO.read('./data/3eyc_parent.fa', 'fasta').seq)
        self.get_alignment_map(self.unaligned_parent, self.aligned_parent, offset=13)
        self.get_aligned_hyperedges()
        
    def get_alignment_map(self, unaligned, aligned, offset):
        '''
        Create an index map from the unaligned parent to the parent sequence
        in the family MSA. The offset value can be used to account for an
        index shift of the unaligned parent, 3EYC begins at 13 for instance.
        '''
        map = [0] * (offset + len(unaligned))
        i = 0
        j = 0
        while i < len(unaligned):
            una = unaligned[i]
            while una != aligned[j]:
                j += 1
            map[i + offset] = j + 1
            j += 1
            i += 1
        self.alignment_map = map
        

    def get_aligned_hyperedges(self):
        aligned_hyperedges = set()
        for he in self.hyperedges:
            hyper_list = []
            for node in he:
                new_index = self.alignment_map[node]
                hyper_list.append(new_index)
            aligned_hyperedges.add(frozenset(hyper_list))
        self.aligned_hyperedges = aligned_hyperedges

    
    def get_hyperedges(self, modelnumber=0, chain_keys=['A'], no_hetero=True,
                         proximity=4.6, constrict=35, exclude_backbone=True,
                         exclude_hydrogens=True):
        '''
        This function operates on a structure to compile a list of hyperedges.
        A hyperedge is defined as the mutual contact between two or more
        residues.
        '''
        model = self.pdb_structure[modelnumber]
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
            for edge in self.find_mutuals(resi, residues, proximity=proximity,
                                          constrict=constrict,
                                          exclude_backbone=exclude_backbone,
                                          exclude_hydrogens=exclude_hydrogens):
                if len(edge) > 0:
                    res_numbers = [i.get_id()[1] for i in edge]
                    hyperedge = frozenset(res_numbers)
                    hyperedges.add(hyperedge)
        self.hyperedges = hyperedges

    def find_mutuals(self, residue, other_residues=[], proximity=4.6,
                 constrict=35, exclude_backbone=True,
                 exclude_hydrogens=True):
        '''
        This provides the recursive core for finding mutual residue contacts. It
        will return a list of all contact paths it found. These represent all
        possible hyperedges of all orders.
        
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
            for resi_atom in residue:  # Iterate over all residue atoms
                makes_contact = False
                sanity = False
                if exclude_hydrogens and resi_atom.get_name()[0] == 'H':
                    continue
                if exclude_backbone and resi_atom.get_name() in PROTEIN_BACKBONE:
                    continue
                for other_atom in other:  # Iterate over all other residue atoms
                    if exclude_hydrogens and resi_atom.get_name()[0] == 'H':
                        continue
                    if exclude_backbone and resi_atom.get_name() in PROTEIN_BACKBONE:
                        continue
                    dist = self.atomic_distance(resi_atom, other_atom)
                    if dist >= constrict:  # Distance sanity check
                        sanity = True
                        break
                    elif dist <= proximity:
                        makes_contact = True
                        break
                if any([makes_contact, sanity]):
                    break
            if makes_contact:
                contacts.append(other)
        if not contacts:
            return [[residue]]
        elif len(contacts) == 1:
            return [[residue], [residue, contacts[0]]]
        else:
            new_contacts = []
            i = 0
            while i < len(contacts):
                new_contacts += self.find_mutuals(contacts[i],contacts[i + 1:])
                i += 1
            contacts = [[residue]]
            for item in new_contacts:
                contacts.append([residue] + item)
                #print(contacts)
                #raw_input()
            return contacts

    def atomic_distance(self, atom1, atom2):
        '''
        Calculates the distance between two PDB Atom objects.
        '''
        #Retrieve atom coordinates
        p1 = atom1.get_coord()
        p2 = atom2.get_coord()
        #Calculate differences
        p_diffs = [p1[i] - p2[i] for i in xrange(3)]
        #Sum the square of differences
        diff_sum = sum(i**2 for i in p_diffs)
        #Return the square root of the sum of squared differences
        return math.sqrt(diff_sum)