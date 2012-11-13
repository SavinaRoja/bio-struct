'''
A class-based approach to my protein hypergraph code.
'''

from Bio.PDB.PDBParser import PDBParser
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
    def __init__(self, structure_file_path, modelnumber, chain_keys):
        self.pdb_structure = parser.get_structure('test', structure_file_path)
        self.hyperedges = set()
        self.hyperresidues = set()

    
        
        
    def get_hyperedges(self, modelnumber=0, chain_keys=[], no_hetero=True,
                         proximity=4.6, constrict=35, exclude_backbone=True,
                         exclude_hydrogens=True):
        '''
        This function operates on a structure to compile a list of hyperedges.
        A hyperedge is defined as the mutual contact between two or more
        residues.
        '''
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
        return hyperedges
    
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