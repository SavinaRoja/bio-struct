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
    
    def hyper_analyze(self, unaligned, aligned, offset, family_msa):
        '''
        A sort of default workflow routine, though obviously all these steps
        could be called independently or modified. This at least a basic
        example of the work that needs to be done.
        '''
        #Parse the structure file to compile a set of all hyperedges
        self.get_hyperedges()
        #Determine how to align the structure sequence to the MSA
        self.unaligned_parent = str(Bio.SeqIO.read(unaligned, 'fasta').seq)
        self.aligned_parent = str(Bio.SeqIO.read(aligned, 'fasta').seq)
        self.get_alignment_map(offset)
        #Align the hyperedges set, using the alignment map
        self.get_aligned_hyperedges()
        #Get all hyperresidues for all hyperedges, as represented by the
        #sequences in the MSA
        self.get_hyperresidues(family_msa)
        #Calculate the residue-wise "q" and "phi" scores
        self.calc_residue_scores()
        #Calculate the edge-wise score: "edge weight"
        self.calc_edge_weights()
        #Calculate the sequence-wise scores, sum weight of connected edges
        self.calc_edge_weights_by_seq(offset)  # Produces both aligned and unaligned

    def dealign(self, sequence):
        return ''.join(sequence.split('-'))

    def calc_edge_weights_by_seq(self, offset):
        self.edge_weights_by_seq = [0] * (len(self.aligned_parent) + offset)
        for i in xrange(len(self.edge_weights_by_seq)):
            for edge in self.edge_weights:
                if i in edge:
                    self.edge_weights_by_seq[i] += self.edge_weights[edge]
        self.unaligned_edge_weights_by_seq = [0] * (len(self.unaligned_parent) + offset)
        for i in xrange(len(self.unaligned_edge_weights_by_seq)):
            j = self.alignment_map[i]
            for edge in self.edge_weights:
                if j in edge:
                    self.unaligned_edge_weights_by_seq[i] += self.edge_weights[edge]

    def calc_edge_weights(self):
        self.edge_weights = {}
        for edge in self.phi_scored_hyperresidues:
            sum = 0
            for residue in self.phi_scored_hyperresidues[edge]:
                normalized_count = self.hyperresidues[edge][residue] / self.seq_num
                resi_phi = self.phi_scored_hyperresidues[edge][residue]
                sum += normalized_count * resi_phi
            self.edge_weights[edge] = sum

    def calc_residue_scores(self):
        #First calculate the q_scores, rho is infinite at the moment
        self.q_scored_hyperresidues = {}
        for edge in self.hyperresidues:
            q_scores = {}
            for residue in self.hyperresidues[edge]:
                q_scores[residue] = self.hyperresidues[edge][residue] / self.seq_num
            self.q_scored_hyperresidues[edge] = q_scores
        #Calculate the phi scores
        self.phi_scored_hyperresidues = {}
        for edge in self.hyperresidues:
            phi_scores = {}
            #May the gods forgive me for the code I am about to write
            for residue in self.hyperresidues[edge]:
                if len(residue) == 1:  # First order
                    phi_score = math.log10(self.q_scored_hyperresidues[edge][residue])
                    phi_scores[residue] = phi_score
                elif len(residue) == 2:  # Second order
                    #Numerator
                    num = self.q_scored_hyperresidues[edge][residue]
                    #I have to decompose the hyperedge to its lesser edges
                    edge1, edge2 = tuple([edge[0]]), tuple([edge[1]])
                    #Same for the residues. This DEPENDS on strict ordering in previous steps
                    resi1, resi2 = tuple([residue[0]]), tuple([residue[1]])
                    #Denominator
                    den = self.q_scored_hyperresidues[edge1][resi1] * self.q_scored_hyperresidues[edge2][resi2]
                    phi_score = math.log(num / den)
                    phi_scores[residue] = phi_score
                elif len(residue) == 3:  # Third order
                    #This third order residue
                    resi_q = self.q_scored_hyperresidues[edge][residue]
                    #I have to decompose the hyperedge to second and first order edges
                    edge11, edge12, edge13 = [i for i in edge]
                    edge21, edge22, edge23 = edge[:2], edge[0:3:2], edge[1:]
                    #Same for the residues. This DEPENDS on strict ordering in previous steps
                    resi11, resi12, resi13 = [i for i in residue]
                    resi21, resi22, resi23 = residue[:2], residue[0:3:2], residue[1:]
                    #Calculate the numerator
                    nq1 = self.q_scored_hyperresidues[tuple([edge11])][tuple([resi11])]
                    nq2 = self.q_scored_hyperresidues[tuple([edge12])][tuple([resi12])]
                    nq3 = self.q_scored_hyperresidues[tuple([edge13])][tuple([resi13])]
                    num = resi_q * nq1 * nq2 * nq3
                    #Calculate the denominator
                    dq1 = self.q_scored_hyperresidues[edge21][resi21]
                    dq2 = self.q_scored_hyperresidues[edge22][resi22]
                    dq3 = self.q_scored_hyperresidues[edge23][resi23]
                    den = dq1 * dq2 * dq3
                    #Finally, the phi score
                    phi_score = math.log10(num / den)
                    phi_scores[residue] = phi_score
                elif len(residue) == 4:  # Fourth order
                    #This is not a simplified algrebraic approach
                    #This fourth order residue
                    resi_q = self.q_scored_hyperresidues[edge][residue]
                    #Decompose to third order
                    edge31 = tuple([edge[0], edge[1], edge[2]])
                    edge32 = tuple([edge[0], edge[1], edge[3]])
                    edge33 = tuple([edge[0], edge[2], edge[3]])
                    edge34 = tuple([edge[1], edge[2], edge[3]])
                    edges3 = [edge31, edge32, edge33, edge34]
                    resi31 = tuple([residue[0], residue[1], residue[2]])
                    resi32 = tuple([residue[0], residue[1], residue[3]])
                    resi33 = tuple([residue[0], residue[2], residue[3]])
                    resi34 = tuple([residue[1], residue[2], residue[3]])
                    resis3 = [resi31, resi32, resi33, resi34]
                    qs3 = [self.q_scored_hyperresidues[i][j] for i,j in zip(edges3, resis3)]
                    #Decompose to second order
                    edge21 = tuple([edge[0], edge[1]])
                    edge22 = tuple([edge[0], edge[2]])
                    edge23 = tuple([edge[0], edge[3]])
                    edge24 = tuple([edge[1], edge[2]])
                    edge25 = tuple([edge[1], edge[3]])
                    edge26 = tuple([edge[2], edge[3]])
                    edges2 = [edge21, edge22, edge23, edge24, edge25, edge26]
                    resi21 = tuple([residue[0], residue[1]])
                    resi22 = tuple([residue[0], residue[2]])
                    resi23 = tuple([residue[0], residue[3]])
                    resi24 = tuple([residue[1], residue[2]])
                    resi25 = tuple([residue[1], residue[3]])
                    resi26 = tuple([residue[2], residue[3]])
                    resis2 = [resi21, resi22, resi23, resi24, resi25, resi26]
                    qs2 = [self.q_scored_hyperresidues[i][j] for i,j in zip(edges2, resis2)]
                    #Decompose to first order
                    edges1 = [tuple([edge[i]]) for i in xrange(4)]
                    resis1 = [tuple([residue[i]]) for i in xrange(4)]
                    qs1 = [self.q_scored_hyperresidues[i][j] for i,j in zip(edges1, resis1)]
                    #Construct our denominator
                    all_qs = qs3 + qs2 + qs1
                    den = 1
                    for q in all_qs:
                        den = den * q
                    phi_score = math.log10(resi_q / den)
                    phi_scores[residue] = phi_score
                else:  # Greater than 4... ewww
                    raise ValueError('Edge order greater than 4!')
            self.phi_scored_hyperresidues[edge] = phi_scores

    def get_hyperresidues(self, msa_file):
        '''
        This function calculates and returns the hyperresidues from a multiply-
        aligned family file based on the aligned hyperedges in the hypergraph.
        '''
        
        resi_type = {'C': 1, 'F': 2, 'Y': 2, 'W': 2, 'H': 3, 'R': 3, 'K': 3,
                     'N': 4, 'D': 4, 'Q': 4, 'E': 4, 'S': 5, 'T': 5, 'P': 5,
                     'A': 5, 'G': 5, 'M': 6, 'I': 6, 'L': 6, 'V': 6, '-': 0}
        
        #A dictionary keyed by hyperedges will hold values of hyperresidues
        hyperresidues = {}
        msa = Bio.AlignIO.read(msa_file, 'fasta')
        sequences = [al.seq.__str__() for al in msa]
        self.seq_num = float(len(sequences))
        for edge in self.aligned_hyperedges:
            edge_list = list(edge)
            edge_list.sort()
            edge_tuple = tuple(edge_list)
            hyperres_dict = {}
            for seq in sequences:
                vertices = tuple([resi_type[seq[i]] for i in edge_tuple])
                if 0 in vertices:  # This edge does not exist in this sequence
                    continue
                elif vertices in hyperres_dict:  # Seen before, add 1
                    hyperres_dict[vertices] += 1
                else:  # New, set to 1
                    hyperres_dict[vertices] = 1
            hyperresidues[edge_tuple] = hyperres_dict
        self.hyperresidues = hyperresidues
        
        
    def get_alignment_map(self, offset):
        '''
        Create an index map from the unaligned parent to the parent sequence
        in the family MSA. The offset value can be used to account for an
        index shift of the unaligned parent, 3EYC begins at 13 for instance.
        '''
        unaligned = self.unaligned_parent
        aligned = self.aligned_parent
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