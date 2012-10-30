#! /usr/bin/python

from hyperedges import get_hyperedges
from hyperresidues import get_hyperresidues
import Bio.SeqIO
import Bio.AlignIO
import align


def main():
    #First, parse the pdb structure for the hyperedges
    with open('./processed/3EYC_hyperedges.txt', 'w') as out:
        hyperedges = get_hyperedges('./data/3EYC.pdb', chain_keys=['A'])
        for he in hyperedges:
            out.write(he.__str__())
            out.write('\n')

    #The hyperedges need to have their sequence numbers adjusted for the
    #multiple alignment. So we compare the original sequence to the aligned one.
    pdb_sequence = str(Bio.SeqIO.read('./data/3eyc_a_from_pdb.fa', 'fasta').seq)
    parent_sequence = str(Bio.SeqIO.read('./data/3eyc_parent.fa', 'fasta').seq)
    aligned_hyperedges = align.get_aligned_hyperedges(hyperedges,
                                                      pdb_sequence,
                                                      parent_sequence,
                                                      initial_number=13)

    #Now that the hyperedges are aligned to the multiple alignment, we can
    #compute the hyperresidues for a full picture of our hypergraph.
    hyperresidues = {}  # This data structure will be 'indexed' by the hyperedge frozenset tuples
    aligned_file = Bio.AlignIO.read('./data/lipocalin_family_aligned.fa', 'fasta')
    sequences = []
    for alseq in aligned_file:
        sequences.append(alseq.seq.__str__())
    hyperresidues = get_hyperresidues(sequences, aligned_hyperedges)
    for he in hyperresidues:
        #print(hyperresidues[he])
        for item in hyperresidues[he]:
            print(hyperresidues[he][item])
            print(hyperresidues[he][item] / 131.0)
    #print(sequences)
    #for alhe in aligned_file:
    
    
    
    
    #with open('./processed/1XKI_hyperedges.txt', 'w') as out:
    #    for he in get_hyperedges('./data/1XKI.pdb', chain_keys=['A']):
    #        out.write(he.__str__())
    #        out.write('\n')

if __name__ == '__main__':
    print('running main')
    main()
