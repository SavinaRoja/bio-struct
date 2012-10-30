#! /usr/bin/python

#resi_type = ((C),(F,Y,W),(H,R,K),(N,D,Q,E),(S,T,P,A,G),(M,I,L,V))
resi_type = {'C': 1, 'F': 2, 'Y': 2, 'W': 2, 'H': 3, 'R': 3, 'K': 3, 'N': 4,
              'D': 4, 'Q': 4, 'E': 4, 'S': 5, 'T': 5, 'P': 5, 'A': 5, 'G': 5,
              'M': 6, 'I': 6, 'L': 6, 'V': 6, '-': 0}

def get_hyperresidues(aligned_sequences, aligned_hyperedges):
    '''
    This function calculates and returns the hyperresidues for the provided
    aligned sequences with the aligned hyperedges.
    '''
    
    hyperresidues = {}
    
    for hyperedge in aligned_hyperedges:  # Iterate over all hyperedges
        #We wish to collect all the hyperresidues for a given hyperedge
        hyperresidue_set = {}  # This will hold that collection
        #Now we iterate over all the sequences in the alignment
        for sequence in aligned_sequences:
            #We need to get the residues at the vertices of the hyperedge
            vertices = []
            for i in hyperedge:  # This will compose the vertices list
                vertices.append(resi_type[sequence[i]])
            vertices = tuple(vertices)
            if 0 in vertices:  # Include none without all residues
                continue
            elif vertices in hyperresidue_set:  # If this has been seen before
                hyperresidue_set[vertices] += 1
            else:  # This has not been seen before
                hyperresidue_set[vertices] = 1
        hyperresidues[hyperedge] = hyperresidue_set
    return hyperresidues