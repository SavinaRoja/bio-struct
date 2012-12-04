#! /usr/bin/python

import hypergraph
import matplotlib.pyplot as plt


def main():
    #Analyze 3EYC, TLC with 1,4-butanediol
    holo = hypergraph.Hypergraph('./data/3EYC.pdb')
    holo.hyper_analyze(unaligned='./data/3eyc_a_from_pdb.fa',
                       aligned='./data/3eyc_parent.fa',
                       offset=13,
                       family_msa='./data/lipocalin_family_aligned.fa')
    #for e in h.edge_weights:
    #    print('{0} : {1}'.format(e, h.edge_weights[e]))
    #for l in h.edge_weights_by_seq:  # Aligned sequence
    #    print(h.edge_weights_by_seq.index(l), l)
    #for l in h.unaligned_edge_weights_by_seq:  # Unaligned sequence
    #    print(l)
    apo = hypergraph.Hypergraph('./data/1XKI.pdb')
    apo.hyper_analyze(unaligned='./data/1xki_a_from_pdb.fa',
                      aligned='./data/1xki_parent.fa',
                      offset=12,
                      family_msa='./data/lipocalin_family_aligned.fa')
    with open('./processed/3EYC_aligned_hyperedges.txt', 'w') as out:
        for i in holo.aligned_hyperedges:
            out.write('{0}\n'.format(i))
    with open('./processed/3EYC_hyperresidues.txt', 'w') as out2:
        for i in holo.hyperresidues:
            for j in i:
                out2.write('{0}, '.format(j))
            out2.write('/n')
    holo_seq = holo.unaligned_edge_weights_by_seq
    apo_seq = apo.unaligned_edge_weights_by_seq
    with open('apo_holo.csv', 'w') as out:
        out.write(',')
        i = 0
        while i < len(apo_seq):
            out.write('{0},'.format(i))
            i +=1
        out.write('\n3EYC,')
        for j in holo_seq:
            out.write('{0},'.format(j))
        out.write('\n1XKI,')
        for k in apo_seq:
            out.write('{0},'.format(k))
    
    plt.plot(holo_seq, 'b-',
             apo_seq, 'r-')
    plt.ylabel('Residue Weight (sum of connected edges)')
    plt.xlabel('TLC Sequence Residues')
    plt.show()

if __name__ == '__main__':
    print('running main')
    main()
