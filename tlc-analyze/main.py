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
    plt.plot(holo.unaligned_edge_weights_by_seq, 'b-',
             apo.unaligned_edge_weights_by_seq, 'r-')
    plt.ylabel('Residue Weight (sum of connected edges)')
    plt.xlabel('TLC Sequence Residues')
    plt.show()

if __name__ == '__main__':
    print('running main')
    main()
