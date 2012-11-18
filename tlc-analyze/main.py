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
    #To compare the holo and apo, structures, artificially start holo at 12
    adjust_holo = [0] + holo.unaligned_edge_weights_by_seq
    apo = apo.unaligned_edge_weights_by_seq
    holo_and_apo = zip(adjust_holo, apo)
    with open('apo_holo.csv', 'w') as out:
        out.write(',')
        i = 0
        while i < len(apo):
            out.write('{0},'.format(i + 12))
            i +=1
        out.write('\n3EYC,')
        for j in adjust_holo:
            out.write('{0},'.format(j))
        out.write('\n1XKI,')
        for k in apo:
            out.write('{0},'.format(k))
    
    plt.plot(holo.unaligned_edge_weights_by_seq, 'b-',
             apo.unaligned_edge_weights_by_seq, 'r-',
             holo_and_apo, 'g-')
    plt.ylabel('Residue Weight (sum of connected edges)')
    plt.xlabel('TLC Sequence Residues')
    plt.show()

if __name__ == '__main__':
    print('running main')
    main()
