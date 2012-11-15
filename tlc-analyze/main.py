#! /usr/bin/python

import hypergraph


def main():
    h = hypergraph.Hypergraph('./data/3EYC.pdb')
    h.hyper_analyze()
    #for e in h.edge_weights:
    #    print('{0} : {1}'.format(e, h.edge_weights[e]))
    #for l in h.edge_weights_by_seq:
    #    print(h.edge_weights_by_seq.index(l), l)

if __name__ == '__main__':
    print('running main')
    main()
