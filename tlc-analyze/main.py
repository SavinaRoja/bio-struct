#! /usr/bin/python

from hyperedges import get_hyperedges


def main():
    with open('./processed/3EYC_hyperedges.txt', 'w') as out:
        for he in get_hyperedges('./data/3EYC.pdb', chain_keys=['A']):
            out.write(he.__str__())
            out.write('\n')
    with open('./processed/1XKI_hyperedges.txt', 'w') as out:
        for he in get_hyperedges('./data/1XKI.pdb', chain_keys=['A']):
            out.write(he.__str__())
            out.write('\n')

if __name__ == '__main__':
    print('running main')
    main()
