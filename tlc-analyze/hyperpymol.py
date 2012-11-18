'''
A small module to facilitate the interface of PyMOL and my hypergraph code. It
will generally be used from the PyMOL interactive prompt. 
'''

from pymol import cmd
import hypergraph
import os.path
from collections import namedtuple

r = namedtuple('Routine',
               'file, unaligned, aligned, offset, family_msa')

routines = {'3EYC': r('./data/3EYC.pdb', './data/3eyc_a_from_pdb.fa',
                      './data/3eyc_parent.fa', 13,
                      './data/lipocalin_family_aligned.fa'),
            '1XKI': r('./data/1XKI.pdb', './data/1xki_a_from_pdb.fa',
                      './data/1xki_parent.fa', 12,
                      './data/lipocalin_family_aligned.fa')}

def get_struct_hypergraph(struct_name, load_struct=True):
    '''
    This method will instantiate a hypergraph object and execute its analysis,
    returning the hypergraph. It will also load the structure into PyMOL unless
    explicitly avoided.
    '''
    #A little safety in case the user adds .pdb or uses lower case
    struct = os.path.splitext(struct_name)[0].upper()
    try:
        routine = routines[struct]
    except KeyError, e:
        print('Routine not defined for this structure')
        raise e
    
    h = hypergraph.Hypergraph(routine.file)
    h.hyper_analyze(*routine[1:])
    
    if load_struct:
        cmd.load(routine.file)
        cmd.select('A', 'chain a')
        cmd.hide('all')
        cmd.show('cartoon', 'A')
        cmd.center('A')
    return h

def heat_map(start_color='black', end_color'white'):
    


