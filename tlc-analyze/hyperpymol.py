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

def get_color_rgb(color):
    colors = cmd.get_color_indices()
    for c in colors:
        if c[0] == color:
            return cmd.get_color_tuple(colors.index(c))
    else:
        print('Could not find a reference for that color')
        return None

def heat_map(seq_scores, start_color='black', end_color='white', offset=0):
    '''
    Generates a heat map of the hyperconservation score by residues in PyMOL
    '''
    start_rgb = get_color_rgb(start_color)
    end_rgb = get_color_rgb(end_color)
    diffs = [end_rgb[i] - start_rgb[i] for i in xrange(3)]
    min_score = min(seq_scores)
    max_score = max(seq_scores)
    diff_score = max_score - min_score
    tick = 0
    for i in seq_scores:
        cmd.select('heat_map', 'chain a and resi {0}'.format(i + offset))
        norm_score = (seq_scores[i] - min_score) / diff_score
        temp_color = [norm_score * diffs[i] + start_rgb[i] for i in xrange(3)]
        cmd.set_color('heat_col{0}'.format(tick), temp_color)
        cmd.color('heat_col{0}'.format(tick), 'heat_map')

