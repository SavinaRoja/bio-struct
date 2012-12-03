'''
    Generates a heat map of the hyperconservation score by residues in PyMOL
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

#Holo structure
holo = get_struct_hypergraph('3eyc', load_struct=True)
#Apo structure
apo = get_struct_hypergraph('1xki', load_struct=False)
#Rudimentary summation of the two sequences
h_seq = holo.unaligned_edge_weights_by_seq
a_seq = apo.unaligned_edge_weights_by_seq
seq_scores = [i + j for i, j in zip(h_seq, a_seq)]

#Modify these to use different colors, use valid PyMOL named colors
start_color= 'yelloworange'
end_color = 'density'

#seq_offset = 12  # The first residue number of the sequence

def get_color_rgb(color):
    try:
        return cmd.get_color_tuple(color)
    except:
        print('Could not find a reference for that color')
        return None

start_rgb = get_color_rgb(start_color)
end_rgb = get_color_rgb(end_color)
diffs = [end_rgb[i] - start_rgb[i] for i in xrange(3)]
min_score = min(seq_scores)
diff_score = max(seq_scores) - min_score
for i in xrange(len(seq_scores)):
    cmd.select('heatmap', 'chain a and resi {0}'.format(i))
    norm_score = (seq_scores[i] - min_score) / diff_score
    temp_color = [start_rgb[j] + norm_score * diffs[j] for j in xrange(3)]
    color_name = 'heat_col{0}'.format(i)
    cmd.set_color(color_name, temp_color)
    cmd.color(color_name, 'heatmap')
