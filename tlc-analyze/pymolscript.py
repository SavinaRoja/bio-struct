import os

os.chdir('./workspace/bio-struct/tlc-analyze')

import hyperpymol


h = hyperpymol.get_struct_hypergraph('3eyc')
hyperpymol.heat_map(h.unaligned_edge_weights_by_seq, offset=13)
