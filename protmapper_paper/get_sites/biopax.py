import sys
import pickle
from indra.sources import biopax
from .util import get_mod_sites


def save_phosphorylation_stmts(owl_file, pkl_file):
    if owl_file.endswith('.gz'):
        bp = biopax.process_owl_gz(owl_file)
    else:
        bp = biopax.process_owl(owl_file)
    sites = get_mod_sites(bp.statements)
    with open(pkl_file, 'wb') as f:
        pickle.dump(sites, f)
    return sites


if __name__ == '__main__':
    owl_file = sys.argv[1]
    pkl_file = sys.argv[2]
    save_phosphorylation_stmts(owl_file, pkl_file)
