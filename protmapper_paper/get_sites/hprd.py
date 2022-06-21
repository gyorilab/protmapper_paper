import sys
import pickle
from indra.sources import hprd
from indra.statements import Phosphorylation
from .util import get_mod_sites


if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    hp = hprd.process_archive(input_file)
    # Keep fully defined phosphorylations only
    stmts = [s for s in hp.statements if isinstance(s, Phosphorylation)
             and s.enz is not None]
    sites = get_mod_sites(stmts)
    with open(output_file, 'wb') as f:
        pickle.dump(sites, f)

