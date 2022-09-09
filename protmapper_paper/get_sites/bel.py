import sys
import pickle
import pybel
from indra.sources import bel
from .util import get_mod_sites

if __name__ == '__main__':
    # Parse the BEL script, takes a few minutes
    if sys.argv[1] == 'parse_belscript':
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        pbg = pybel.from_bel_script(input_file)
        pybel.to_pickle(pbg, output_file)
    elif sys.argv[1] == 'get_pybel_stmts_by_site':
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        pbg = pybel.from_pickle(input_file)
        pbp = bel.process_pybel_graph(pbg)
        sites = get_mod_sites(pbp.statements)
        with open(output_file, 'wb') as f:
            pickle.dump(sites, f)
    else:
        sys.exit(1)
