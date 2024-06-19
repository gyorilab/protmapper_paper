import sys
import pickle
from indra.sources import bel
from .util import get_mod_sites


if __name__ == '__main__':
    # Parse the BEL script, takes a few minutes
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    pbp = bel.process_belscript(input_file)
    sites = get_mod_sites(pbp.statements)
    with open(output_file, 'wb') as f:
        pickle.dump(sites, f)
