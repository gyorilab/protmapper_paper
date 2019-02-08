import sys
import pickle
from indra.sources import signor
from .util import get_mod_sites

if __name__ == '__main__':
    output_file = sys.argv[1]
    sp = signor.process_from_web()
    sites = get_mod_sites(sp.statements)
    with open(output_file, 'wb') as f:
        pickle.dump(sites, f)

