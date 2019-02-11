import sys
import pickle
import pybel
from pybel.struct.filters import has_protein_modification
from indra.sources import bel
from indra.sources.bel.processor import get_agent
from .util import get_mod_sites

if __name__ == '__main__':
    # Parse the BEL script, takes a few minutes
    if sys.argv[1] == 'parse_belscript':
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        pbg = pybel.from_path(input_file)
        pybel.to_pickle(pbg, output_file)
    # Get all variant sites from the graph
    #elif sys.argv[1] == 'get_pybel_mod_agents':
    #    pbg = pybel.from_pickle('output/large_corpus_pybel.pkl')
    #    mod_nodes = [get_agent(n) for n in pbg.nodes()
    #                 if has_protein_modification(n)]
    #    with open('output/bel_mod_agents.pkl', 'wb') as f:
    #        pickle.dump(mod_nodes, f)
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
