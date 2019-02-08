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
        pbg = pybel.from_path('data/large_corpus.bel')
        pybel.to_pickle(pbg, 'output/large_corpus_pybel.pkl')
    # Get all variant sites from the graph
    elif sys.argv[1] == 'get_pybel_mod_agents':
        pbg = pybel.from_pickle('output/large_corpus_pybel.pkl')
        mod_nodes = [get_agent(n) for n in pbg.nodes()
                     if has_protein_modification(n)]
        with open('output/bel_mod_agents.pkl', 'wb') as f:
            pickle.dump(mod_nodes, f)
    elif sys.argv[1] == 'get_pybel_stmts_by_site':
        pbg = pybel.from_pickle('output/large_corpus_pybel.pkl')
        pbp = bel.process_pybel_graph(pbg)
        sites = get_mod_sites(pbp.statements)
        with open('output/bel_mod_stmts_by_site.pkl', 'wb') as f:
            pickle.dump(sites, f)
    else:
        sys.exit(1)
