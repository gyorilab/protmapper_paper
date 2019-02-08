import sys
from indra.sources import signor
from indra.tools import assemble_corpus as ac
from indra.statements import Phosphorylation, Dephosphorylation, ModCondition


if __name__ == '__main__':
    pkl_file = sys.argv[1]
    sp = signor.process_from_web()
    
    for stmt in sp.statements:
    phos = ac.filter_by_type(sp.statements, Phosphorylation)
    dephos = ac.filter_by_type(sp.statements, Dephosphorylation)
    respos_stmts = [s for s in (phos + dephos) if s.residue and s.position]
    modified_agents = []
    for s in respos_stmts:
        sub_agent = s.sub
        if 'UP' in sub_agent.db_refs:
            
    ac.dump_statements(respos_stmts, pkl_file)
    
