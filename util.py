from collections import defaultdict
from indra.sources import biopax
from indra.statements import ActiveForm, Phosphorylation, Dephosphorylation


def get_mod_sites(indra_stmts):
    lhs_sites = defaultdict(list)
    rhs_sites = defaultdict(list)

    def add_agent_mods(stmt, agent):
        phos_mods = [mc for mc in agent.mods
                     if mc.mod_type == 'phosphorylation' and
                        mc.residue and mc.position]
        for pm in set(phos_mods):
            site = (agent.db_refs['UP'], pm.residue, pm.position)
            lhs_sites[site].append(stmt)

    for stmt in indra_stmts:
        # Skip complexes
        if len(stmt.agent_list()) > 2:
            continue
        # If it's an active form, put the site in the lhs_site category
        if isinstance(stmt, ActiveForm):
            if 'UP' in stmt.agent.db_refs:
                add_agent_mods(stmt, stmt.agent)
        # For other statement types:
        else:
            assert len(stmt.agent_list()) == 2
            # First, check the subject agent for mod conditions; if so, add
            # LHS annotation
            subj = stmt.agent_list()[0]
            if subj.mods:
                add_agent_mods(stmt, subj)
            # Second, check if the statement is a Phos/Dephos stmt
            if (isinstance(stmt, Phosphorylation) or \
               isinstance(stmt, Dephosphorylation)) and \
               'UP' in stmt.sub.db_refs:
                site = (stmt.sub.db_refs['UP'], stmt.residue, stmt.position)
                rhs_sites[site].append(stmt)
    return lhs_sites, rhs_sites

