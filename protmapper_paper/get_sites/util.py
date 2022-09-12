from indra.statements import ActiveForm, Phosphorylation, Dephosphorylation


def get_mod_sites(indra_stmts):
    sites = {}

    def add_site(site, side, stmt):
        if site not in sites:
            sites[site] = {'lhs': [], 'rhs': []}
        sites[site][side].append(stmt)

    def add_agent_mods(stmt, agent):
        phos_mods = [mc for mc in agent.mods
                     if mc.mod_type == 'phosphorylation' and
                        mc.residue and mc.position]
        for pm in set(phos_mods):
            site = (agent.db_refs['UP'], pm.residue, pm.position)
            add_site(site, 'lhs', stmt)

    for stmt in indra_stmts:
        # Skip large complexes
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
            if subj and subj.mods and 'UP' in subj.db_refs:
                add_agent_mods(stmt, subj)
            # Second, check if the statement is a Phos/Dephos stmt
            if (isinstance(stmt, Phosphorylation) or
                isinstance(stmt, Dephosphorylation)) and \
                    ('UP' in stmt.sub.db_refs or 'UPISO' in stmt.sub.db_refs):
                if 'UPISO' in stmt.sub.db_refs:
                    up_id = stmt.sub.db_refs['UPISO']
                else:
                    up_id = stmt.sub.db_refs['UP']
                site = (up_id, stmt.residue, stmt.position)
                add_site(site, 'rhs', stmt)
    return sites

