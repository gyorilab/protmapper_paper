"""This script assembles a set of INDRA Statements that represent
exactly the set of annotations shown in the Venn diagram. An annotation
is defined as a specific substrate identity, residue, and position,
as well as a specific controller, both subject to some filtering conditions."""

import sys
import tqdm
import pickle
from collections import defaultdict
from indra.databases import uniprot_client
from indra.statements import Agent, Phosphorylation

# In these cases name normalization is not available, so we map them manually
name_mapping = {
    ('UP', 'P04637'): 'HBX',
    ('UP', 'P0CW71'): 'MCP1'
}


def group_by_controller(phos_stmts):
    stmts_by_controller = defaultdict(list)
    names = {}
    for stmt in phos_stmts:
        if not stmt.enz:
            continue
        if 'UP' in stmt.enz.db_refs:
            db_ns = 'UP'
            db_id = stmt.enz.db_refs['UP']
        else:
            db_ns, db_id = stmt.enz.get_grounding()
        if not db_ns or not db_id:
            continue
        # We collect names separately for the rare corner case where
        # the name is not normalized by grounding
        names[(db_ns, db_id)] = name_mapping[(db_ns, db_id)] \
            if (db_ns, db_id) in name_mapping else stmt.enz.name
        stmts_by_controller[(db_ns, db_id)].append(stmt)
    stmts_by_controller = {
        (names[(db_ns, db_id)], db_ns, db_id): stmts
        for (db_ns, db_id), stmts in stmts_by_controller.items()
    }
    return stmts_by_controller


def get_base_agent(agent):
    return Agent(agent.name, db_refs=agent.db_refs)


def squash_statements(stmt_group, key):
    ctrl_key, (up_id, res, pos) = key
    # Normalize the name here based on the key to handle some corner cases
    stmt_group[0].enz.name = ctrl_key[0]
    # We create a base enzyme agent to make sure any special state
    # is avoided
    enz_base = get_base_agent(stmt_group[0].enz)
    sub_base = get_base_agent(stmt_group[0].sub)
    sub_base.db_refs['UP'] = up_id
    # The substrate and the residue position are shared in this group
    # so we can use those directly
    stmt = Phosphorylation(enz_base, sub_base, res, pos)
    # We can now combine evidences for all Statements from all source
    evs = []
    for stmt_in_group in stmt_group:
        evs += stmt_in_group.evidence
    stmt.evidence = evs
    return stmt


def filter_statements(stmts):
    # Make sure the enzyme has protein-like grounding
    stmts = [s for s in stmts
             if s.enz.get_grounding()[0] in {'UP', 'FPLX', 'HGNC'}]
    # Make sure the substrate is a human protein
    stmts = [s for s in stmts if s.sub.db_refs.get('UP')
             and uniprot_client.is_human(s.sub.db_refs.get('UP'))]
    return stmts


def get_statements(sites_file, mappings_file, output_pkl=None):
    # Read in the sites collected and the result of mappings
    with open(sites_file, 'rb') as fh:
        sites = pickle.load(fh)
    with open(mappings_file, 'rb') as fh:
        mappings = pickle.load(fh)

    # We now collect all the statements based on the mapped site
    stmts_by_mapped_site = defaultdict(list)
    for site, data in tqdm.tqdm(sites.items()):
        mapping_result = mappings[site]
        # If there is no mapping result or the site is invalid and
        # couldn't be mapped then we skip these statements
        if not mapping_result or (not mapping_result.valid and
                                  not mapping_result.mapped_id):
            continue
        # If a mapping was done, we take the mapped values here, otherwise
        # we take the original value
        up_id = mapping_result.up_id if mapping_result.mapped_id is None \
            else mapping_result.mapped_id
        pos = mapping_result.orig_pos if mapping_result.valid \
            else mapping_result.mapped_pos
        res = mapping_result.orig_res if mapping_result.valid \
            else mapping_result.mapped_res

        # Skip corner case of non-human->human mapping to be consistent
        # with Venn diagram
        if not uniprot_client.is_human(mapping_result.up_id) and \
                mapping_result.mapped_id is not None and \
                uniprot_client.is_human(mapping_result.mapped_id):
            continue

        # We can now take all Statements from all sources and add them to
        # the list of Statements for the given mapped site
        for source, stmts in data['rhs'].items():
            # Add metadata on mapping to evidence
            for stmt in stmts:
                for ev in stmt.evidence:
                    ev.annotations['site_mapping'] = mapping_result.to_json()
                    ev.source_api = source
            stmts_by_mapped_site[(up_id, res, pos)] += stmts

    # We next group Statements for a given mapped site by the controller
    stmts_by_controller_and_mapped_site = defaultdict(list)
    for (up_id, res, pos), stmts in stmts_by_mapped_site.items():
        # Group statements by controller
        stmt_groups = group_by_controller(stmts)
        for ctrl_key, stmt_group in stmt_groups.items():
            stmts_by_controller_and_mapped_site[
                (ctrl_key, (up_id, res, pos))
            ] += stmt_group

    # We can finally squash Statements within each group, make sure the final
    # statements take on mapped residue and position values, and then
    # put them in a flat list.
    statements = []
    for key, stmt_group in tqdm.tqdm(stmts_by_controller_and_mapped_site.items()):
        stmt = squash_statements(stmt_group, key)
        # We now set the original/valid or mapped residue and position
        # at the Statement level
        stmt.residue = key[1][1]
        stmt.position = key[1][2]
        statements.append(stmt)

    # We now filter Statements to the subset of interest
    statements = filter_statements(statements)

    # Dump Statements into a pickle
    if output_pkl:
        with open(output_pkl, 'wb') as fh:
            pickle.dump(statements, fh)

    return statements


if __name__ == '__main__':
    sites_file, mappings_file, output_file = sys.argv[1:4]
    stmts = get_statements(sites_file, mappings_file, output_file)
