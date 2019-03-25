import sys
import pickle
import itertools
from collections import Counter
from indra.tools import assemble_corpus as ac
from indra_db.util import get_primary_db, get_raw_stmts_frm_db_list


def get_db_phos_stmts(filename):
    from indra_db.client import get_statements_by_gene_role_type
    phos_stmts = get_statements_by_gene_role_type(
                        stmt_type='Phosphorylation', fix_refs=False,
                        preassembled=False,
                        with_evidence=True, with_support=False)
    with open(filename, 'wb') as f:
        pickle.dump(phos_stmts, f)
    return phos_stmts


def get_db_agent_mod_stmts(filename):
    def has_mod_agents(stmt):
        mod_agents = []
        for agent in stmt.agent_list():
            if agent is not None:
                for mc in agent.mods:
                    if has_site_pos(mc):
                        return True
        return False

    def has_site_pos(mc):
        return mc.position is not None and mc.residue is not None
    batch_size = 100000
    db = get_primary_db()
    site_stmts = []
    for idx, db_stmt_batch in db.select_all_batched(
        batch_size, db.RawStatements, db.RawStatements.reading_id.isnot(None)):
        stmt_tuples = get_raw_stmts_frm_db_list(db, db_stmt_batch,
                                                fix_refs=False)
        stmts = [s[1] for s in stmt_tuples]
        for stmt in stmts:
            if has_mod_agents(stmt):
                site_stmts.append(stmt)
        print('Finished batch %d' % idx)
        print('Currently have %d site statements' % len(site_stmts))
        with open(filename, 'wb') as f:
            pickle.dump(site_stmts, f)
    return site_stmts


def preprocess_db_stmts(stmts, output_file, filter_stmt_site):
    """Take the statements from the database and grounding map them; """
    print("Mapping grounding")
    gmap_stmts = ac.map_grounding(stmts)
    #ac.dump_statements(gmap_stmts, prefix + '_gmap.pkl')
    print("Sorting and filtering")
    # Next, eliminate exact duplicates
    stmts_by_deep_hash = [(s.get_hash(shallow=False), s) for s in gmap_stmts]
    stmts_by_deep_hash.sort(key=lambda x: x[0])
    uniq_stmts = []
    for k, group in itertools.groupby(stmts_by_deep_hash, key=lambda x: x[0]):
        uniq_stmts.append(list(group)[0][1])
    if filter_stmt_site:
        # Filter to statements with residue and position
        site_stmts = [s for s in uniq_stmts if s.residue and s.position]
    else:
        site_stmts = uniq_stmts
    # Organize into a dictionary indexed by site
    ac.dump_statements(site_stmts, output_file)
    return site_stmts


def get_reader_agent_mod_stmts_by_site(agent_mod_stmts, reader, filename):
    # First filter statements to those that have objects with uniprot IDs
    stmts_by_site = {}
    # Filter to stmts for this reader
    reader_stmts = [s for s in agent_mod_stmts
                    if s.evidence[0].source_api == reader]
    for s in reader_stmts:
        for agent in s.agent_list():
            if agent is None:
                continue
            up_id = agent.db_refs.get('UP')
            if not up_id:
                continue
            for mc in agent.mods:
                if mc.residue is None or mc.position is None or \
                    mc.residue not in ('S', 'T', 'Y'):
                    continue
                site = (up_id, mc.residue, mc.position)
                if site not in stmts_by_site:
                    stmts_by_site[site] = {'lhs': [], 'rhs': []}
                stmts_by_site[site]['lhs'].append(s)
    with open(filename, 'wb') as f:
        pickle.dump(stmts_by_site, f)


def get_reader_stmts_by_site(phos_stmts, reader, filename):
    # First filter statements to those that have objects with uniprot IDs
    stmts_by_site = {}
    # Filter to stmts for this reader
    reader_stmts = [s for s in phos_stmts
                    if s.evidence[0].source_api == reader]
    for s in reader_stmts:
        up_id = s.sub.db_refs.get('UP')
        # Filter to stmts with substrate UP ID, residue and position
        if up_id is None or \
           s.residue is None or s.position is None or \
           s.residue not in ('S', 'T', 'Y'):
            continue
        site = (up_id, s.residue, s.position)
        if site not in stmts_by_site:
            stmts_by_site[site] = {'lhs': [], 'rhs': []}
        stmts_by_site[site]['rhs'].append(s)
    with open(filename, 'wb') as f:
        pickle.dump(stmts_by_site, f)


# TODO: generalize this to agent mod sites
def get_reader_sites(input_file):
    input_stmts = ac.load_statements(input_file)
    readers = ('reach', 'sparser', 'rlimsp')
    pm = ProtMapper(use_cache=True, cache_path=CACHE_PATH)
    sites_by_reader = {}
    # For all readers
    for reader in readers:
        sites = []
        # Filter to stmts for this reader
        reader_stmts = [s for s in input_stmts
                        if s.evidence[0].source_api == reader]
        for s in reader_stmts:
            up_id = s.sub.db_refs.get('UP')
            # Filter to stmts with substrate UP ID, residue and position
            if up_id is None or s.residue is None or s.position is None:
                continue
            if s.residue not in ('S', 'T', 'Y'):
                continue
            site = (up_id, s.residue, s.position)
            # Get the mapped site for the residue
            ms = pm.map_to_human_ref(up_id, 'uniprot', s.residue, s.position)
            sites.append(ms)
        # Group, tabulate frequency
        site_ctr = Counter(sites)
        # Store in dict
        sites_by_reader[reader] = site_ctr
    # Save sites
    with open('output/reader_sites.pkl', 'wb') as f:
        pickle.dump(sites_by_reader, f)
    # Save cache
    pm.save_cache()


if __name__ == '__main__':
    # Get statements from INDRA database
    if sys.argv[1] == 'get_phos_stmts':
        get_db_phos_stmts(sys.argv[2])
    elif sys.argv[1] == 'get_agent_mod_stmts':
        get_db_agent_mod_stmts(sys.argv[2])
    # Map grounding, remove identical statements
    elif sys.argv[1] == 'preprocess_stmts':
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        filter_stmt_site = True if sys.argv[4] == 'true' else False
        input_stmts = ac.load_statements(input_file)
        preproc_stmts = preprocess_db_stmts(input_stmts, output_file,
                                            filter_stmt_site)
    elif sys.argv[1] == 'stmts_by_site':
        input_file = sys.argv[2]
        reader = sys.argv[3]
        filename = sys.argv[4]
        input_stmts = ac.load_statements(input_file)
        get_reader_stmts_by_site(input_stmts, reader, filename)
    elif sys.argv[1] == 'agent_mod_stmts_by_site':
        input_file = sys.argv[2]
        reader = sys.argv[3]
        filename = sys.argv[4]
        input_stmts = ac.load_statements(input_file)
        get_reader_agent_mod_stmts_by_site(input_stmts, reader, filename)
    elif sys.argv[1] == 'reader_sites':
        input_file = sys.argv[2]
        get_reader_sites(input_file)
    else:
        print("Unrecognized arguments.")
