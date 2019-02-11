import sys
import pickle
import itertools
from collections import Counter
from indra.tools import assemble_corpus as ac


def get_db_phos_stmts(filename):
    from indra_db.client import get_statements_by_gene_role_type
    phos_stmts = get_statements_by_gene_role_type(
                        stmt_type='Phosphorylation', fix_refs=False,
                        preassembled=False,
                        with_evidence=True, with_support=False)
    with open(filename, 'wb') as f:
        pickle.dump(phos_stmts, f)
    return phos_stmts


def preprocess_db_stmts(stmts, output_file):
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
    # Filter to statements with residue and position
    site_stmts = [s for s in uniq_stmts if s.residue and s.position]
    # Organize into a dictionary indexed by site
    ac.dump_statements(site_stmts, output_file)
    return site_stmts


def get_stmts_by_site(phos_stmts, basename):
    # First filter statements to those that have objects with uniprot IDs
    readers = ['reach', 'sparser']
    for reader in readers:
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
        filename = '%s_%s.sites.pkl' % (basename, reader)
        with open(filename, 'wb') as f:
            pickle.dump(stmts_by_site, f)


def get_reader_sites(input_file):
    input_stmts = ac.load_statements(input_file)
    readers = ('reach', 'sparser')
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
    # Map grounding, remove identical statements
    elif sys.argv[1] == 'preprocess_stmts':
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        input_stmts = ac.load_statements(input_file)
        preproc_stmts = preprocess_db_stmts(input_stmts, output_file)
    elif sys.argv[1] == 'stmts_by_site':
        input_file = sys.argv[2]
        basename = sys.argv[3]
        input_stmts = ac.load_statements(input_file)
        get_stmts_by_site(input_stmts, basename)
    elif sys.argv[1] == 'reader_sites':
        input_file = sys.argv[2]
        get_reader_sites(input_file)
    else:
        print("Unrecognized arguments.")

