import csv
import codecs
import sys
import gzip
import json
import pickle
import itertools
from collections import defaultdict, Counter
import pandas
import pystow
import tqdm
from indra.statements import stmts_from_json, Phosphorylation
from indra.tools import assemble_corpus as ac
from indra.util import batch_iter
from indra_db.util import get_primary_db, get_raw_stmts_frm_db_list


def get_all_indra_phos_stmts(filename, db_pkl, rlimsp_pkl):
    with open(db_pkl, 'rb') as fh:
        db_stmts = pickle.load(fh)
    with open(filename, 'wb') as fh:
        pickle.dump(db_stmts, fh)


class StatementJSONDecodeError(Exception):
    pass


def load_statement_json(json_str: str, attempt: int = 1, max_attempts: int = 5):
    try:
        return json.loads(json_str)
    except json.JSONDecodeError:
        if attempt < max_attempts:
            json_str = codecs.escape_decode(json_str)[0].decode()
            return load_statement_json(
                json_str, attempt=attempt + 1, max_attempts=max_attempts
            )
    raise StatementJSONDecodeError(
        f"Could not decode statement JSON after "
        f"{attempt} attempts: {json_str}"
    )


def get_db_phos_stmts(phos_fname, agent_fname):
    base_folder = pystow.module('indra', 'db')
    drop_readings_fname = base_folder.join(name='drop_readings.pkl')
    raw_stmts_fname = base_folder.join(name='raw_statements.tsv.gz')
    text_refs_fname = base_folder.join(name='text_refs_principal.tsv.gz')
    reading_text_content_fname = base_folder.join(
        name='reading_text_content_meta.tsv.gz')

    print("Loading text content metadata")
    df = pandas.read_csv(reading_text_content_fname, header=None,
                         sep='\t')
    df.sort_values(2, inplace=True)

    def load_text_refs_by_trid(fname):
        text_refs = {}
        for line in tqdm.tqdm(
                gzip.open(fname, "rt", encoding="utf-8"),
                desc="Processing text refs into a lookup dictionary",
        ):
            ids = line.strip().split("\t")
            id_names = ["TRID", "PMID", "PMCID", "DOI", "PII", "URL",
                        "MANUSCRIPT_ID"]
            d = {}
            for id_name, id_val in zip(id_names, ids):
                if id_val != "\\N":
                    if id_name == 'TRID':
                        id_val = int(id_val)
                    d[id_name] = id_val
            text_refs[int(ids[0])] = d
        return text_refs

    print("Loading dropped readings info")
    with open(drop_readings_fname, 'rb') as fh:
        drop_readings = pickle.load(fh)

    text_refs = load_text_refs_by_trid(text_refs_fname)
    reading_id_to_text_ref_id = dict(zip(df[0], df[3]))

    phos_stmts = []
    agent_site_stmts = []

    def agent_has_psite(agent):
        for mod in agent.mods:
            if mod.mod_type == 'phosphorylation' and mod.residue and mod.position:
                return True
        return False

    with gzip.open(raw_stmts_fname, 'rt') as fh:
        reader = csv.reader(fh, delimiter='\t')
        for lines in tqdm.tqdm(batch_iter(reader, 10000), total=7000):
            stmts_jsons = []
            for raw_stmt_id, db_info_id, reading_id, stmt_json_raw in lines:
                refs = None
                # Skip non-phosphorylations - quick pre-filter at the JSON
                # string level to find Phosphorylation statements and also
                # Agents with a phosphorylation mod condition
                if 'phosphorylation' not in stmt_json_raw.lower():
                    continue
                # We are focusing on reading here so we need to recover
                # text references
                if reading_id != '\\N':
                    # Skip if this is for a dropped reading
                    if int(reading_id) in drop_readings:
                        continue
                    text_ref_id = reading_id_to_text_ref_id.get(int(reading_id))
                    if text_ref_id:
                        refs = text_refs.get(text_ref_id)
                # We can now load the JSON and fix the text references
                stmt_json = load_statement_json(stmt_json_raw)
                if refs:
                    stmt_json['evidence'][0]['text_refs'] = refs
                    if refs.get('PMID'):
                        stmt_json['evidence'][0]['pmid'] = refs['PMID']
                stmts_jsons.append(stmt_json)
            # Deserialize the set of statements and filter
            stmts = stmts_from_json(stmts_jsons)
            stmts = ac.filter_evidence_source(stmts, ['reach', 'sparser', 'rlimsp'])
            stmts = ac.fix_invalidities(stmts, in_place=True)
            # Keep Phosphorylation statements with residue and position
            phos_stmts += [s for s in stmts if isinstance(s, Phosphorylation)
                           and s.residue and s.position]
            # Separately extract agent site statements
            agent_site_stmts += [s for s in stmts
                                 if any(agent_has_psite(a)
                                        for a in s.real_agent_list())]

    with open(phos_fname, 'wb') as f:
        pickle.dump(phos_stmts, f)
    with open(agent_fname, 'wb') as f:
        pickle.dump(agent_site_stmts, f)


def preprocess_db_stmts(stmts, output_file):
    """Take the statements from the database and grounding map them; """
    print("Mapping grounding")
    gmap_stmts = ac.map_grounding(stmts)
    #ac.dump_statements(gmap_stmts, prefix + '_gmap.pkl')
    print("Sorting and filtering")
    uniq_stmts = list({s.get_hash(shallow=False): s
                       for s in gmap_stmts}.values())
    # Organize into a dictionary indexed by site
    ac.dump_statements(uniq_stmts, output_file)
    return uniq_stmts


def get_reader_agent_mod_stmts_by_site(agent_mod_stmts, reader, filename):
    # First filter statements to those that have objects with uniprot IDs
    stmts_by_site = {}
    # Filter to stmts for this reader
    reader_stmts = [s for s in agent_mod_stmts
                    if s.evidence[0].source_api == reader]
    for s in reader_stmts:
        for agent in s.real_agent_list():
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


if __name__ == '__main__':
    # Get statements from INDRA database
    if sys.argv[1] == 'get_db_phos_stmts':
        get_db_phos_stmts(sys.argv[2], sys.argv[3])
    # Map grounding, remove identical statements
    elif sys.argv[1] == 'preprocess_stmts':
        input_file = sys.argv[2]
        output_file = sys.argv[3]
        input_stmts = ac.load_statements(input_file)
        preproc_stmts = preprocess_db_stmts(input_stmts, output_file)
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
    else:
        print("Unrecognized arguments.")
