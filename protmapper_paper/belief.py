import sys
import json
import pandas
import pickle
import matplotlib
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from indra import belief
from indra.belief.skl import CountsScorer
from indra.ontology.bio import bio_ontology
from indra.preassembler import Preassembler
from indra.statements import stmts_from_json_file, stmts_to_json_file


def get_curation_data(stmts_file, curation_data_file,
                      training_corpus_file, training_corpus_curations_file,
                      training_corpus_extra_evidence_file):
    """This function can be used to produce
       - protmapper_belief_training_corpus.json
       - protmapper_belief_training_corpus_curations.json
       - protmapper_belief_training_corpus_extra_evidences.pkl
       from inputs:
       - bioexp_asmb_preassembled.pkl
       - curation_dataset_with_bg_psp.pkl
    The ebove files are made available in version control directly so aren't
    produced as part of the Makefile.
    """
    # Load pickle of bioexp statements.
    with open(stmts_file, 'rb') as fh:
        stmts = pickle.load(fh)

    # Get dataset of curated statements along with correctness values
    def stmts_for_df(df, stmts_by_hash):
        stmt_list = []
        for row in df.itertuples():
            stmt_hash = row.stmt_hash
            if stmt_hash not in stmts_by_hash:
                continue
            stmt_list.append(stmts_by_hash[stmt_hash])
        return stmt_list

    def load_curation_data(filename):
        with open(filename, 'rb') as f:
            dataset = pickle.load(f)
            df = pandas.DataFrame.from_records(dataset)
            df = df.fillna(0)
        # Every column except agent names and stmt type should be int
        dtype_dict = {col: 'int64' for col in df.columns
                      if col not in ('agA_name', 'agA_ns', 'agA_id', 'stmt_type',
                                     'agB_name', 'agB_ns', 'agB_id')}
        df = df.astype(dtype_dict)
        return df

    stmts_by_hash = {s.get_hash(): s for s in stmts}
    # Load the curated data
    kge_df = load_curation_data(curation_data_file)
    # Get statements from curation dataframe
    kge_stmts = stmts_for_df(kge_df, stmts_by_hash)
    kge_hashes = [s.get_hash() for s in kge_stmts]

    # Curated correctness values
    y_arr = kge_df['correct'].values
    y_arr = [int(val) for val in y_arr]  # Convert to int for JSON serialization

    # Get the more specific evidences for the curated statements
    extra_evidence = belief.get_ev_for_stmts_from_supports(kge_stmts)

    # Save the curated stmts
    stmts_to_json_file(kge_stmts, training_corpus_file)

    # Save the curations
    with open(training_corpus_curations_file, 'w') as f:
        json.dump(dict(zip(kge_hashes, y_arr)), f)

    # Save the extra evidences
    with open(training_corpus_extra_evidence_file, 'wb') as f:
        pickle.dump(extra_evidence, f)


if __name__ == '__main__':
    # This is to make sure preassembly happens with respect to the
    # groundings as used to find unique annotations
    from indra.statements import agent
    agent.default_ns_order = ['UP', 'HGNC', 'FPLX']

    annotation_stmts, training_corpus, curations, \
        curation_extra_evs, belief_output = sys.argv[1:6]
    with open(annotation_stmts, 'rb') as fh:
        protmap_stmts = pickle.load(fh)
    training_stmts = stmts_from_json_file(training_corpus)
    with open(curations, 'r') as fh:
        curations_by_hash = json.load(fh)
        # Convert hashes from strings to int after JSON deserialization
        curations_by_hash = {int(k): v for k, v in curations_by_hash.items()}
    with open(curation_extra_evs, 'rb') as fh:
        extra_evidence = pickle.load(fh)

    # Match the correctness values to the statements by hash
    correctness_arr = [curations_by_hash[s.get_hash()] for s in training_stmts]

    all_sources = list({ev.source_api
                        for stmt in training_stmts for ev in stmt.evidence})

    # Use a Random Forest classifier
    clf = RandomForestClassifier(n_estimators=2000, max_depth=13)

    # Provide all available features as well as the more specific evidences.
    # However, limit the source evidence to only what is obtained from readers
    scorer = CountsScorer(
        clf, all_sources, use_stmt_type=True, use_num_pmids=True,
        use_promoter=True, use_avg_evidence_len=True,
        include_more_specific=True,
        use_residue_position=True)

    # Train the model on the curated data
    scorer.fit(training_stmts, correctness_arr, extra_evidence)
    dbs = {'psp', 'signor', 'hprd', 'pid', 'bel', 'reactome'}

    # Assign every db source_api to biopax as a proxy
    for stmt in protmap_stmts:
        for ev in stmt.evidence:
            if ev.source_api in dbs:
                ev.source_api = 'biopax'

    pa = Preassembler(bio_ontology, protmap_stmts)
    assembled_stmts = pa.combine_related(return_toplevel=False)
    with_db_ev = []
    no_db_ev = []
    for s in assembled_stmts:
        sources = set([e.source_api for e in s.evidence])
        for supp_stmt in s.supports:
            for supp_ev in supp_stmt.evidence:
                sources.add(supp_ev.source_api)
        if 'biopax' in sources:
            with_db_ev.append(s)
        else:
            no_db_ev.append(s)
    for stmt in with_db_ev:
        stmt.belief = 1

    belief_engine = belief.BeliefEngine(scorer)
    belief_engine.set_hierarchy_probs(no_db_ev)

    beliefs_by_hash = {stmt.get_hash(): stmt.belief for stmt in assembled_stmts}
    with open(belief_output, 'w') as fh:
        json.dump(beliefs_by_hash, fh, indent=1)

    # Now set beliefs on the original statements
    for stmt in protmap_stmts:
        stmt.belief = beliefs_by_hash[stmt.get_hash()]

    with open('output/annotation_statements_with_belief.pkl', 'wb') as fh:
        pickle.dump(protmap_stmts, fh)

    # We now create two plots for all annotations and just from text mining
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42
    matplotlib.rcParams['font.sans-serif'] = "Arial"
    matplotlib.rcParams['font.family'] = "sans-serif"

    # Plot all Statements
    fig = plt.figure()
    beliefs = [s.belief for s in protmap_stmts]
    hist_res = plt.hist(beliefs, alpha=0.5, bins=20)
    plt.xlim(0, 1.02)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel('Number of annotations', fontsize=16)
    plt.xlabel('INDRA belief score for annotation', fontsize=16)
    fig.subplots_adjust(left=0.17, bottom=0.15, right=0.97, top=0.92)
    plt.savefig('plots/all_belief_hist.pdf')

    # Plot ones just with text mining support
    fig = plt.figure()
    rf_beliefs = [s.belief for s in no_db_ev]
    hist_res = plt.hist(rf_beliefs, alpha=0.5, bins=20)
    plt.xlim(0, 1.02)
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.ylabel('Number of annotations', fontsize=16)
    plt.xlabel('INDRA belief score for annotation', fontsize=16)
    fig.subplots_adjust(left=0.17, bottom=0.15, right=0.97, top=0.92)
    plt.savefig('plots/reading_only_belief_hist.pdf')

    print('Len: %d' % len(rf_beliefs))
    print('Min: %f' % min(rf_beliefs))
    print('Max: %f' % max(rf_beliefs))
    import numpy
    print('Median: %f' % numpy.median(rf_beliefs))
