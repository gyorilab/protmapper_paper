import sys
import json
import pickle
from matplotlib import pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from indra import belief
from indra.belief.skl import CountsScorer
from indra.ontology.bio import bio_ontology
from indra.preassembler import Preassembler
from indra.statements import stmts_from_json_file


if __name__ == '__main__':
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
    protmap_stmts = pa.combine_related(return_toplevel=False)
    with_db_ev = []
    no_db_ev = []
    for s in protmap_stmts:
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

    beliefs_by_hash = {stmt.get_hash(): stmt.belief for stmt in protmap_stmts}
    with open(belief_output, 'w') as fh:
        json.dump(beliefs_by_hash, fh, indent=1)

    # Plot reading-only Statements
    rf_beliefs = [s.belief for s in no_db_ev]
    hist_res = plt.hist(rf_beliefs, alpha=0.5, bins=20)
    plt.xlim(0, 1.02)
    plt.ylabel('Count')
    plt.xlabel('Belief')
