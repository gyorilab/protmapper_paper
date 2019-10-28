import tqdm
import pickle
from fuzzywuzzy import process
from collections import defaultdict

with open('output/indra_db_stmts_new.pkl', 'rb') as fh:
    new_stmts = pickle.load(fh)
    new_stmtsd = defaultdict(list)
    new_stmtst = {}
    for stmt in new_stmts:
        new_stmtsd[stmt.get_hash(shallow=True, refresh=True)].append(stmt)
        for ev in stmt.evidence:
            if ev.text and ev.pmid:
                new_stmtst[ev.text] = ev.pmid

with open('output/all_sites.pkl', 'rb') as fh:
    all_sites = pickle.load(fh)

unmatched_ev = []
unmatched = []
for _, std in all_sites.items():
    rhs = std.get('rhs', [])
    for source, old_stmts in rhs.items():
        if source not in {'reach', 'sparser'}:
            continue
        for old_stmt in old_stmts:
            matched = False
            old_hash = old_stmt.get_hash(shallow=True, refresh=True)
            if old_hash in new_stmtsd:
                for stmt in new_stmtsd[old_hash]:
                    if stmt.evidence[0].source_api == old_stmt.evidence[0].source_api:
                        if stmt.evidence[0].text == old_stmt.evidence[0].text:
                            old_stmt.evidence[0].pmid = stmt.evidence[0].pmid
                            matched = True
                            print('Match found for %s' % old_stmt)
                            break
                if not matched:
                    print('No evidence match found for %s' % old_stmt)
                    unmatched.append(old_stmt)
            else:
                print('No match found for %s' % old_stmt)
                unmatched.append(old_stmt)

really_unmatched = []
all_texts = list(new_stmtst.keys())
for unmatched_stmt in tqdm.tqdm(unmatched):
    pmid = new_stmtst.get(unmatched_stmt.evidence[0].text)
    if pmid:
        unmatched_stmt.evidence[0].pmid = pmid
        print('Text match found for %s' % unmatched_stmt)
    else:
        fuzz_txts = []
        for txt in all_texts:
            if txt.startswith(unmatched_stmt.evidence[0].text[:10]):
                fuzz_txts.append(txt)
        if not fuzz_txts:
            really_unmatched.append(unmatched_stmt)
            print('No text match found for %s' % unmatched_stmt)
        else:
            scores = process.extract(unmatched_stmt.evidence[0].text, fuzz_txts, limit=5)
            top_text = scores[0][0]
            pmid = new_stmtst[top_text]
            print('Scored matched: %s' % (','.join([str(s[1]) for s in scores])))
            if len(scores) >= 2 and scores[0][1]-scores[1][1] < 1:
                print('Score difference too low, no match foundd for %s' % unmatched_stmt)
                really_unmatched.append(unmatched_stmt)
            else:
                unmatched_stmt.evidence[0].pmid = pmid

with open('output/all_sites_updated.pkl', 'wb') as fh:
    pickle.dump(all_sites, fh)

