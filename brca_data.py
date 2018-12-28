import re
import csv
import numpy as np
from indra.databases import uniprot_client

BRCA_DATA = 'data/breast_phosphosites.txt'
UP_MAPPINGS = 'data/HUMAN_9606_idmapping.dat'
BRCA_MAPPED = 'output/brca_up_mappings.txt'

# Build Refseq -> Uniprot Mappings
rs_up = {}
with open(UP_MAPPINGS, 'rt') as f:
    csvreader = csv.reader(f, delimiter='\t')
    for up_id, other_db, other_id in csvreader:
        if other_db == 'RefSeq':
            if other_id in rs_up:
                rs_up[other_id].append(up_id)
            else:
                rs_up[other_id] = [up_id]
# Parse the BRCA data and add UP IDs where available
brca_data = []
with open(BRCA_DATA, 'rt') as f:
    lines = f.readlines()
    for ix, line in enumerate(lines[1:]): # Skip the header line
        if ix % 100 == 0:
            print(ix)
        m = re.match(r'(.*?)\.(.*):([sty])(\d+)', line)
        assert len(m.groups()) == 4
        gene_name, rs_id, res, pos = m.groups()
        res = res.upper()
        # Get Uniprot ID
        up_ids = rs_up.get(rs_id)
        if up_ids is None:
            print(f"No Uniprot ID for RefSeq %s" % rs_id)
            valid = None
            site = (gene_name, rs_id, up_ids, res, pos, valid)
        else:
            # Shouldn't matter which Uniprot ID we pick since the sequences
            # will be the same
            valid = []
            for up_id in up_ids:
                valid.append(uniprot_client.verify_location(up_id, res, pos))
            up_id_str = '|'.join(up_ids)
            valid_str = '|'.join([str(v) for v in valid])
            site = (gene_name, rs_id, up_id_str, res, pos, valid_str)
            if not np.all(valid):
                print("%d: Site not valid: %s" % (ix, str(site)))
        brca_data.append(site)
with open(BRCA_MAPPED, 'wt') as f:
    csvwriter = csv.writer(f, delimiter='\t')
    csvwriter.writerows()

"""
Stats
- Number of proteins in Refseq with no UP ID mapping
- Number of cases where site is not valid despite mapping (!)
"""
