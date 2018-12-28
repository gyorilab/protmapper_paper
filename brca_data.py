import re
import csv
import sys
import pickle
import numpy as np
import pandas as pd
from indra.databases import uniprot_client


BRCA_DATA = 'data/breast_phosphosites.txt'
UP_MAPPINGS = 'data/HUMAN_9606_idmapping.dat'
BRCA_MAPPED = 'output/brca_up_mappings.txt'

if __name__ == '__main__':
    if sys.argv[1] == 'map_uniprot':
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
                    valid = []
                    for up_id in up_ids:
                        valid.append(uniprot_client.verify_location(up_id, res,
                                                                    pos))
                    up_id_str = '|'.join(up_ids)
                    valid_str = '|'.join([str(v) for v in valid])
                    site = (gene_name, rs_id, up_id_str, res, pos, valid_str)
                    if not np.all(valid):
                        print("%d: Site not valid: %s" % (ix, str(site)))
                brca_data.append(site)
        with open(BRCA_MAPPED, 'wt') as f:
            csvwriter = csv.writer(f, delimiter='\t')
            csvwriter.writerows()
    elif sys.argv[1] == 'site_stats':
        # Load INDRA statements, sorted by site
        indra_sites_file = sys.argv[2]
        with open(indra_sites_file, 'rb') as f:
            stmts_by_site = pickle.load(f)
        df = pd.read_csv(BRCA_MAPPED, delimiter='\t')
        total_sites = len(df)
        no_up_id = 0
        no_valid_up = 0
        for gene_name, rs_id, up_id, res, pos, valid in df.values:
            if up_id is np.nan:
                no_up_id += 1
                assert valid is np.nan
            elif valid is not np.nan:
                valids = [bool(v) for v in valid.split('|')]
                up_ids = up_id.split('|')
                if not np.any(valids):
                    no_valid_up += 1
        text = ("No UP ID: %d / %d (%.1f)\n" %
                (no_up_id, total_sites, (100*no_up_id / total_sites)))
        text += ("No valid UP sequence: %d / %d (%.1f)\n" %
                 (no_valid_up, total_sites, (100*no_valid_up / total_sites)))
        print(text)
        with open(sys.argv[2], 'wt') as f:
            f.write(text)
        # There are multiple uniprot IDs associated with each refseq ID.
        # To figure out whether it has any annotations, we turn each BRCA
        # site into a tuple (up_id, res, pos), and look this up in a dict
        # of sites with their annotations. If there are matches based on
        # ANY of the uniprot IDs we count that as a match.
    else:
        print("Argument must be one of map_uniprot, site_stats")
        sys.exit(1)
    """
    Stats
    - Number of proteins in Refseq with no UP ID mapping
    - Number of cases where site is not valid despite mapping (!)
    """
