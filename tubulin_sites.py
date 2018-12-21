"""
Read the Tubulin stabilizer dataset from Excel, parse into a Pandas dataframe,
and generate a list of all sites observed in the data with Uniprot IDs and
motif.
"""

from os.path import join
import pickle
from collections import defaultdict
import pandas as pd
from indra.util import write_unicode_csv
from indra.sources import indra_db_rest as idr
from indra.databases import hgnc_client, uniprot_client

site_list = []

filename_tc  = "Mature_Neuron_MT_pMS_Stabilizer1_2_time_course_ys.xlsx"
df_tc = pd.read_excel(join('ms_data', filename_tc))

p = lambda x: x.split("|")[1]

df_tc["UniprotId"] = df_tc["Protein Id"].apply(p)

import ipdb; ipdb.set_trace()

df_cols = df_tc[['site_id', 'UniprotId', 'Site Position', 'Motif']]

sites = []
for site_id, up_id, pos, motif in df_cols.values:
    if ';' in str(pos):
        pos_list = str(pos).split(';')
        motif_list = motif.split(';')
        assert len(pos_list) == len(motif_list)
        for ix in range(len(pos_list)):
            res = motif_list[ix][6]
            if res not in ('S', 'T', 'Y'):
                continue
            sites.append((site_id, up_id, res, pos_list[ix], motif_list[ix]))
    else:
        res = motif[6]
        if res not in ('S', 'T', 'Y'):
            continue
        sites.append((site_id, up_id, res, str(pos), motif))

write_unicode_csv('tubulin_sites.txt', sites, delimiter=',')

iso_sites = [s for s in sites if '-' in s[1]]
iso_only = []
pos_valid_in_ref = []
mappable_in_ref = []

for site in iso_sites[0:200]:
    site_id, up_id, res, pos, motif = site
    ref_iso = up_id.split('-')[0]
    ref_seq = uniprot_client.get_sequence(ref_iso)
    assert len(motif) == 13
    ref_motif_start = ref_seq.find(motif)
    if ref_motif_start == -1:
        iso_only.append(site)
    else:
        ref_pos = str(ref_motif_start + 6)
        if ref_pos == pos:
            pos_valid_in_ref.append(site)
        else:
            mappable_in_ref.append(site)



"""
with open('ms_data/tubulin_sorted_site_list.txt', 'rt') as f:
    for line in f.readlines():
        gene_name, respos = line.strip().split('_')
        res = respos[0]
        pos = respos[1:]
        site_list.append((gene_name, res, pos))

with open('../indra_apps/tubulin/work/phospho_stmts.pkl', 'rb') as f:
    print("Loading statements")
    stmts = pickle.load(f)

stmts_by_site = defaultdict(list)
for s in stmts:
    site = (s.sub.name, s.residue, s.position)
    stmts_by_site[site].append(s)

no_stmts = []
for gene, res, pos in site_list[0:50]:
    hgnc_id = hgnc_client.get_hgnc_id(gene)
    if hgnc_id is None:
        print("No HGNC ID for %s, skipping" % gene)
        continue
    if (gene, res, pos) not in stmts_by_site:
        print("Getting sequence for %s" % gene)
        up_id = hgnc_client.get_uniprot_id(hgnc_id)
        sequence = uniprot_client.get_sequence(up_id)
        seq_start = int(pos) - 8
        if seq_start < 0:
            seq_start = 0
        seq_end = int(pos) + 7
        if seq_end > len(sequence):
            seq_end = len(sequence)
        no_stmts.append((gene, res, pos, sequence[seq_start:seq_end]))

"""
