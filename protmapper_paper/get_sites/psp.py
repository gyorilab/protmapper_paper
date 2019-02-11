import sys
import pickle
import pandas as pd
from indra.statements import *
from indra.databases import uniprot_client

# Read the data
data_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(data_file, delimiter='\t', skiprows=3)

agents = []
relations_by_site = {}

for kin_up, sub_up, respos in df[['KIN_ACC_ID', 'SUB_ACC_ID',
                                  'SUB_MOD_RSD']].values:
    residue = respos[0]
    position = respos[1:]
    mc = ModCondition('phosphorylation', residue, position)
    gene_name = uniprot_client.get_gene_name(sub_up)
    if gene_name is None:
        continue
    sub_ag = Agent(gene_name, mods=[mc], db_refs={'UP': sub_up})
    agents.append(sub_ag)

    site_key = (sub_up, residue, position)
    if site_key in relations_by_site:
        relations_by_site[site_key].append(kin_up)
    else:
        relations_by_site[site_key] = [kin_up]

with open(output_file, 'wb') as f:
    pickle.dump(agents, f)

#with open('output/psp_relations_by_site.pkl', 'wb') as f:
#    pickle.dump(relations_by_site, f)

