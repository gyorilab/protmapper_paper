import pickle
import pandas as pd
from indra.statements import *
from indra.databases import uniprot_client

# Read the data
df = pd.read_csv('data/Kinase_Substrate_Dataset', delimiter='\t', skiprows=3)

agents = []

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

with open('output/psp_kinase_substrate_tsv.pkl', 'wb') as f:
    pickle.dump(agents, f)
