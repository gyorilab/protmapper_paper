import numpy
import pickle
import pandas

df = pandas.read_csv('../output/site_info.csv')
df = df[(df['SOURCE'] == 'reach') | (df['SOURCE'] == 'sparser') | 
        (df['SOURCE'] == 'rlimsp')]
df = df[df['VALID'] == False]
freq = numpy.array(df.FREQ) / numpy.sum(df.FREQ)
samples = [list(numpy.random.multinomial(1, freq)).index(1)
           for _ in range(100)]
df_sample = df.iloc[samples]

with open('../output/indra_reach.sites.pkl', 'rb') as fh:
    reach_sites = pickle.load(fh)
with open('../output/indra_sparser.sites.pkl', 'rb') as fh:
    sparser_sites = pickle.load(fh)
with open('../output/indra_rlimsp.sites.pkl', 'rb') as fh:
    rlimsp_sites = pickle.load(fh)

sentences = []
for _, row in df_sample.iterrows():
    std = reach_sites if row['SOURCE'] == 'reach' else sparser_sites
    key = (row['UP_ID'], row['ORIG_RES'], str(int(row['ORIG_POS'])))
    stmts = std[key]['rhs']
    sentences.append(stmts[0].evidence[0].text)

df_sample['EVIDENCE'] = sentences
df_sample.to_excel('reading_sites_sample.xlsx')
