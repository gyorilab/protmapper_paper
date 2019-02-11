import os
import sys
import pickle

FILE_SOURCE_MAP = {
    'PathwayCommons10.hprd.BIOPAX.sites.pkl': 'hprd',
    'PathwayCommons10.kegg.BIOPAX.sites.pkl': 'kegg',
    'PathwayCommons10.panther.BIOPAX.sites.pkl': 'panther',
    'PathwayCommons10.pid.BIOPAX.sites.pkl': 'pid',
    'PathwayCommons10.reactome.BIOPAX.sites.pkl': 'reactome',
    'PathwayCommons10.wp.BIOPAX.sites.pkl': 'wp',
     #'psp_kinase_substrate_tsv.sites.pkl': 'psp',
    'Kinase_substrates.sites.pkl': 'psp',
    'bel_large_corpus.sites.pkl': 'bel',
    'signor.sites.pkl': 'signor',
    'indra_reach.sites.pkl': 'reach',
    'indra_sparser.sites.pkl': 'sparser',
}


if __name__ == '__main__':
    # List of site dict files to combine
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    all_sites = {}
    for site_file in input_files:
        source = FILE_SOURCE_MAP[os.path.basename(site_file)]
        with open(site_file, 'rb') as f:
            site_dict = pickle.load(f)
            for site in site_dict:
                if site not in all_sites:
                    all_sites[site] = {'lhs': {}, 'rhs': {}}
                all_sites[site]['lhs'][source] = site_dict[site]['lhs']
                all_sites[site]['rhs'][source] = site_dict[site]['rhs']
    with open(output_file, 'wb') as f:
        pickle.dump(all_sites, f)
