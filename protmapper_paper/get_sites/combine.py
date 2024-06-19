import os
import sys
import pickle

FILE_SOURCE_MAP = {
    'PathwayCommons12.pid.BIOPAX.sites.pkl': 'pid',
    'PathwayCommons12.reactome.BIOPAX.sites.pkl': 'reactome',
    'Kinase_substrates.sites.pkl': 'psp',
    'bel_large_corpus.sites.pkl': 'bel',
    'signor.sites.pkl': 'signor',
    'hprd.sites.pkl': 'hprd',
    'indra_reach.sites.pkl': 'reach',
    'indra_sparser.sites.pkl': 'sparser',
    'indra_rlimsp.sites.pkl': 'rlimsp',
    'indra_reach_agent_mod.sites.pkl': 'reach',
    'indra_sparser_agent_mod.sites.pkl': 'sparser',
    'indra_rlimsp_agent_mod.sites.pkl': 'rlimsp',
}


if __name__ == '__main__':
    # List of site dict files to combine
    output_file = sys.argv[1]
    input_files = sys.argv[2:]
    all_sites = {}
    for site_file in input_files:
        source = FILE_SOURCE_MAP[os.path.basename(site_file)]
        print('Loading %s from %s' % (source, site_file))
        with open(site_file, 'rb') as f:
            site_dict = pickle.load(f)
            for site in site_dict:
                if site not in all_sites:
                    all_sites[site] = {'lhs': {}, 'rhs': {}}
                # Make sure we can handle separate agent mod statements
                # for a given reader so we check if there is already
                # something there for the given site and only overwrite
                # if there isn't.
                if not all_sites[site]['lhs'].get(source):
                    all_sites[site]['lhs'][source] = site_dict[site]['lhs']
                if not all_sites[site]['rhs'].get(source):
                    all_sites[site]['rhs'][source] = site_dict[site]['rhs']
    with open(output_file, 'wb') as f:
        print('Writing into %s' % output_file)
        pickle.dump(all_sites, f)
