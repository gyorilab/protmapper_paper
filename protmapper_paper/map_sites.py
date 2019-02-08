import sys
import pickle
from protmapper import ProtMapper

def map_sites(sites_dict):
    """Tabulate valid, invalid, and mapped sites from a set of Agents."""
    site_map = {}
    pm = ProtMapper()
    for site_ix, site in enumerate(sites_dict.keys()):
        if site_ix % 1000 == 0:
            print('%d of %d' % (site_ix, len(sites_dict)))
        up_id, res, pos = site
        try:
            ms = pm.map_to_human_ref(up_id, 'uniprot', res, pos)
            site_map[site] = ms
        except Exception as e:
            print("Error: %s" % str(e))
            print("up_id: %s, res %s, pos %s" % (up_id, res, pos))
    # Now that we've collected a list of all the sites, tabulate frequencies
    return site_map

if __name__ == '__main__':
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    # Load the sites dict
    with open(input_file, 'rb') as f:
        sites_dict = pickle.load(f)
    # Create a map from each site to its MappedSite result
    site_map = map_sites(sites_dict)
    with open(output_file, 'wb') as f:
        pickle.dump(site_map, f)

