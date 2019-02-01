import sys
import pickle
from protmapper import phosphosite_client as pc
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
from indra.databases import uniprot_client


def get_reader_sites(reader_sites, reader):
    sites = set()
    for ms in reader_sites[reader]:
        if ms.error_code is not None:
            continue
        if ms.valid is True:
            sites.add((ms.up_id, ms.orig_res, ms.orig_pos))
        elif ms.valid is False:
            if ms.mapped_res and ms.mapped_pos:
                sites.add((ms.up_id, ms.mapped_res, ms.mapped_pos))
    return sites


def filter_human(sites):
    return set([s for s in sites if uniprot_client.is_human(s[0])])


if __name__ == '__main__':
    psp_sites = pc.sites_only(exclude_isoforms=True)
    psp_sites = filter_human(psp_sites)
    # Load the reader sites
    with open('output/reader_sites.pkl', 'rb') as f:
        rs = pickle.load(f) # Contains MappedSite objects
        reach_sites = filter_human(get_reader_sites(rs, 'reach'))
        sparser_sites = filter_human(get_reader_sites(rs, 'sparser'))
        reader_sites = reach_sites | sparser_sites
    # Load the sites from PSP
    #with open(sys.argv[2], 'rb') as f:
    #    psp_sites = pickle.load(f) # List of INDRA Agents
    plt.ion()
    plt.figure()
    venn2((psp_sites, reader_sites),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser'))
    plt.savefig('plots/psp_reader_site_overlap.pdf')
    plt.figure()
    venn2((reach_sites, sparser_sites),
          set_labels=('REACH', 'Sparser'))
    plt.savefig('plots/reach_sparser_site_overlap.pdf')

