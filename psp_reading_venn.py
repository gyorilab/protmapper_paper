import sys
import pickle
from collections import defaultdict
from protmapper import phosphosite_client as pc
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
from indra.databases import uniprot_client
from indra.util import read_unicode_csv

def get_kinases(stmts):
    kinases = set()
    for s in stmts:
        if s.enz is None:
            continue
        if 'UP' in s.enz.db_refs:
            up_id = s.enz.db_refs['UP']
            if uniprot_client.is_human(up_id):
                kinases.add(s.enz.db_refs['UP'])
        elif 'FPLX' in s.enz.db_refs:
            kinases.add(s.enz.db_refs['FPLX'])
    return kinases


def filter_human(sites):
    return {k: v for k, v in sites.items() if uniprot_client.is_human(k[0])}


def filter_human_list(sites):
    return [s for s in sites if uniprot_client.is_human(s[0])]


def get_reader_sites(reader_sites, readers, stmts_by_site):
    sites = defaultdict(set)
    for reader in readers:
        for ms in reader_sites[reader]:
            if ms.error_code is not None:
                continue
            # Get the statements for the site using origin res/pos
            try:
                kinases_for_site = get_kinases(stmts_by_site[
                                         (ms.up_id, ms.orig_res, ms.orig_pos)])
            except KeyError:
                kianses_for_site = set()
            if ms.valid is True:
                sites[(ms.up_id, ms.orig_res, ms.orig_pos)] |= kinases_for_site
            elif ms.valid is False:
                if ms.mapped_res and ms.mapped_pos:
                    sites[(ms.up_id, ms.mapped_res, ms.mapped_pos)] |= \
                                                              kinases_for_site
    return filter_human(sites)


def get_site_kinase_tuples(site_dict):
    return [site + (kin,) for site, kin_list in site_dict.items()
                          for kin in kin_list]


if __name__ == '__main__':
    # Load raw sites from PSP
    psp_sites = pc.sites_only(exclude_isoforms=True)
    psp_sites = filter_human_list(psp_sites)
    # Load the reader sites
    with open('output/indra_stmts_by_site.pkl', 'rb') as f:
        indra_stmts_by_site = pickle.load(f)
    with open('output/reader_sites.pkl', 'rb') as f:
        rs = pickle.load(f) # Contains MappedSite objects
        reach_sites = get_reader_sites(rs, ['reach'], indra_stmts_by_site)
        sparser_sites = get_reader_sites(rs, ['sparser'], indra_stmts_by_site)
        reader_sites = get_reader_sites(rs, ['reach', 'sparser'],
                                        indra_stmts_by_site)

    plt.ion()
    # Sites: PSP vs. Readers
    plt.figure()
    venn2((set(psp_sites), set(reader_sites.keys())),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser'))
    plt.savefig('plots/psp_reader_site_overlap.pdf')

    # Sites: REACH vs. Sparser
    plt.figure()
    venn2((set(reach_sites.keys()), set(sparser_sites.keys())),
          set_labels=('REACH', 'Sparser'))
    plt.savefig('plots/reach_sparser_site_overlap.pdf')

    # Load the annotated sites from PSP
    with open('output/psp_relations_by_site.pkl', 'rb') as f:
        ks_sites = pickle.load(f) # Dict of sites to lists of enzymes
    ks_human = filter_human(ks_sites)

    # Compile site-regulator tuples
    ks_annots = set(get_site_kinase_tuples(ks_human))
    reader_annots = set(get_site_kinase_tuples(reader_sites))
    reach_annots = set(get_site_kinase_tuples(reach_sites))
    sparser_annots = set(get_site_kinase_tuples(sparser_sites))

    # Annotations: PSP vs. Readers
    plt.figure()
    venn2((ks_annots, reader_annots),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser'))
    plt.savefig('plots/psp_reader_annotation_overlap.pdf')

    # Annotations: REACH vs. Sparser
    plt.figure()
    venn2((reach_annots, sparser_annots),
          set_labels=('REACH', 'Sparser'))
    plt.savefig('plots/reach_sparser_annotation_overlap.pdf')

    # Read the kinases list
    kinases = [r[0] for r in read_unicode_csv('data/kinases.tsv', skiprows=1,
                                              delimiter='\t')]
    psp_kin = set([s for s in ks_annots if s[3] in kinases])
    reader_kin = set([s for s in reader_annots if s[3] in kinases])
    reach_kin = set([s for s in reach_annots if s[3] in kinases])
    sparser_kin = set([s for s in sparser_annots if s[3] in kinases])

    # Kinases: PSP vs. Readers
    plt.figure()
    venn2((psp_kin, reader_kin),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser'))
    plt.savefig('plots/psp_reader_kinase_overlap.pdf')

    # Kinases: REACH vs. Sparser
    plt.figure()
    venn2((reach_kin, sparser_kin),
          set_labels=('REACH', 'Sparser'))
    plt.savefig('plots/reach_sparser_kinase_overlap.pdf')

