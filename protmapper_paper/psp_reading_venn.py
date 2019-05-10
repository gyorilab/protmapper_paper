import pandas
import matplotlib
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3
from indra.statements import Agent
from indra.tools.expand_families import Expander
from indra.databases import hgnc_client, uniprot_client

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

expander = Expander()


def get_source_annots(df, sources):
    # Filter to sources
    df = df[df.SOURCE.isin(sources)]
    annot_sites = {}
    for idx, row in df.iterrows():
        if row.VALID:
            annot_sites[(row.CTRL_ID, row.CTRL_NS, row.UP_ID, row.ORIG_RES,
                         row.ORIG_POS)] = row.CTRL_FREQ
        elif row.MAPPED_POS and row.MAPPED_RES:
            annot_sites[(row.CTRL_ID, row.CTRL_NS, row.MAPPED_ID,
                         row.MAPPED_RES, row.MAPPED_POS)] = row.CTRL_FREQ
    return annot_sites


def get_source_sites(df, sources):
    # Filter to sources
    df = df[df.SOURCE.isin(sources)]
    sites = {}
    for idx, row in df.iterrows():
        # Note that we use ORIG_RES and ORIG_POS here just because that gives
        # us unique sites whether or not the sites were mapped
        if row.VALID:
            sites[(row.UP_ID, row.ORIG_RES, row.ORIG_POS)] = row.FREQ
        elif row.MAPPED_POS and row.MAPPED_RES:
            sites[(row.MAPPED_ID, row.MAPPED_RES, row.MAPPED_POS)] = row.FREQ
    return sites


def get_venn_dict_weighted(sites1, sites2):
    venn_dict = {'10': 0, '01': 0, '11': 0}
    for site, count in sites1.items():
        if site in sites2:
            venn_dict['11'] += count
        else:
            venn_dict['10'] += count
    for site, count in sites2.items():
        if site not in sites1:
            venn_dict['01'] += count
    return venn_dict


def get_venn_dict_unweighted(site_dicts):
    return [set(sites.keys()) for sites in site_dicts]


def filter_all_annots(df):
    # Filter the annotations table
    # Filter out all the rows with error code, typically for missing site or
    # residue
    df = df[df.ERROR_CODE.isna()]
    # Keep only rows where the site is VALID or a mapping was found
    df = df[df.DESCRIPTION != 'NO_MAPPING_FOUND']
    # Filter to protein controllers only
    df = df[df.CTRL_IS_PROTEIN == True]
    # Filter to human protein substrates only
    df = df[df.apply(lambda x: uniprot_client.is_human(x['UP_ID']), axis=1)]
    return df


def filter_kinase_annots(annot_sites, include_fplx=True):
    kinase_sites = {}
    for k, v in annot_sites.items():
        ctrl_id, ctrl_ns, _, _, _ = k
        if ctrl_ns == 'HGNC':
            # If genes with HGNC IDs aren't known to be kinases, they will
            # be filtered out here
            if hgnc_client.is_kinase(ctrl_id):
                kinase_sites[k] = v
        elif include_fplx and ctrl_ns == 'FPLX':
            children = expander.get_children(Agent(ctrl_id,
                                                   db_refs={'FPLX': ctrl_id}))
            for _, hgnc_name in children:
                if hgnc_client.is_kinase(hgnc_name):
                    kinase_sites[k] = v
                    break
        # The rest of the entries here typically have UP IDs that correspond
        # to non-human proteins or aren't proteins at all.
    return kinase_sites


def filter_sites(df):
    # Filter out all the rows with error code, typically for missing site or
    # residue
    df = df[df.ERROR_CODE.isna()]
    # Keep only rows where the site is VALID or a mapping was found
    df = df[df.DESCRIPTION != 'NO_MAPPING_FOUND']
    # Filter to human protein substrates only
    df = df[df.apply(lambda x: uniprot_client.is_human(x['UP_ID']), axis=1)]
    return df


def print_reading_contribs(reader_sites, psp_sites):
    reader_only = set(reader_sites.keys()) - set(psp_sites.keys())
    for ctrl_id, ctrl_ns, up_id, residue, pos in \
            sorted(reader_only, key=lambda x: (x[0], x[2], x[4])):
        target_name = uniprot_client.get_gene_name(up_id, web_fallback=False)
        if target_name is None:
            print('Could not get gene name for %s' % up_id)
        print('%s -> %s-%s%s' % (ctrl_id, target_name, residue, int(pos)))

if __name__ == '__main__':
    # SITES
    df = pandas.read_csv('output/site_info.csv')
    dfs = filter_sites(df)
    psp_sites = get_source_sites(dfs, ['psp'])
    reach_sites = get_source_sites(dfs, ['reach'])
    sparser_sites = get_source_sites(dfs, ['sparser'])
    rlimsp_sites = get_source_sites(dfs, ['rlimsp'])
    reader_sites = get_source_sites(dfs, ['reach', 'sparser', 'rlimsp'])

    plt.ion()
    # Sites: PSP vs. Readers
    plt.figure()
    venn2(get_venn_dict_unweighted([psp_sites, reader_sites]),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser/RLIMS-P'))
    plt.savefig('plots/psp_reader_sites_overlap_distinct.pdf')

    # Sites: REACH vs. Sparser vs. RLIMS-P
    plt.figure()
    venn3(get_venn_dict_unweighted([reach_sites, sparser_sites, rlimsp_sites]),
          set_labels=('REACH', 'Sparser', 'RLIMS-P'))
    plt.savefig('plots/reader_sites_overlap_distinct.pdf')

    # ANNOTATIONS
    # Load the overall annotations table
    df = pandas.read_csv('output/annotations.csv')

    dfa = filter_all_annots(df)
    psp_annots = get_source_annots(dfa, ['psp'])
    db_annots = get_source_annots(dfa, ['psp', 'hprd', 'signor', 'pid',
                                        'reactome', 'bel'])
    reach_annots = get_source_annots(dfa, ['reach'])
    sparser_annots = get_source_annots(dfa, ['sparser'])
    rlimsp_annots = get_source_annots(dfa, ['rlimsp'])
    reader_annots = get_source_annots(dfa, ['reach', 'sparser', 'rlimsp'])

    # Annotations: PSP vs. Readers
    plt.figure()
    venn2(get_venn_dict_unweighted([psp_annots, reader_annots]),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser/RLIMS-P'))
    plt.savefig('plots/psp_reader_annotation_overlap_distinct.pdf')

    # Annotations: REACH vs. Sparser vs. RLIMS-P
    plt.figure()
    venn3(get_venn_dict_unweighted([reach_annots, sparser_annots,
                                    rlimsp_annots]),
          set_labels=('REACH', 'Sparser', 'RLIMS-P'))
    plt.savefig('plots/reader_annotation_overlap_distinct.pdf')

    # Kinases only
    psp_annotsk = filter_kinase_annots(psp_annots)
    db_annotsk = filter_kinase_annots(db_annots)
    reach_annotsk = filter_kinase_annots(reach_annots)
    sparser_annotsk = filter_kinase_annots(sparser_annots)
    rlimsp_annotsk = filter_kinase_annots(rlimsp_annots)
    reader_annotsk = filter_kinase_annots(reader_annots)

    # Annotations: PSP vs. Readers
    plt.figure()
    venn2(get_venn_dict_unweighted([psp_annotsk, reader_annotsk]),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser/RLIMS-P'))
    plt.savefig('plots/psp_reader_annotation_overlap_distinct_kinase.pdf')

    # Annotations: REACH vs. Sparser vs. RLIMS-P
    plt.figure()
    venn3(get_venn_dict_unweighted([reach_annotsk, sparser_annotsk,
                                    rlimsp_annotsk]),
          set_labels=('REACH', 'Sparser', 'RLIMS-P'))
    plt.savefig('plots/reader_annotation_overlap_distinct_kinase.pdf')

    # Non-FamPlex Kinases only
    psp_annotskn = filter_kinase_annots(psp_annots, False)
    db_annotskn = filter_kinase_annots(db_annots, False)
    reach_annotskn = filter_kinase_annots(reach_annots, False)
    sparser_annotskn = filter_kinase_annots(sparser_annots, False)
    rlimsp_annotskn = filter_kinase_annots(rlimsp_annots, False)
    reader_annotskn = filter_kinase_annots(reader_annots, False)

    # Annotations: PSP vs. Readers
    plt.figure()
    venn2(get_venn_dict_unweighted([psp_annotskn, reader_annotskn]),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser/RLIMS-P'))
    plt.savefig('plots/psp_reader_annotation_overlap_distinct_kinase_nofamplex.pdf')

    # Annotations: DBs vs. Readers
    plt.figure()
    venn2(get_venn_dict_unweighted([db_annotskn, reader_annotskn]),
          set_labels=('Databases', 'REACH/Sparser/RLIMS-P'))
    plt.savefig('plots/db_reader_annotation_overlap_distinct_kinase_nofamplex.pdf')

    # Annotations: REACH vs. Sparser vs. RLIMS-P
    plt.figure()
    venn3(get_venn_dict_unweighted([reach_annotskn, sparser_annotskn,
                                    rlimsp_annotskn]),
          set_labels=('REACH', 'Sparser', 'RLIMS-P'))
    plt.savefig('plots/reader_annotation_overlap_distinct_kinase_nofamplex.pdf')

    print_reading_contribs(reader_annotskn, psp_annotskn)
