import pandas
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3


def get_source_annots(df, sources):
    # Filter to sources
    df = df[df.SOURCE.isin(sources)]
    annot_sites = {}
    for idx, row in df.iterrows():
        # Note that we use ORIG_RES and ORIG_POS here just because that gives
        # us unique sites whether or not the sites were mapped
        annot_sites[(row.CTRL_ID, row.UP_ID, row.ORIG_RES, row.ORIG_POS)] = \
            row.CTRL_FREQ
    return annot_sites


def get_source_sites(df, sources):
    # Filter to sources
    df = df[df.SOURCE.isin(sources)]
    sites = {}
    for idx, row in df.iterrows():
        sites[(row.UP_ID, row.ORIG_RES, row.ORIG_POS)] = row.FREQ
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
    return df


def filter_kinase_annots(df):
    # Filter the annotations table
    # Filter out all the rows with error code, typically for missing site or
    # residue
    df = df[df.ERROR_CODE.isna()]
    # Keep only rows where the site is VALID or a mapping was found
    df = df[df.DESCRIPTION != 'NO_MAPPING_FOUND']
    # Filter to protein controllers only
    df = df[df.CTRL_IS_PROTEIN == True]
    return df


def filter_sites(df):
    # Filter out all the rows with error code, typically for missing site or
    # residue
    df = df[df.ERROR_CODE.isna()]
    # Keep only rows where the site is VALID or a mapping was found
    df = df[df.DESCRIPTION != 'NO_MAPPING_FOUND']
    return df


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
    plt.savefig('plots/psp_reader_sites_overlap_distincs.pdf')

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
    reach_annots = get_source_annots(dfa, ['reach'])
    sparser_annots = get_source_annots(dfa, ['sparser'])
    rlimsp_annots = get_source_annots(dfa, ['rlimsp'])
    reader_annots = get_source_annots(dfa, ['reach', 'sparser', 'rlimsp'])

    # Annotations: PSP vs. Readers
    plt.figure()
    venn2(get_venn_dict_unweighted([psp_annots, reader_annots]),
          set_labels=('PhosphoSitePlus', 'REACH/Sparser/RLIMS-P'))
    plt.savefig('plots/psp_reader_annotation_overlap_distincs.pdf')

    # Annotations: REACH vs. Sparser vs. RLIMS-P
    plt.figure()
    venn3(get_venn_dict_unweighted([reach_annots, sparser_annots,
                                    rlimsp_annots]),
          set_labels=('REACH', 'Sparser', 'RLIMS-P'))
    plt.savefig('plots/reader_annotation_overlap_distinct.pdf')

