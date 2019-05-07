import sys
import itertools
import numpy as np
import pandas as pd
from protmapper import ProtMapper
from protmapper_paper.cptac.brca import get_sites, up_for_hgnc


pd.set_option('display.max_colwidth', 200)
pd.set_option('expand_frame_repr', False)
pd.set_option('display.width', 500)


def get_analysis_conditions():
    # sm = sites mapped
    # dm = data mapped
    # sources = [psp_only, psp_dbs, all]
    sm_options = [False, True]
    dm_options = [False, True]
    source_options = ['psp_only', 'psp_dbs', 'all']
    return list(itertools.product(sm_options, dm_options, source_options))


def build_annot(ann_df, sm_opt, source_opt):
    # Build an index of sites by original position, linked to controllers
    ann_cts_by_site = {}
    if source_opt == 'psp_only':
        ann_source = ann_df[ann_df.SOURCE == 'psp']
    elif source_opt == 'psp_dbs':
        ann_source = ann_df[(ann_df.SOURCE != 'rlimsp') &
                            (ann_df.SOURCE != 'reach') &
                            (ann_df.SOURCE != 'signor')]
    elif source_opt == 'all':
        ann_source = ann_df
    else:
        raise ValueError("Invalid source_opt")

    def get_annot_counts_for_df(df, group_key):
        for group_key, group_df in df.groupby(group_key):
            num_ctrls = len(group_df['CTRL_ID'].unique())
            ann_cts_by_site[group_key] = num_ctrls
        return ann_cts_by_site

    if sm_opt is True:
        # Use mapped sites
        valid_sites = ann_source[ann_source.VALID == True]
        valid_cts = get_annot_counts_for_df(valid_sites,
                               ['GENE_NAME', 'ORIG_RES', 'ORIG_POS'])
        mapped_sites = ann_source[~ann_source.MAPPED_POS.isna()]
        mapped_cts = get_annot_counts_for_df(mapped_sites,
                               ['GENE_NAME', 'MAPPED_RES', 'MAPPED_POS'])
        ann_cts = valid_cts
        ann_cts.update(mapped_cts)
    else:
        valid_sites = ann_source[ann_source.VALID == True]
        valid_cts = get_annot_counts_for_df(valid_sites,
                               ['GENE_NAME', 'ORIG_RES', 'ORIG_POS'])
        ann_cts = valid_cts
    return ann_cts


def build_data(data_sites, dm_opt):
    pm = ProtMapper()
    filt_sites = []
    for site in data_sites:
        (refseq, gene, res, pos, pep, respos) = site
        hgnc_name, up_id = up_for_hgnc(gene)
        if dm_opt is True and up_id is not None:
            ms = pm.map_peptide_to_human_ref(up_id, 'uniprot', pep,
                                             int(respos) + 1)
            if ms.valid and ms.mapped_res and ms.mapped_pos:
                res = ms.mapped_res
                pos = ms.mapped_pos
            site = (refseq, hgnc_name, res, pos, pep, respos)
        filt_sites.append(site)
    return filt_sites


def count_annotations(annot, data, sm_opt, dm_opt, source_opt):
    results = {}
    for refseq, gene, res, pos, pep, respos in data:
        site_key = (gene, res, pos)
        if site_key in annot:
            results[site_key] = annot[site_key]
    ctrl_counts = list(results.values())
    median_ctrls = np.median(ctrl_counts)
    mean_ctrls = np.mean(ctrl_counts)
    return {'num_sites': len(results), 'mean_ctrls': mean_ctrls,
            'median_ctrls': median_ctrls}


if __name__ == '__main__':
    annot_file = sys.argv[1]
    datafile = sys.argv[2]
    #datafile = ('../../data/CPTAC2_Breast_Prospective_Collection_BI_'
    #            'Phosphoproteome.phosphosite.tmt10.tsv')
    # Load the data
    data_sites = get_sites(datafile)
    # Load the annotations
    ann_df = pd.read_csv(annot_file,
                         dtype={'ORIG_POS': str, 'MAPPED_POS': str})
    # Filter out any errors
    ann_df = ann_df[ann_df.ERROR_CODE.isna()]

    # Get the different types of analysis conditions
    conditions = get_analysis_conditions()
    results = []
    for sm_opt, dm_opt, source_opt in conditions:
        print(sm_opt, dm_opt, source_opt)
        ann = build_annot(ann_df, sm_opt, source_opt)
        data = build_data(data_sites, dm_opt)
        count = count_annotations(ann, data, sm_opt, dm_opt, source_opt)
        result = (sm_opt, dm_opt, source_opt, count['num_sites'],
                  count['mean_ctrls'], count['median_ctrls'])
        results.append(result)
    result_df = pd.DataFrame.from_records(results, columns=
                    ['SITES_MAPPED', 'DATA_MAPPED', 'SOURCES', 'NUM_SITES',
                     'MEAN_CTRL_COUNT', 'MEDIAN_CTRL_COUNT'])
    print(results)
