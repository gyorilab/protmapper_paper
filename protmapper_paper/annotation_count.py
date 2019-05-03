import sys
import itertools
import pandas as pd
from protmapper_paper.brca import get_sites

pd.set_option('display.max_colwidth', 200)
pd.set_option('expand_frame_repr', False)
pd.set_option('display.width', 500)


def get_analysis_conditions():
    # sm = sites mapped
    # dm = data mapped
    # sources = [psp_only, psp_dbs, all]
    sm_options = [True, False]
    dm_options = [True, False]
    source_options = ['psp_only', 'psp_dbs', 'all']
    return list(itertools.product(sm_options, dm_options, source_options))


def build_annot(ann_df, sm_opt):
    # Build an index of sites by original position, linked to controllers
    ann_cts_by_site = {}
    if True: # sm_opt is False:
        for group_key, group_df in \
                     ann_df.groupby('GENE_NAME', 'ORIG_RES', 'ORIG_POS'):
            num_ctrls = len(group_df['CTRL_ID'].unique())
            ann_cts_by_site[group_key] = num_ctrls
    return ann_cts_by_site


def count_annotations(annot, sm_opt, dm_opt, source_opt):
    pass


if __name__ == '__main__':
    annot_file = sys.argv[1]
    #datafile = sys.argv[2]
    datafile = ('data/CPTAC2_Breast_Prospective_Collection_BI_'
                'Phosphoproteome.phosphosite.tmt10.tsv')
    # Load the data
    data_df = get_sites(datafile)
    # Load the annotations
    ann_df = pd.read_csv(annot_file, dtype=str)
    # Filter out any errors
    ann_df = ann_df[ann_df.ERROR_CODE.isna()]

    # Get the different types of analysis conditions
    conditions = get_analysis_conditions()

    """
    results = {}
    # Run the analysis for each condition
    for sm_opt, dm_opt, source_opt in conditions:
        count_annotations(ann_df, sm_opt, dm_opt, source_opt)
    """
