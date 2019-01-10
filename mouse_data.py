import pickle
import numpy as np
import pandas as pd
from protmapper import phosphosite_client as pspc

def get_ipi_mappings(ipi_file):
    names = ('DB_NAME', 'DB_ID', 'IPI_ID', 'UP_IDS', 'MGI_ID')
    df = pd.read_csv(ipi_file, delimiter='\t', skiprows=1,
                     usecols=(0, 1, 2, 4, 10), names=names)
    return df


def get_data(filename):
    col_names = [
       'ProteinIpiSequest', 'GeneSymbol', 'Sequence',
       'Residue', 'Site', 'Ascore', 'NumTissues', 'Tissues', 'TotalCount',
       'GeneIds', 'Entropy', 'MotifPeptide', 'RedundancyCount',
       'NonredundantIpi', 'MatchingProteins']
    df = pd.read_csv(filename, delimiter='\t',
                     usecols=(1, 2, 5, 6, 7, 9, 10, 11, 14, 24, 25, 26, 27,
                              28, 29),
                     names=col_names, skiprows=1)

    def get_id_from_string(id_type, outer_delimiter):
        def func(id_string):
            if id_string is np.nan:
                return id_string
            id_list = [id.strip() for id in id_string.split(outer_delimiter)]
            sp_ids = [id.strip().split(':', maxsplit=1)[1] for id in id_list
                      if id.startswith(id_type)]
            if len(sp_ids) > 1:
                print("Warning: multiple IDs: %s" % sp_ids)
            return '|'.join(sp_ids) if sp_ids else np.nan
        return func

    df["MgiId"] = df['GeneIds'].apply(get_id_from_string('MG', ';'))
    df["SwissProtId"] = df['ProteinIpiSequest'].apply(
                                    get_id_from_string('SWISS-PROT', '|'))
    df["TremblId"] = df['ProteinIpiSequest'].apply(
                                    get_id_from_string('TREMBL', '|'))
    return df


def count_annotations(df, site_list):
    site_data = df[['SwissProtId', 'TremblId', 'Residue', 'Site']].values
    site_list_set = set(site_list)
    annot = []
    no_annot = []
    for ix, (sp_id, tr_id_str, res, pos) in enumerate(site_data):
        last_len = len(annot)
        print("site %d: annot %d, no annot %d" %
                (ix, len(annot), len(no_annot)))
        # If both SwissProt and Trembl IDs are NaN, then we're not going to
        # find any annotations
        site_entry = (ix, sp_id, tr_id_str, res, str(pos))
        # Try the Swiss Prot ID
        site_found = False
        if sp_id is not np.nan:
            if (sp_id, res, str(pos)) in site_list_set:
                site_found = True
        if tr_id_str is not np.nan:
            for tr_id in tr_id_str.split(';'):
                if (tr_id, res, str(pos)) in site_list_set:
                    site_found = True
                    break
        # If we got here, it means none of the IDs were found to be annotated
        if site_found:
            annot.append(site_entry)
        else:
            no_annot.append(site_entry)
    return annot, no_annot


def count_kin_sub(df, site_dict):
    site_annot = 0
    total_annot = 0
    site_data = df[['SwissProtId', 'Residue', 'Site']].values
    for up_id, res, pos in site_data:
        site_key = (up_id, res, str(pos))
        if site_key in site_dict:
            site_annot += 1
            num_kinases = len(site_dict[site_key])
            total_annot += num_kinases
        print("annot %d, total annot %d" % (site_annot, total_annot))
    return site_annot, total_annot


if __name__ == '__main__':
    df = get_data('data/cell5459mmc2.txt')
    # Get baseline for the frequency with which the ascribed Uniprot ID plus
    # site and position is known to phosphosite
    psp_sites = pspc.sites_only()
    annot_sites = count_annotations(df, psp_sites)

    #with open('output/psp_relations_by_site.pkl', 'rb') as f:
    #    psp_kin_sub = pickle.load(f)
    #site_annot, total_annot = count_kin_sub(df, psp_kin_sub)
