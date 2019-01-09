import numpy as np
import pandas as pd


def get_ipi_mappings(ipi_file):
    names = ('DB_NAME', 'DB_ID', 'IPI_ID', 'UP_IDS', 'MGI_ID')
    df = pd.read_csv(ipi_file, delimiter='\t', skiprows=1,
                     usecols=(0, 1, 2, 4, 10), names=names)
    return df


if __name__ == '__main__':
    col_names = [
       'ProteinIpiSequest', 'GeneSymbol', 'Sequence',
       'Residue', 'Site', 'Ascore', 'NumTissues', 'Tissues', 'TotalCount',
       'GeneIds', 'Entropy', 'MotifPeptide', 'RedundancyCount',
       'NonredundantIpi', 'MatchingProteins']
    #df = get_ipi_mappings('data/ipi.MOUSE.xrefs')
    data_file = 'data/cell5459mmc2.txt'
    df = pd.read_csv(data_file, delimiter='\t',
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



