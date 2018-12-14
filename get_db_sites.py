import pickle
from indra.util import read_unicode_csv, write_unicode_csv
from indra.databases import uniprot_client
from indra_db.client import get_statements_by_gene_role_type

def get_db_sites(reload=False, filename='output/gene_sites.csv'):
    if reload is False:
        return list(read_unicode_csv(filename))
    else:
        phos_stmts = get_statements_by_gene_role_type(
                            stmt_type='Phosphorylation', fix_refs=False,
                            with_evidence=False, with_support=False)

        phos_filename = 'output/db_phos_stmts.pkl'

        with open(phos_filename, 'wb') as f:
            pickle.dump(phos_stmts, f)

        filter_no_enzyme = [s for s in phos_stmts
                            if s.agent_list()[0] is not None]

        gene_sites = []
        for ix, s in enumerate(filter_no_enzyme):
            if ix % 1000 == 0:
                print(ix)
            up_id = s.sub.db_refs.get('UP')
            if not up_id:
                continue
            up_mnemonic = uniprot_client.get_mnemonic(up_id)
            if not s.residue or not s.position:
                continue
            gene_name = uniprot_client.get_gene_name(up_id)
            gene_sites.append([up_id, gene_name, up_mnemonic, s.residue,
                               s.position])

        write_unicode_csv(filename, gene_sites)
        return gene_sites

if __name__ == '__main__':
    gene_sites = get_db_sites(reload=False)

