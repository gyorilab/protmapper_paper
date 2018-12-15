import pickle
from collections import Counter, defaultdict
from indra.databases import uniprot_client
from indra.util import read_unicode_csv, write_unicode_csv
from indra_db.client import get_statements_by_gene_role_type
from sitemapper import SiteMapper


def get_db_phos_stmts(reload=False, filename='output/db_phos_stmts.pkl'):
    if reload is False:
        with open(filename, 'rb') as f:
            phos_stmts = pickle.load(f)
    else:
        phos_stmts = get_statements_by_gene_role_type(
                            stmt_type='Phosphorylation', fix_refs=False,
                            preassembled=False,
                            with_evidence=True, with_support=False)
        with open(filename, 'wb') as f:
            pickle.dump(phos_stmts, f)
    return phos_stmts


def get_db_sites(reload=False, filename='output/gene_sites.csv',
                 phos_stmts=None):
    if reload is False:
        return list(set([tuple(r) for r in read_unicode_csv(filename)]))
    else:
        if phos_stmts is None:
            raise ValueError('phos_stmts must be provided.')

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
            gene_sites.append((up_id, gene_name, up_mnemonic, s.residue,
                               s.position))

        write_unicode_csv(filename, gene_sites)
        return gene_sites


def site_cache_stats():
    sm = SiteMapper(use_cache=True)
    ms_desc_list = []
    examples = defaultdict(list)
    for site_key, mapped_site in sm._cache.items():
        if mapped_site is not None:
            ms_desc_list.append(mapped_site.description)
            if mapped_site.description in (
                           'INFERRED_ALTERNATIVE_ISOFORM',
                           'INFERRED_MOUSE_SITE',
                           'INFERRED_RAT_SITE',
                           'INFERRED_METHIONINE_CLEAVAGE',
                           'NO_MAPPING_FOUND'):
                examples[mapped_site.description].append(mapped_site)
        else:
            ms_desc_list.append('NO_MAPPING_FOUND')

    ctr = Counter(ms_desc_list)
    ctr = sorted([(k, v) for k, v in ctr.items()],
                 key=lambda x: x[1], reverse=True)
    return ctr, examples



if __name__ == '__main__':
    """
    phos_stmts = get_db_phos_stmts(reload=True)
    gene_sites = get_db_sites(reload=True, phos_stmts=phos_stmts)
    """
    """
    gene_sites = get_db_sites(reload=False)
    import random
    random.seed(1)
    random.shuffle(gene_sites)
    gene_sites = gene_sites[0:1000]
    sm = SiteMapper(use_cache=True)
    site_list = [(gs[0], 'uniprot', gs[3], gs[4]) for gs in gene_sites]
    mapped_sites = sm.map_sitelist_to_human_ref(site_list)
    with open('sm_cache.pkl', 'wb') as f:
        import pickle
        pickle.dump(sm._cache, f)
    """
    site_ctr, examples = site_cache_stats()

