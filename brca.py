import re
import numpy as np
import pandas as pd
from protmapper import ProtMapper
from protmapper import uniprot_client
from indra.databases import hgnc_client, uniprot_client



def get_sites(df):
    def get_max_peptide(seq):
        peptides = seq.split(';')
        max_pep_ix = np.argmax([len(p) for p in peptides])
        return peptides[max_pep_ix]

    def parse_sites(site_text, peptides):
        rest = site_text
        sites = []
        while rest:
            m = re.match('([sty]\d+)', rest)
            if m:
                one_site = m.groups()[0]
                res = one_site[0]
                pos = one_site[1:]
                rest = rest[len(one_site):]
                sites.append((res.upper(), pos))
        # Get peptide positions of
        peptide = get_max_peptide(peptides)
        positions = []
        for ix, c in enumerate(peptide):
            if c.islower():
                positions.append(ix)
        assert len(positions) == len(sites)
        return list(zip(sites, positions, [peptide.upper()] * 3))


    site_gene = df[['Phosphosite', 'Gene', 'Peptide']]

    sites = []
    for rs_site_text, gene, peptides in site_gene.values:
        refseq, site_text = rs_site_text.split(':')
        sites_for_id = parse_sites(site_text, peptides)
        for site, respos, peptide in sites_for_id:
            assert site[0] == peptide[respos]
            sites.append((refseq, gene, site[0], site[1], peptide, respos))

    hgnc_sitelist = []
    for refseq, gene, res, pos, pep, respos in sites:
        hgnc_id = hgnc_client.get_current_hgnc_id(gene)
        if hgnc_id is None:
            print("Couldn't find current HGNC ID for gene %s" % gene)
            hgnc_name = gene
            up_id = None
        elif isinstance(hgnc_id, list):
            print("More than one HGNC ID for gene %s" % gene)
            hgnc_name = gene
            up_id = None
        else:
            hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
            up_id_str = hgnc_client.get_uniprot_id(hgnc_id)
            if up_id_str is None:
                print("No Uniprot ID for HGNC ID %s, gene %s" % (hgnc_id, gene))
                up_id = None
            elif ',' in up_id_str:
                up_ids = [u.strip() for u in up_id_str.split(',')]
                up_id = up_ids[0]
            else:
                up_id = up_id_str
        hgnc_sitelist.append((hgnc_name, up_id, 'uniprot', res, pos, pep, respos))

    refseq_sitelist = []
    for refseq, gene, res, pos, pep, respos in sites:
        up_ids = uniprot_client.get_ids_from_refseq(refseq, reviewed_only=True)
        if not up_ids:
            up_id = None
        elif len(up_ids) > 1:
            print("More than one up id for rs: up %s, rs %s" %
                    (str(up_ids), refseq))
            up_id = up_ids[0]
        else:
            up_id = up_ids[0]
        refseq_sitelist.append((refseq, up_id, 'uniprot', res, pos, pep, respos))

    return (hgnc_sitelist, refseq_sitelist)


def valid_counts(sitelist, id_type, pm):
    no_up_id = []
    iso_specific = []
    valid = []
    mappable = []
    for ix, (ext_id, up_id, _, res, pos, pep, respos) in \
                                                enumerate(sitelist):
        if ix % 10000 == 0:
            print(ix)
        if up_id is None:
            no_up_id.append(True)
            iso_specific.append(False)
            valid.append(False)
            mappable.append(False)
        elif up_id is not None:
            no_up_id.append(False)
            if '-' in up_id and up_id.split('-')[1] != '1':
                iso_specific.append(True)
            else:
                iso_specific.append(False)
            is_valid = uniprot_client.verify_location(up_id, res, pos)
            valid.append(is_valid)
            base_id = up_id.split('-')[0]
            ms = pm.map_peptide_to_human_ref(base_id, 'uniprot', pep,
                                             int(respos) + 1)
            if not ms.valid and is_valid:
                #import ipdb; ipdb.set_trace()
                pass
            if ms.valid and ms.mapped_res and ms.mapped_pos:
                mappable.append(True)
            else:
                mappable.append(False)
    ext_ids, up_ids, _, res, pos, pep, respos = list(zip(*sitelist))
    ext_id_col = 'refseq' if id_type == 'refseq' else 'hgnc'
    data_dict = {ext_id_col: ext_ids, 'up_id': up_ids, 'res': res, 'pos': pos,
                 'no_up_id': no_up_id, 'iso_specific': iso_specific,
                 'valid': valid, 'peptide': pep, 'respos': respos,
                 'mappable': mappable}
    df = pd.DataFrame.from_dict(data_dict, orient='columns')
    return df


def print_valid_stats(df):
    print("Total Sites: %d" % len(df))
    no_up_id = len(df[df.no_up_id == True])
    has_up = df[df.no_up_id == False]
    print("No Uniprot ID: %d" % no_up_id)
    print("Isoform-specific ID: %d" % len(has_up[has_up.iso_specific == True]))
    valid = has_up[has_up.valid == True]
    invalid = has_up[has_up.valid == False]
    vm = valid[valid.mappable == True]
    vnm = valid[valid.mappable == False]
    ivm = invalid[invalid.mappable == True]
    ivnm = invalid[invalid.mappable == False]
    print("Valid: %d" % len(valid))
    print("Invalid: %d (%.1f)" % (len(invalid), (100 * len(invalid) / len(df))))
    print("Valid mappable: %d" % len(vm))
    print("Valid not mappable: %d" % len(vnm))
    print("Invalid mappable: %d" % len(ivm))
    print("Invalid not mappable: %d" % len(ivnm))

# Using HGNC IDs

# Stats of sites themselves (before mapping to)

# Build up list of uniprot IDs with sites

# Provide dict of sites


if __name__ == '__main__':
    df = pd.read_csv('data/CPTAC2_Breast_Prospective_Collection_BI_'
                     'Phosphoproteome.phosphosite.tmt10.tsv', delimiter='\t')

    hgnc_sitelist, refseq_sitelist = get_sites(df)
    pm = ProtMapper()
    print("-- HGNC --")
    hgnc_valid = valid_counts(hgnc_sitelist, 'hgnc', pm)
    print_valid_stats(hgnc_valid)
    print("-- RefSeq ID --")
    rs_valid = valid_counts(refseq_sitelist, 'refseq', pm)
    print_valid_stats(rs_valid)

    #hgnc_valid = valid_counts(hgnc_sitelist, pm)

"""
sites = df[~df.Site.isna()].Site
site_parse = []
for site_text in sites:
    gene, site_raw = site_text.rsplit('-', maxsplit=1)
    hgnc_id = hgnc_client.get_hgnc_id(gene)
    up_id_str = hgnc_client.get_uniprot_id(hgnc_id)
    if up_id_str is None:
        print("No Uniprot ID for %s" % gene)
        continue
    if ',' in up_id_str:
        up_ids = [u.strip() for u in up_id_str.split(',')]
        up_id = up_ids[0]
    else:
        up_id = up_id_str

    while site_raw:
        m = re.match('([STY]\d+[sty])', site_raw)
        if m:
            one_site = m.groups()[0]
            res = one_site[0]
            pos = one_site[1:-1]
            site_raw = site_raw[len(one_site):]
            site_parse.append((up_id, 'uniprot', res, pos))

pm = ProtMapper()

mapped = pm.map_sitelist_to_human_ref(site_parse)
"""
