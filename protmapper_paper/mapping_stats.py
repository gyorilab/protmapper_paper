"""Script to generate mapping statistics for CPTAC data sets."""
import re
import sys
import numpy as np
import pandas as pd
from protmapper import ProtMapper
from protmapper import uniprot_client
from indra.databases import hgnc_client, uniprot_client


def get_sites(datafile):
    """Return a list of sites for a given data set; each site is a tuple."""
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

    df = pd.read_csv(datafile, delimiter='\t')

    site_gene = df[['Phosphosite', 'Gene', 'Peptide']]

    sites = []
    for rs_site_text, gene, peptides in site_gene.values:
        refseq, site_text = rs_site_text.split(':')
        sites_for_id = parse_sites(site_text, peptides)
        for site, respos, peptide in sites_for_id:
            assert site[0] == peptide[respos]
            sites.append((refseq, gene, site[0], site[1], peptide, respos))

    return sites


def up_for_hgnc(gene):
    """Return HGNC symbol and UniProt ID for a potentially outdated gene
    name."""
    hgnc_id = hgnc_client.get_current_hgnc_id(gene)
    if hgnc_id is None:
        #print("Couldn't find current HGNC ID for gene %s" % gene)
        hgnc_name = gene
        up_id = None
    elif isinstance(hgnc_id, list):
        #print("More than one HGNC ID for gene %s" % gene)
        hgnc_name = gene
        up_id = None
    else:
        hgnc_name = hgnc_client.get_hgnc_name(hgnc_id)
        up_id_str = hgnc_client.get_uniprot_id(hgnc_id)
        if up_id_str is None:
            #print("No Uniprot ID for HGNC ID %s, gene %s" % (hgnc_id, gene))
            up_id = None
        elif ',' in up_id_str:
            up_ids = [u.strip() for u in up_id_str.split(',')]
            up_id = up_ids[0]
        else:
            up_id = up_id_str
    return hgnc_name, up_id


def up_id_for_rs(refseq):
    """Return a UniProt ID for a RefSeq ID."""
    up_ids = uniprot_client.get_ids_from_refseq(refseq, reviewed_only=True)
    if not up_ids:
        up_id = None
    elif len(up_ids) > 1:
        #print("More than one up id for rs: up %s, rs %s" %
        #        (str(up_ids), refseq))
        up_id = up_ids[0]
    else:
        up_id = up_ids[0]
    return up_id


def map_peptide(up_id, res, pos, pep, respos, suffix, pm):
    """Map a given peptide to reference and return dict of mapping status."""
    result = {}
    if up_id is None:
        result['mappable_%s' % suffix] = None
        result['valid_%s' % suffix] = None
    else:
        is_valid = uniprot_client.verify_location(up_id, res, pos)
        base_id = up_id.split('-')[0]
        ms = pm.map_peptide_to_human_ref(base_id, 'uniprot', pep,
                                         int(respos) + 1)
        result['valid_%s' % suffix] = is_valid
        if ms.valid and ms.mapped_res and ms.mapped_pos:
            result['mappable_%s' % suffix] = True
        else:
            result['mappable_%s' % suffix] = False
    return result


def iso_specific(up_id):
    """Return True if the given UniProt ID corresponds to the -1 isoform."""
    return ('-' in up_id and up_id.split('-')[1] != '1')


def valid_counts(sitelist):
    results = []
    pm = ProtMapper()
    for ix, (refseq, gene, res, pos, pep, respos) in enumerate(sitelist):
        result = {'refseq': refseq, 'gene': gene, 'res': res, 'pos': pos,
                  'peptide': pep, 'respos': respos}
        if ix % 10000 == 0:
            print(ix)
        up_rs = up_id_for_rs(refseq)
        hgnc_name, up_hgnc = up_for_hgnc(gene)
        result['up_hgnc'] = up_hgnc
        result['up_rs'] = up_rs
        if up_rs is None:
            result['up_rs_iso_specific'] = None
        else:
            result['up_rs_iso_specific'] = iso_specific(up_rs)
        hgnc_map_result = map_peptide(up_hgnc, res, pos, pep, respos, 'hgnc',
                                      pm)
        rs_map_result = map_peptide(up_rs, res, pos, pep, respos, 'rs',
                                    pm)
        result.update(hgnc_map_result)
        result.update(rs_map_result)
        results.append(result)

    df = pd.DataFrame.from_dict(results, orient='columns')
    return df


def print_valid_stats(df):
    no_up_id_hgnc = df[df.up_hgnc.isna()]
    has_up_id_hgnc = df[~df.up_hgnc.isna()]
    no_up_id_rs = df[df.up_rs.isna()]
    has_up_id_rs = df[~df.up_rs.isna()]
    iso_spec = df[df.up_rs_iso_specific == True]
    print()
    print("Total Sites: %d" % len(df))
    print("No Uniprot ID from HGNC: %d" % len(no_up_id_hgnc))
    print("No Uniprot ID from RS: %d" % len(no_up_id_rs))
    print("Isoform-specific RS ID: %d (%.1f)" % (len(iso_spec),
                                           100 * len(iso_spec) / len(df)))
    print("Isoform-specific RS but mappable: %d" %
            (len(iso_spec[iso_spec['mappable_hgnc'] == True])))
    print("No Uniprot ID from RS but mappable: %d"
            % (len(no_up_id_rs[no_up_id_rs.mappable_hgnc == True])))
    for id in ('hgnc', 'rs'):
        valid = df[df['valid_%s' % id] == True]
        invalid = df[df['valid_%s' % id] == False]
        vm = valid[valid['mappable_%s' % id] == True]
        vnm = valid[valid['mappable_%s' % id] == False]
        ivm = invalid[invalid['mappable_%s' % id] == True]
        ivnm = invalid[invalid['mappable_%s' % id] == False]
        print("%s Valid: %d" % (id, len(valid)))
        print("%s Invalid: %d (%.1f)" %
            (id, len(invalid), (100 * len(invalid) / len(df))))
        print("%s Valid mappable: %d" % (id, len(vm)))
        print("%s Valid not mappable: %d" % (id, len(vnm)))
        print("%s Invalid mappable: %d" % (id, len(ivm)))
        print("%s Invalid not mappable: %d" % (id, len(ivnm)))
        print("%s total mappable: %d" %
                (id, len(df[df['mappable_%s' % id] == True])))


if __name__ == '__main__':
    # Args
    datafile = sys.argv[1]
    output_file = sys.argv[2]

    # Process the dataset
    sites = get_sites(datafile)

    # Get counts of valid and mappable sites
    site_results = valid_counts(sites)

    # Save dataframe as CSV
    site_results.to_csv(output_file)

    # Print results
    print_valid_stats(site_results)

