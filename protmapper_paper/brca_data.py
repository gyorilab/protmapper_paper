import re
import csv
import sys
import pickle
import numpy as np
import pandas as pd
from indra.tools import assemble_corpus as ac
from indra.databases import uniprot_client, hgnc_client
from protmapper import ProtMapper


BRCA_DATA = 'data/breast_phosphosites.txt'
UP_MAPPINGS = 'data/HUMAN_9606_idmapping.dat'
BRCA_MAPPED = 'output/brca_up_mappings.txt'


if __name__ == '__main__':
    if sys.argv[1] == 'map_uniprot':
        # Build Refseq -> Uniprot Mappings
        rs_up = {}
        with open(UP_MAPPINGS, 'rt') as f:
            csvreader = csv.reader(f, delimiter='\t')
            for up_id, other_db, other_id in csvreader:
                if other_db == 'RefSeq':
                    if other_id in rs_up:
                        rs_up[other_id].append(up_id)
                    else:
                        rs_up[other_id] = [up_id]
        # Parse the BRCA data and add UP IDs where available
        brca_data = []
        with open(BRCA_DATA, 'rt') as f:
            lines = f.readlines()
            for ix, line in enumerate(lines[1:]): # Skip the header line
                if ix % 100 == 0:
                    print(ix)
                m = re.match(r'(.*?)\.(.*):([sty])(\d+)', line)
                assert len(m.groups()) == 4
                gene_name, rs_id, res, pos = m.groups()
                res = res.upper()
                # Get Uniprot ID
                up_ids = rs_up.get(rs_id)
                if up_ids is None:
                    print(f"No Uniprot ID for RefSeq %s" % rs_id)
                    valid = None
                    site = (gene_name, rs_id, up_ids, res, pos, valid)
                else:
                    valid = []
                    for up_id in up_ids:
                        valid.append(uniprot_client.verify_location(up_id, res,
                                                                    pos))
                    up_id_str = '|'.join(up_ids)
                    valid_str = '|'.join([str(v) for v in valid])
                    site = (gene_name, rs_id, up_id_str, res, pos, valid_str)
                    if not np.all(valid):
                        print("%d: Site not valid: %s" % (ix, str(site)))
                brca_data.append(site)
        with open(BRCA_MAPPED, 'wt') as f:
            csvwriter = csv.writer(f, delimiter='\t')
            csvwriter.writerows(brca_data)
    elif sys.argv[1] == 'site_stats':
        pm = ProtMapper()
        # Load INDRA statements, sorted by site
        indra_stmts_filename = sys.argv[2]
        output_filename = sys.argv[3]
        with open(indra_stmts_filename, 'rb') as f:
            stmts_by_site = pickle.load(f)
        # Load the BRCA peptides with uniprot IDs
        df = pd.read_csv(BRCA_MAPPED, delimiter='\t')
        no_up_id = 0
        no_valid_up = 0
        has_stmts = 0
        has_mapped_stmts = 0
        has_canonical_iso_count = 0
        iso_id_but_pep_in_canon = 0
        counter = 0
        values = list(df.values)
        #import random
        #random.shuffle(values)
        #values = values # values[0:2000]
        total_sites = len(values)
        for gene_name, rs_id, up_id, res, pos, valid in values:
            counter += 1
            if (counter % 100) == 0:
                print(counter)
            # If there is no uniprot ID, then we can't get the sequence,
            # and can't look the motif in the canonical seq
            if up_id is np.nan:
                no_up_id += 1
                assert valid is np.nan
            elif valid is not np.nan:
                valids = [True if v == 'True' else False
                          for v in valid.split('|')]
                up_ids = up_id.split('|')
                # UP mappings but none with the residue at the given position
                # in the sequence; this actually doesn't seem to happen at all
                if not np.any(valids):
                    no_valid_up += 1
                # At least one UP ID with sequence and matching res/pos
                else:
                    # Is at least one of the UP IDs the same as the canonical
                    # isoform for the gene?
                    has_canonical_iso = False
                    # Are any of the matching UP IDs in the INDRA site set?
                    has_annot_up_id = False
                    # Is the mapped site in the INDRA site?
                    has_annot_map_up_id = False
                    # pep_in_canon
                    pep_in_canon = False
                    # Check all the uniprot IDs
                    valid_up_ids = []
                    hgnc_id = hgnc_client.get_hgnc_id(gene_name)
                    canonical_iso = hgnc_client.get_uniprot_id(hgnc_id)
                    if canonical_iso is None or \
                       len(canonical_iso.split(',')) > 1:
                        continue
                    for u, v in zip(up_ids, valids):
                        # Skip any where the site position is not valid
                        if not v:
                            continue
                        # Make a note of the valid UP ID
                        valid_up_ids.append(u)
                        # Check to see if we have an isoform-specific ID
                        if len(u.split('-')) == 2:
                            base_up, isoform_num = u.split('-')
                        else:
                            base_up = u
                            isoform_num = None
                            # No isoform dash, andc
                            if base_up == canonical_iso:
                                has_canonical_iso = True
                        # Check to see if this ID is the canonical one for
                        # the gene
                        if ((isoform_num is None and base_up == canonical_iso)
                            or (isoform_num == '1' and
                                                base_up == canonical_iso)):
                            has_canonical_iso = True
                        # No canonical isoform in matching set of IDs; check
                        # to see if the peptide is present in the canonical
                        # isoform
                        if (base_up, res, str(pos)) in stmts_by_site:
                            if isoform_num is not None and isoform_num != '1':
                                print("Site matches indra even though "
                                      "isoform is not 1", u)
                            has_annot_up_id = True
                    if has_canonical_iso == True:
                        has_canonical_iso_count += 1
                    else:
                        # For now, just take first valid ID
                        for valid_id in valid_up_ids:
                            motif, site_pos = pm.motif_from_position(valid_id,
                                                                     pos)
                            ms = pm.map_peptide_to_human_ref(
                                    gene_name, 'hgnc', motif, site_pos)
                            if ms.mapped_pos is not None:
                                pep_in_canon = True
                                #print(ms)
                                if (ms.up_id, ms.mapped_res, \
                                        str(ms.mapped_pos)) in stmts_by_site:
                                    has_annot_map_up_id = True
                    if pep_in_canon:
                        iso_id_but_pep_in_canon += 1

                    if has_annot_up_id:
                        has_stmts += 1
                        print("annot", gene_name, res, pos)
                    elif has_annot_map_up_id:
                        has_mapped_stmts += 1
                        print("mapped annot", gene_name, res, pos)

        no_can_id_ct = total_sites - has_canonical_iso_count
        text = ("No UP ID: %d / %d (%.1f)\n" %
                (no_up_id, total_sites, (100*no_up_id / total_sites)))
        text += ("No valid UP sequence: %d / %d (%.1f)\n" %
                 (no_valid_up, total_sites, (100*no_valid_up / total_sites)))
        text += ("Has canonical iso UP ID: %d / %d (%.1f)\n" %
                 (has_canonical_iso_count, total_sites,
                     (100*has_canonical_iso_count / total_sites)))
        non_canonical_iso_count = total_sites - has_canonical_iso_count
        text += ("Non-canonical iso UP ID: %d / %d (%.1f)\n" %
                 (non_canonical_iso_count, total_sites,
                     (100*non_canonical_iso_count / total_sites)))
        text += ("Annotated in INDRA DB: %d / %d (%.1f)\n" %
                 (has_stmts, total_sites, (100*has_stmts / total_sites)))
        text += ("Mapped site annotated in INDRA DB: %d / %d (%.1f)\n" %
                 (has_mapped_stmts, total_sites,
                    (100*has_mapped_stmts / total_sites)))
        text += ("No can. ID but motif in can. ID: %d / %d (%.1f)\n" %
                 (iso_id_but_pep_in_canon, no_can_id_ct,
                  (100*iso_id_but_pep_in_canon / no_can_id_ct)))

        print(text)
        with open(output_filename, 'wt') as f:
            f.write(text)
        # There are multiple uniprot IDs associated with each refseq ID.
        # To figure out whether it has any annotations, we turn each BRCA
        # site into a tuple (up_id, res, pos), and look this up in a dict
        # of sites with their annotations. If there are matches based on
        # ANY of the uniprot IDs we count that as a match.
    else:
        print("Argument must be one of map_uniprot, site_stats")
        sys.exit(1)
    """
    Stats
    - Number of proteins in Refseq with no UP ID mapping
    - Number of cases where site is not valid despite mapping (!)
    """
