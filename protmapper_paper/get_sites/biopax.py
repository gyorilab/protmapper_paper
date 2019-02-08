import sys
import pickle
from indra.tools import assemble_corpus as ac
from indra.statements import Agent, ModCondition, Phosphorylation, \
                             Dephosphorylation
from indra.sources import biopax
from indra.sources.biopax import processor as bpc
from indra.sources.biopax import pathway_commons_client as pcc

def save_modified_agents(owl_file, output_file):
    print('Reading %s...' % owl_file)
    model = pcc.owl_to_model(owl_file)
    mf_class = bpc._bpimpl('ModificationFeature')

    objs = model.getObjects().toArray()

    agents = []
    for obj in objs:
        if not isinstance(obj, mf_class):
            continue
        try:
            mc = bpc.BiopaxProcessor._extract_mod_from_feature(obj)
        except Exception as e:
            print('ERROR: ' + str(e))
            continue
        if not mc or not mc.residue or not mc.position:
            continue

        proteins = obj.getFeatureOf().toArray()
        if not proteins:
            continue
        for protein in proteins:
            name = bpc.BiopaxProcessor._get_element_name(protein)
            db_refs = bpc.BiopaxProcessor._get_db_refs(protein)
            agent = Agent(name, mods=[mc], db_refs=db_refs)
            reactions = protein.getParticipantOf().toArray()
            if not reactions:
                upstream = protein.getMemberPhysicalEntityOf().toArray()
                for u in upstream:
                    reactions += u.getParticipantOf().toArray()
            for reaction in reactions:
                controls = reaction.getControlledOf().toArray()
                if not controls:
                    agents.append(agent)
                for contr in controls:
                    agents.append(agent)

    with open(output_file, 'wb') as fh:
        pickle.dump(agents, fh)


def save_phosphorylation_stmts(owl_file, pkl_file):
    bp = biopax.process_owl(owl_file)
    phos = ac.filter_by_type(bp.statements, Phosphorylation)
    dephos = ac.filter_by_type(bp.statements, Dephosphorylation)
    stmts = phos + dephos
    respos_stmts = [s for s in stmts if s.residue and s.position]
    ac.dump_statements(respos_stmts, pkl_file)


if __name__ == '__main__':
    owl_file = sys.argv[1]
    pkl_file = sys.argv[2]
    #save_phosphorylation_stmts(owl_file, pkl_file)
    #save_modified_agents(owl_file, pkl_file)


