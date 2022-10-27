import sys

# setting path
sys.path.append('..')

# importing
from utils.biomed_apis import *
from utils.other_functions import *

import pandas as pd
import numpy as np
import os


def extract_go_hierarchy(go_dict):
    go_id_to_links = {}

    for go_term, go_fields in go_dict.items():

        parents = []
        if 'is_a' in go_fields:
            parents += [p.split(" ! ")[0] for p in go_fields['is_a']]
        if 'relationship' in go_fields:
            parents += [p.split(" ! ")[0].strip('part_of ') for p in go_fields['relationship'] if 'part_of' in p]

        # add to dict
        for g in [go_term] + parents:
            if g not in go_id_to_links:
                go_id_to_links[g] = dict()
                go_id_to_links[g]['parents'] = []
                go_id_to_links[g]['children'] = []
        for p in parents:
            go_id_to_links[p]['children'] += [go_term]
            go_id_to_links[go_term]['parents'] += [p]
    return go_id_to_links


def load_mappings(output_folder):
    # TODO have these passed in
    global go_term_cell_comp2protein
    global go_id_to_links
    global go2protein
    data_folder = "../data/"

    '''Protein to GO'''
    # !rm 'data/goa_human.gaf'
    # !wget -N -P data/ http://geneontology.org/gene-associations/goa_human.gaf.gz
    # !gunzip 'data/goa_human.gaf.gz'
    protein2go, go2protein = map_protein2go_ids(output_folder=data_folder)

    '''GO to GO'''
    # !wget -N -P data/ http://purl.obolibrary.org/obo/go/go-basic.obo

    # GO Term information
    go_dict = get_go_term_information(output_folder=data_folder)

    # GO Term ID -> GO Term Name
    go_id2term = map_go_id2term(go_dict)

    # Separate GO terms into 3 lists for each type
    go_bio_proc, go_cell_comp, go_mol_func = separate_go_terms_into_three_lists_for_each_tree(go_dict)

    # Cellular Compartment Names
    cell_comp_name = list()
    cell_comp_id2name = dict()

    for ID in go_cell_comp:
        name = go_id2term[ID]
        cell_comp_name.append(name)
        cell_comp_id2name[ID] = name

    # GO Term Cellular Compartments -> Protein
    go_term_cell_comp2protein, protein2go_term_cell_comp = get_protein2go_term(protein2go, cell_comp_id2name)

    # extract go hierarchy (is_a and part_of)
    go_id_to_links = extract_go_hierarchy(go_dict)





def get_interacting_partners(proteins, k=1, score_thresh=0.975,
                             output_folder="./output/kg", debug=False):
    '''
    Identifies the interacting partners of the input protein list
    '''

    # convert UniProt IDs to STRING IDs
    proteins_as_string, str_to_up, up_to_str = convert_uniprot_to_string_ids(proteins,
                                                                             debug=debug, return_mappings=True)

    # Get interacting partners
    interacting_partners_df, counts = k_hop_interactors_string_api(proteins_as_string,
                                                                   k=k,
                                                                   score_thresh=score_thresh,
                                                                   debug=debug,
                                                                   return_counts=True)
    interacting_partners = set(interacting_partners_df['query_ensp']).union(
        set(interacting_partners_df['partner_ensp']))

    # convert back to UniProt IDs
    added_string_ids = interacting_partners.difference(proteins_as_string)
    added_uniprot_ids, up_to_str_1, str_to_up_1 = convert_string_to_uniprot_ids(added_string_ids, debug=debug,
                                                                                return_mappings=True)

    # combine mappings
    string_id_to_uniprot_id = {**str_to_up, **str_to_up_1}
    uniprot_id_to_string_id = {**up_to_str, **up_to_str_1}

    # assemble kg triples
    head_list = [string_id_to_uniprot_id[s] for s in interacting_partners_df['query_ensp'] if
                 s in string_id_to_uniprot_id]
    tail_list = []
    for s in interacting_partners_df['partner_ensp']:
        if s in string_id_to_uniprot_id:
            tail_list += [string_id_to_uniprot_id[s]]
        else:
            tail_list += [""]
    interacting_partners_df['h'] = head_list
    interacting_partners_df['r'] = 'binds_to'
    interacting_partners_df['t'] = tail_list
    # drop missing tail rows (string could not map to uniprot)
    string_triples = interacting_partners_df[['h', 'r', 't']].replace('', np.nan).dropna(subset=['t'])

    # output edges to kg
    if output_folder:
        output_file = os.path.join(output_folder, "string_ppi.csv")
        string_triples.to_csv(output_file, index=False)

    # TODO remake added_uniprot_ids to only include those in the kg edges
    added_uniprot_ids_1 = set(string_triples['t'])

    if debug:
        print(added_uniprot_ids)
        print("%d additional UniProt proteins from protein-protein interactions with k=%d and score_thresh=%f"
              % (len(added_uniprot_ids_1), k, score_thresh))

    return added_uniprot_ids

#TODO move to prepare_kg_data.py or utils
def convert_uniprot_to_string_ids(proteins, debug=False, return_mappings=False):
    protein_set = set(proteins)

    # Convert UniProt IDs to STRING IDs
    job_id = submit_id_mapping_uniprot_api(
        from_db='UniProtKB_AC-ID',
        to_db='STRING',
        ids=protein_set)

    # Wait until result finished.
    if check_id_mapping_results_ready_uniprot_api(job_id):
        link = get_id_mapping_results_link_uniprot_api(job_id)
        results = get_id_mapping_results_search_UniProtAPI(link)

    # Parse results
    stringISuniprot, uniprotISstring = dict(), dict()
    for u2s in results['results']:
        uniprot = u2s['from']
        string = u2s['to']

        if '9606' in string:  # TODO fix hard-coded. TODO this is human specific
            uniprotISstring[uniprot] = string
            stringISuniprot[string] = uniprot

    if debug:
        print(len(stringISuniprot), '/', len(protein_set), 'UniProt proteins aligned with STRING')

    added_proteins = list(stringISuniprot.keys())
    if return_mappings:
        return added_proteins, stringISuniprot, uniprotISstring
    return added_proteins

#TODO move to prepare_kg_data.py or utils
def convert_string_to_uniprot_ids(proteins, debug=False, return_mappings=False):
    protein_set = set(proteins)

    # Convert UniProt IDs to STRING IDs
    job_id = submit_id_mapping_uniprot_api(
        from_db='STRING',
        to_db='UniProtKB',
        ids=protein_set)

    # Wait until result finished.
    if check_id_mapping_results_ready_uniprot_api(job_id):
        link = get_id_mapping_results_link_uniprot_api(job_id)
        results = get_id_mapping_results_search_UniProtAPI(link)

    # Parse results
    stringISuniprot, uniprotISstring = dict(), dict()
    for u2s in results['results']:
        string = u2s['from']
        uniprot = u2s['to']['primaryAccession']

        if '9606' in string:  # TODO fix hard-coded. TODO this is human specific
            uniprotISstring[uniprot] = string
            stringISuniprot[string] = uniprot

    if debug:
        print(len(uniprotISstring), '/', len(protein_set), 'STRING proteins aligned with UniProt')

    added_proteins = list(uniprotISstring.keys())
    if return_mappings:
        return added_proteins, uniprotISstring, stringISuniprot
    return added_proteins


def get_pathway_partners(proteins, count_thresh=sys.maxsize, proportion_thresh=0.25, debug=False,
                         output_folder="../output/kg/"):
    protein2pathway, pathway2protein = load_protein2pathway_data('../data/UniProt2Reactome.txt',
                                                                 output_folder="../data/") #TODO move to prepare_kg_data.py or utils

    # Pick pathways that have a substantial amount of mitochondrial proteins
    # How do we define 'substantial'?
    pathways_of_interest = list()
    for pathway, pathway_proteins in pathway2protein.items():
        prots_in_pw = set(proteins).intersection(set(pathway_proteins))
        coverage = len(prots_in_pw) / len(pathway_proteins)
        if len(prots_in_pw) >= count_thresh or coverage >= proportion_thresh:
            pathways_of_interest.append(pathway)

    # Proteins in mitochondrial proteins' pathways
    proteins_in_pathways = set()
    for pw in pathways_of_interest:
        proteins = set(pathway2protein[pw])
        proteins_in_pathways = proteins_in_pathways.union(proteins)

    # output to kg
    if output_folder:
        output_file = os.path.join(output_folder, "reactome_edges.csv")
        filtered_pathway2protein = {k: v for k, v in pathway2protein.items() if k in pathways_of_interest}
        h = []
        r = []
        t = []
        for pw, proteins in filtered_pathway2protein.items():
            for p in proteins:
                h += [p]
                r += ['participates_in_pathway']
                t += [pw]
        out_df = pd.DataFrame({'h': h, 'r': r, 't': t})
        out_df.to_csv(output_file, index=False)

    added_proteins = proteins_in_pathways.difference(proteins)

    if debug:
        print("Settings: count_thresh=%d; proportion_thresh=%f" % (count_thresh, proportion_thresh))
        print("%d pathways identified" % len(pathways_of_interest))
        print("%d proteins added" % (len(added_proteins)))

    return added_proteins


#TODO make it accept multiple GO ID's
def get_proteins_from_go(go_id, go_id_to_links, go2protein):
    # go through the hierarchy and get proteins
    go_queue = [go_id]
    extracted_proteins = set()
    while len(go_queue) > 0:
        g = go_queue.pop()
        if g in go2protein:
            prots = go2protein[g]
            extracted_proteins = extracted_proteins.union(prots)
            go_queue += go_id_to_links[g]['children']

    print("%d proteins extracted" % len(extracted_proteins))
    return extracted_proteins

def prepare_subcellular_compartment_proteins(parameters,
                                             output_folder="./data",
                                             debug=False):
    # parse parameters
    go_term = parameters['go-term']
    include_ppi = parameters['include_ppi']
    ppi_k = parameters['ppi_k']
    ppi_score_thresh = parameters['ppi_score_thresh']
    include_pathways = parameters['include_pathways']
    pw_count_thresh = parameters['pw_count_thresh']
    pw_proportion_thresh = parameters['pw_proportion_thresh']
    include_transcription_factor_dependence = parameters['include_transcription_factor_dependence']

    # load data
    load_mappings(output_folder)

    # Get organelle-specific proteins
    #     organelle_proteins = get_organelle_proteins(go_term_cell_comp2protein, organelle_name = organelle_name)
    organelle_proteins = get_proteins_from_go(go_term, go_id_to_links, go2protein)

    print("%d proteins relevant to go term %s" % (len(organelle_proteins), go_term))
    proteins_of_interest = set(organelle_proteins)
    if include_ppi:
        ppi_proteins = get_interacting_partners(organelle_proteins,
                                                k=ppi_k,
                                                score_thresh=ppi_score_thresh,
                                                debug=debug)
        print(type(proteins_of_interest))
        print(type(ppi_proteins))
        proteins_of_interest = proteins_of_interest.union(set(ppi_proteins))
        print("%d proteins added from protein-protein interaction" % len(ppi_proteins))
    if include_pathways:
        pathway_proteins = get_pathway_partners(organelle_proteins,
                                                count_thresh=pw_count_thresh,
                                                proportion_thresh=pw_proportion_thresh,
                                                debug=debug)
        proteins_of_interest = proteins_of_interest.union(set(pathway_proteins))
        print("%d proteins added with common pathways" % len(pathway_proteins))
    if include_transcription_factor_dependence:
        tfd_proteins = get_transcription_factor_dependence_partners(organelle_proteins)
        proteins_of_interest = proteins_of_interest.union(set(tfd_proteins))
        print("%d proteins from transcription factor dependence" % len(tfd_proteins))

    print("In total, %d proteins of interest assembled" % (len(proteins_of_interest)))

    if output_folder:
        out_protein_list_file = os.path.join(output_folder,"proteins_of_interest.txt")
        with open(out_protein_list_file,"w") as out_file:
            out_file.write("\n".join(proteins_of_interest))
            print("Written to file %s"%out_protein_list_file)
    return proteins_of_interest

def get_transcription_factor_dependence_partners(proteins,data_folder = '../data/GRNdb'):
    #######
    # TFD #
    #######
    global tf_gene_name_2_target_gene_name

    tf_gene_name_2_target_gene_name = dict()  # TF gene name to target gene name
    gene_names = set()  # Gene Names (Tfs and targets)
    tf_gene_names = set()  # Transcription Factor gene names
    target_gene_names = set()  # Target gene names
    for file in os.listdir(data_folder):
        if 'txt' not in file:
            continue

        for line in open(os.path.join(data_folder, file)):
            line = line.strip().split('\t')
            if "There is something wrong" in line[0]:
                break

            confidence = line[5]
            if confidence == 'High':
                # Gene names
                tf_gene_name = line[0]
                targ_gene_name = line[1]

                # Save gene names
                gene_names.add(tf_gene_name)
                gene_names.add(targ_gene_name)
                tf_gene_names.add(tf_gene_name)
                target_gene_names.add(targ_gene_name)

                # TF Gene -targets-> Target Gene
                tf_gene_name_2_target_gene_name.setdefault(tf_gene_name, set()).add(targ_gene_name)

    # Change the values from set into a list
    tf_gene_name_2_target_gene_name = switch_dictset_to_dictlist(tf_gene_name_2_target_gene_name)

    print(len(tf_gene_name_2_target_gene_name))

    # TODO better API for mapping gene to protein

    '''Dictionary (all known mappings)'''
    protein_ids_2_gene_ids = json.load(
        open('../data/all_uniprot2entrez.json', 'r'))
    gene_ids_2_protein_ids = json.load(
        open('../data/all_entrez2uniprot.json', 'r'))

    '''Dictionary'''
    gene_name_2_gene_id = dict()

    # Gene Name
    for gene_name, protein_ids in gene_ids_2_protein_ids.items():

        # Protein IDs
        for protein_id in protein_ids:

            # ProteinID -is- Gene IDs
            try:
                gene_ids = protein_ids_2_gene_ids[protein_id]
                for gene_id in gene_ids:
                    # Gene Name -is- Gene ID
                    gene_name_2_gene_id.setdefault(gene_name, set()).add(gene_id)
            except:
                continue
    '''Check mappings'''
    multiple_gene_ids_per_gene_name, one_gene_id_per_gene_name = set(), set()

    for k, v in gene_name_2_gene_id.items():
        if len(v) > 1:
            multiple_gene_ids_per_gene_name.add(k)
        else:
            one_gene_id_per_gene_name.add(k)

    print(len(multiple_gene_ids_per_gene_name), 'gene names with multiple gene IDs (possibly bad)')
    print(len(one_gene_id_per_gene_name), 'gene names with one gene ID (good)')

    '''Remove unclear mappings'''
    for gene_name, gene_ids in gene_ids_2_protein_ids.copy().items():
        if len(gene_ids) > 1:
            gene_ids_2_protein_ids.pop(gene_name)

    '''Checking that no gene names correspond to multiple protein IDs'''
    something_went_wrong = False
    for gene_name, protein_ids in gene_ids_2_protein_ids.copy().items():
        if len(protein_ids) > 1:
            print(gene_name, protein_ids)
            something_went_wrong = True
    if not something_went_wrong:
        print('All good')

    tf_protein_id_2_target_gene_id = dict()
    tf_protein_id_2_target_protein_id = dict()
    target_protein_id_2_tf_protein_id = dict()

    for tf_gene_name, target_gene_names in tf_gene_name_2_target_gene_name.items():

        # TF Gene Name -is- TF Protein ID
        try:
            tf_protein_ids = gene_ids_2_protein_ids[tf_gene_name]
            tf_protein_ids = list(tf_protein_ids)
        except:
            continue

        # Target Gene Names -is- Target Gene ID
        for target_gene_name in target_gene_names:
            try:
                target_gene_id = gene_name_2_gene_id[target_gene_name]
                target_gene_id = list(target_gene_id)[0]
                for tf_protein_id in tf_protein_ids:
                    tf_protein_id_2_target_gene_id.setdefault(tf_protein_id, set()).add(target_gene_id)
            except:
                continue

                # Protein ID's Gene -is targeted by-> TF Gene
            try:
                target_gene_protein_ids = gene_ids_2_protein_ids[target_gene_id]
                for target_gene_protein_id in target_gene_protein_ids:
                    for tf_protein_id in tf_protein_ids:
                        tf_protein_id_2_target_protein_id.setdefault(tf_protein_id, set()).add(target_gene_protein_id)
                        target_protein_id_2_tf_protein_id.setdefault(target_gene_protein_id, set()).add(tf_protein_id)
            except:
                continue

    ''' Output the Protein-Gene relationships'''
    tf_protein_id_2_target_gene_id = switch_dictset_to_dictlist(tf_protein_id_2_target_gene_id)
    target_protein_id_2_tf_protein_id = switch_dictset_to_dictlist(target_protein_id_2_tf_protein_id)
    tf_protein_id_2_target_protein_id = switch_dictset_to_dictlist(tf_protein_id_2_target_protein_id)

    # TODO move to prepare_kg_data.py or utils
    json.dump(tf_protein_id_2_target_gene_id, open('../data/tf_protein_id_2_target_gene_id.json', 'w'))
    json.dump(target_protein_id_2_tf_protein_id, open('../data/target_protein_id_2_tf_protein_id.json', 'w'))
    json.dump(tf_protein_id_2_target_protein_id, open('../data/tf_protein_id_2_target_protein_id.json', 'w'))

    organelle_tf2target, organelle_target2tf = dict(), dict()

    # TODO where do we get this variable from?
    # protein_id2names = json.load(open('data/id2syns_not_case_varied.json'))
    print(len(proteins))
    print(tf_protein_id_2_target_protein_id)
    for protein in proteins:
        # organelle Proteins' Targets
        try:
            target_proteins = tf_protein_id_2_target_protein_id[protein]
            for target_protein in target_proteins:
                organelle_tf2target.setdefault(protein, set()).add(target_protein)

        except:
            pass

        # organellechondrial Proteins' Transcription Factors
        try:
            their_tfs = target_protein_id_2_tf_protein_id[protein]
            for tf in their_tfs:
                organelle_target2tf.setdefault(protein, set()).add(tf)
        except:
            pass
    print(len(organelle_tf2target))
    print(len(organelle_target2tf))
    all_tfs = list()
    for organelle_prot, tfs in organelle_target2tf.items():
        all_tfs += tfs
    all_tfs = list(set(all_tfs))

    non_organelle_prot_are_tfs = set()
    for protein in all_tfs:
        if protein not in proteins:
            non_organelle_prot_are_tfs.add(protein)
    print(len(all_tfs), 'TFs target organelle\'s proteins\' genes', '(' + \
          str(len(non_organelle_prot_are_tfs)), 'non-organelle\'s protein TFs)')

    all_targets = list()
    for organelle_prot, targets in organelle_tf2target.items():
        all_targets += targets
    all_targets = list(set(all_targets))

    non_organelle_prot_are_targs = set()
    for protein in all_targets:
        if protein not in proteins:
            non_organelle_prot_are_targs.add(protein)

    print(len(organelle_tf2target.keys()), 'organelle\'s proteins are TFs')
    print(len(all_targets), 'proteins\' genes are targeted by organelle\'s proteins', '(' + \
          str(len(non_organelle_prot_are_targs)), 'non-organelle\'s protein targets)')

    # added_proteins = proteins_of_interest.difference(my_protein_list)
    # #     print("Added proteins: %d" % (len(added_proteins)))
    # #     return added_proteins
    return set()

parameters = {'go-term': 'GO:0005739',
              'include_ppi': False,
              'ppi_k': 1, 'ppi_score_thresh': 0.99,
              'include_pathways': True,
              'pw_count_thresh': 4,
              'pw_proportion_thresh': 0.50,
              'include_transcription_factor_dependence': True}

# parameters = {'go-term':'GO:0005739',
#              'include_ppi':True,
#              'ppi_k':1, 'ppi_score_thresh': 0.99,
#              'include_pathways':True,
#              'pw_count_thresh':sys.maxsize,
#              'pw_proportion_thresh':0.50,
#              'include_transcription_factor_dependence':False}
# TODO make below into a function, accepting parameters object
print(os.getcwd())
output_folder = "../output"
if not os.path.exists(output_folder):
   os.makedirs(output_folder)
   kg_output_folder = os.path.join(output_folder,"kg")
   print(kg_output_folder)
   os.makedirs(kg_output_folder)
proteins = prepare_subcellular_compartment_proteins(parameters, output_folder=output_folder, debug=False)
