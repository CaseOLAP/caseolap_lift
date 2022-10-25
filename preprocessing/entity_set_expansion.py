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
    interacting_partners_df, counts = k_hop_interactors_STRINGAPI(proteins_as_string,
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


def convert_uniprot_to_string_ids(proteins, debug=False, return_mappings=False):
    protein_set = set(proteins)

    # Convert UniProt IDs to STRING IDs
    job_id = submit_id_mapping_UniProtAPI(
        from_db='UniProtKB_AC-ID',
        to_db='STRING',
        ids=protein_set)

    # Wait until result finished.
    if check_id_mapping_results_ready_UniProtAPI(job_id):
        link = get_id_mapping_results_link_UniProtAPI(job_id)
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


def convert_string_to_uniprot_ids(proteins, debug=False, return_mappings=False):
    protein_set = set(proteins)

    # Convert UniProt IDs to STRING IDs
    job_id = submit_id_mapping_UniProtAPI(
        from_db='STRING',
        to_db='UniProtKB',
        ids=protein_set)

    # Wait until result finished.
    if check_id_mapping_results_ready_UniProtAPI(job_id):
        link = get_id_mapping_results_link_UniProtAPI(job_id)
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
                                                                 output_folder="../data/")

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


def get_transcription_factor_dependence_partners(my_protein_list):
    # load target_protein_id_2_tf_protein_id
    global target_protein_id_2_tf_protein_id
    global tf_protein_id_2_target_protein_id

    # get tf dependence proteins
    proteins_of_interest = set()

    # query from_list:
    for p, prots in target_protein_id_2_tf_protein_id.items():
        if p in my_protein_list:
            proteins_of_interest = proteins_of_interest.union(prots)
    # query to_list
    for p, prots in tf_protein_id_2_target_protein_id.items():
        if p in my_protein_list:
            proteins_of_interest = proteins_of_interest.union(prots)

    added_proteins = proteins_of_interest.difference(my_protein_list)
    print("Added proteins: %d" % (len(added_proteins)))
    return added_proteins


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


def output_kg_edges(core_proteins, added_proteins, output_folder='./output/kg'):
    '''
    core proteins are identified by go term
    added protein is a dict{origin->proteins}, for origins as 'ppi','pathway','tfd'
    '''
    # make folder if does not exist


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
    return proteins_of_interest


parameters = {'go-term': 'GO:0005739',
              'include_ppi': False,
              'ppi_k': 1, 'ppi_score_thresh': 0.99,
              'include_pathways': True,
              'pw_count_thresh': 4,
              'pw_proportion_thresh': 0.50,
              'include_transcription_factor_dependence': False}

# parameters = {'go-term':'GO:0005739',
#              'include_ppi':True,
#              'ppi_k':1, 'ppi_score_thresh': 0.99,
#              'include_pathways':True,
#              'pw_count_thresh':sys.maxsize,
#              'pw_proportion_thresh':0.50,
#              'include_transcription_factor_dependence':False}
print(os.getcwd())
output_folder = "../output"
if not os.path.exists(output_folder):
   os.makedirs(output_folder)
   kg_output_folder = os.path.join(output_folder,"kg")
   print(kg_output_folder)
   os.makedirs(kg_output_folder)
proteins = prepare_subcellular_compartment_proteins(parameters, output_folder=output_folder, debug=False)