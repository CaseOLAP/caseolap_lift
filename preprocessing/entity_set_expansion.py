import sys

# setting path
sys.path.append('..')

# importing
from utils.biomed_apis import *
from utils.other_functions import *
import pandas as pd
import numpy as np
import os


def load_mappings(resource_to_mapping_files, mappings_folder="./parsed_mappings/"):

    resource_to_mappings = {resource: [] for resource in resource_to_mapping_files.keys()}

    for resource, file_list in resource_to_mapping_files.items():
        resource_folder = os.path.join(mappings_folder,resource)

        for mapping_file in file_list:
            # check the required files exist
            file_path = os.path.join(resource_folder,mapping_file)
            print(file_path)
            if not os.path.exists(file_path):
                print("Required %s file from %s resource does not exist! Exit" % (file_path,resource))
                sys.exit(1)

            # load file and store
            mapping = json.load(open(file_path,'r'))
            resource_to_mappings[resource] += [mapping]

    return resource_to_mappings


def get_interacting_partners(proteins, k=1, score_thresh=0.975,
                             output_folder="../output/kg", debug=False):
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
        results = get_id_mapping_results_search_uniprot_api(link)

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
        results = get_id_mapping_results_search_uniprot_api(link)

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


def get_pathway_partners(proteins, pathway2protein, count_thresh=sys.maxsize, proportion_thresh=0.25, debug=False,
                         output_folder="../output/kg/"):
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


def prepare_resource_mappings(include_ppi, include_pathways, include_transcription_factor_dependence):
    # Right now, PPI does not use any pre-existing mappings
    if include_ppi:
        pass

    resources_to_include = ['GO'] # requires GO by default
    # only add mappings if they are flagged as True
    if include_pathways:
        resources_to_include += ['Reactome']
    if include_transcription_factor_dependence:
        resources_to_include += ['Transcription_Factor_Dependence']

    # all mappings
    resource_to_mapping_files = {'GO':['go_id_to_links.json','go2protein.json'],
                            'Reactome':['pathway2protein.json'],
                            'Transcription_Factor_Dependence':['all_entrez2uniprot.json','all_uniprot2entrez.json','id2synonyms_not_case_varied.json','gene_name_2_protein_id.json']#TODO
                            }
    return {resource: mapping_files for resource,mapping_files in resource_to_mapping_files.items() if resource in resources_to_include}


def prepare_subcellular_compartment_proteins(parameters,
                                             output_folder="../output",
                                             mapping_folder='../parsed_mappings',
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
    resource_mapping_files = prepare_resource_mappings(include_ppi, include_pathways, include_transcription_factor_dependence)
    print(mapping_folder)
    resource_mappings = load_mappings(resource_mapping_files, mappings_folder=mapping_folder)
    go_id_to_links, go2protein = resource_mappings['GO']

    # Get organelle-specific proteins
    organelle_proteins = get_proteins_from_go(go_term, go_id_to_links, go2protein)

    print("%d proteins relevant to go term %s" % (len(organelle_proteins), go_term))
    proteins_of_interest = set(organelle_proteins)
    if include_ppi:
        ppi_proteins = get_interacting_partners(organelle_proteins,
                                                k=ppi_k,
                                                score_thresh=ppi_score_thresh,
                                                output_folder=output_folder,
                                                debug=debug)
        print(type(proteins_of_interest))
        print(type(ppi_proteins))
        proteins_of_interest = proteins_of_interest.union(set(ppi_proteins))
        print("%d proteins added from protein-protein interaction" % len(ppi_proteins))
    if include_pathways:
        pathway2protein = resource_mappings['Reactome'][0]
        pathway_proteins = get_pathway_partners(organelle_proteins, pathway2protein,
                                                count_thresh=pw_count_thresh,
                                                proportion_thresh=pw_proportion_thresh,
                                                debug=debug)
        proteins_of_interest = proteins_of_interest.union(set(pathway_proteins))
        print("%d proteins added with common pathways" % len(pathway_proteins))
    if include_transcription_factor_dependence:
        protein_ids_2_gene_ids,gene_ids_2_protein_ids,gene_name_2_protein_id,tf_gene_name_2_target_gene_name = resource_mappings['Transcription_Factor_Dependence']
        tfd_proteins = get_transcription_factor_dependence_partners(organelle_proteins,protein_ids_2_gene_ids,gene_ids_2_protein_ids,gene_name_2_protein_id,tf_gene_name_2_target_gene_name)
        proteins_of_interest = proteins_of_interest.union(set(tfd_proteins))
        print("%d proteins from transcription factor dependence" % len(tfd_proteins))

    print("In total, %d proteins of interest assembled" % (len(proteins_of_interest)))

    if output_folder:
        out_protein_list_file = os.path.join(output_folder,"proteins_of_interest.txt")
        with open(out_protein_list_file,"w") as out_file:
            out_file.write("\n".join(proteins_of_interest))
            print("Written to file %s"%out_protein_list_file)
        out_organelle_list_file = os.path.join(output_folder,"core_proteins.txt")
        with open(out_organelle_list_file,"w") as out_file:
            out_file.write("\n".join(organelle_proteins))
            print("Written to file %s"%out_organelle_list_file)
    return proteins_of_interest

def get_transcription_factor_dependence_partners(proteins,protein_ids_2_gene_ids,gene_ids_2_protein_ids,gene_name_2_protein_id,tf_gene_name_2_target_gene_name, output_folder = "../output/kg"):

    '''Dictionary'''
    gene_name_2_gene_id = dict()

    # Gene Name
    for gene_name, protein_ids in gene_name_2_protein_id.items():

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
    print(len(tf_gene_name_2_target_gene_name))
    print("Gene name to target gene name")
    print([(k,v) for k,v in tf_gene_name_2_target_gene_name.items()][:10])
    print("Gene name to protein id")
    print([(k,v) for k,v in gene_name_2_protein_id.items()][:10])

    for tf_gene_name, target_gene_names in tf_gene_name_2_target_gene_name.items():
    
        # TF Gene Name -is- TF Protein ID
        try:
            tf_protein_ids = gene_name_2_protein_id[tf_gene_name]
            tf_protein_ids = list(tf_protein_ids)
            print(len(tf_protein_ids))
        except:
            #print("%s gene not found!"%(tf_gene_name))
            continue
        # Target Gene Names -is- Target Gene ID
        #print([(k,v) for k,v in gene_name_2_gene_id.items()][:10])
        #print(target_gene_names[:10])
        for target_gene_name in target_gene_names:
            try:
                #print(target_gene_name in gene_name_2_gene_id)
                target_gene_id = gene_name_2_gene_id[target_gene_name]
 #               print(target_gene_id)
                target_gene_id = list(target_gene_id)[0]
                for tf_protein_id in tf_protein_ids:
                    tf_protein_id_2_target_gene_id.setdefault(tf_protein_id, set()).add(target_gene_id)
  #              print(len(tf_protein_ids))
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
    json.dump(tf_protein_id_2_target_gene_id, open(os.path.join(output_folder,"tf_protein_id_2_target_gene_id.json"), 'w'))
    json.dump(target_protein_id_2_tf_protein_id, open(os.path.join(output_folder,"target_protein_id_2_tf_protein_id.json"), 'w'))
    json.dump(tf_protein_id_2_target_protein_id, open(os.path.join(output_folder,"tf_protein_id_2_target_protein_id.json"), 'w'))
    # get tf dependence proteins
    proteins_of_interest = set()
    
    # query from_list:
    for p, prots in target_protein_id_2_tf_protein_id.items():
        if p in proteins:
            proteins_of_interest = proteins_of_interest.union(prots)
    # query to_list
    for p, prots in tf_protein_id_2_target_protein_id.items():
        if p in proteins:
            proteins_of_interest = proteins_of_interest.union(prots)
    
    added_proteins = proteins_of_interest.difference(proteins)
    print("Added proteins: %d"%(len(added_proteins)))
    return added_proteins

parameters = {'go-term': 'GO:0005739',
              'include_ppi': False,
              'ppi_k': 1, 'ppi_score_thresh': 0.99,
              'include_pathways': False,
              'pw_count_thresh': 4,
              'pw_proportion_thresh': 0.50,
              'include_transcription_factor_dependence': True}

root_directory = '/caseolap_lift_shared_folder'
root_directory = '../'
mapping_folder = os.path.join(root_directory,'parsed_mappings')
output_folder = os.path.join(root_directory,'output')
if not os.path.exists(output_folder):
   os.makedirs(output_folder)
   kg_output_folder = os.path.join(output_folder,"kg")
   print(kg_output_folder)
   os.makedirs(kg_output_folder)
proteins = prepare_subcellular_compartment_proteins(parameters, mapping_folder=mapping_folder, output_folder=output_folder, debug=False)
