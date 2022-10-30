import os
import json
from .biomedkg_utils import switch_dictset_to_dictlist, switch_dictlist_to_dictset


def load_protein2pathway_data_PE(uniprot2reactome_file_path, output_folder="./data"):
        
    protein2pathway, pathway2protein = dict(), dict()

    for line in open(uniprot2reactome_file_path):
        line = line.strip().split('\t')

        # Protein -involved in-> Pathway
        try:  
            species, protein, pathway = line[7], line[0], line[3]

            if species.lower() == 'homo sapiens':
                protein2pathway.setdefault(protein, list()).append(pathway)
                pathway2protein.setdefault(pathway, list()).append(protein)
        except: pass
        
     
    # Export
    protein2pathway = switch_dictlist_to_dictset(protein2pathway)
    protein2pathway = switch_dictset_to_dictlist(protein2pathway)
    
    pathway2protein = switch_dictlist_to_dictset(pathway2protein)
    pathway2protein = switch_dictset_to_dictlist(pathway2protein)

    json.dump(protein2pathway, open(os.path.join(output_folder,'protein2pathway'),'w'))
    json.dump(pathway2protein, open(os.path.join(output_folder,'pathway2protein'),'w'))
    
    return protein2pathway, pathway2protein


def map_top_proteins2go(protein_DIS, protein2go, go_tree_of_interest, go_id2term):
    protein2go_DIS, go2protein_DIS, go_ids_DIS, go_terms_DIS = dict(), dict(), dict(), dict()

    # Each protein
    for protein, values in protein_DIS.items():
        if len(values) == 0: continue
            
        # Each GO term for this type of tree
        try: go_ids = protein2go[protein]
        except: continue
        for go_id in go_ids:
            if go_id in go_tree_of_interest:
                
                # Protein -> GO, GO Term ID, GO Term Name
                protein2go_DIS.setdefault(protein,list()).append(go_id)
                go2protein_DIS.setdefault(go_id,list()).append(protein)
                go_ids_DIS[go_id] = go_ids_DIS.get(go_id, 0) + 1
                go_term = go_id2term[go_id]
                go_terms_DIS[go_term] = go_terms_DIS.get(go_term, 0) + 1
                
          
    go_ids_DIS = dict(sorted(go_ids_DIS.items(), key = lambda x:x[1], reverse=True))
    go_terms_DIS = dict(sorted(go_terms_DIS.items(), key = lambda x:x[1], reverse = True))
    
    return protein2go_DIS, go2protein_DIS, go_ids_DIS, go_terms_DIS

def get_go_term_information(output_folder="./data/"):
    # Convert GO obo file to dict
    ID = ''
    go_dict = dict()
    for line in open(os.path.join(output_folder,'go-basic.obo')):
        if line.startswith('id: '):
            ID = line.split('id: ')[1].strip('\n')
            continue
        if ': ' in line and ID != '':
            go_term = line.split(': ')[0]
            go_info = line.split(': ')[1].strip('\n')
            go_dict.setdefault(ID,dict()).setdefault(go_term,[]).append(go_info)
            
    return go_dict
      
    
def map_go_id2term(go_dict):
    # GO Term ID -> GO Term Name
    go_id2term = dict()
    for go_id in go_dict:
        go_id2term[go_id] = go_dict[go_id]['name'][0]
        
    return go_id2term


def separate_go_terms_into_three_lists_for_each_tree(go_dict):
    go_bio_proc, go_cell_comp, go_mol_func = list(), list(), list()

    for go, go_info in go_dict.items():
        tree_name = go_info['namespace'][0]
        if tree_name == 'biological_process':
            go_bio_proc.append(go)
        elif tree_name == 'cellular_component':
            go_cell_comp.append(go)
        elif tree_name == 'molecular_function':
            go_mol_func.append(go)    
            
    return go_bio_proc, go_cell_comp, go_mol_func



def get_protein2go_term(protein2go, go_tree_id2name):
    '''
    FUNCTION:
    - Map proteins to the GO Term, not just GO ID
    '''
    protein2go_term = dict()
    go_term2protein = dict()

    for protein, gos in protein2go.items():

        # Each GO Term
        for go in gos:

            try:
                # Protein -> GO Term
                go_term = go_tree_id2name[go]
                protein2go_term.setdefault(protein, list()).append(go_term)
                go_term2protein.setdefault(go_term, list()).append(protein)
            except:
                # Not in relevant GO tree/ontology
                continue

    go_term2protein = switch_dictset_to_dictlist(switch_dictlist_to_dictset(go_term2protein))
    protein2go_term = switch_dictset_to_dictlist(switch_dictlist_to_dictset(protein2go_term))

    return go_term2protein, protein2go_term



def get_organelle_proteins(go_term2protein, organelle_name):
    '''
    FUNCTION:
    - Gets the proteins in the cellular compartment you specify
    
    PARAMS:
    - organelle_name (str): The cellular compartment whose proteins
      you want. For example, 'mitochondri' finds all cellular compartments
      containing 'mitochondri' including 'mitochondrion' and 
      'mitochondrial membrane'
    - go_term2protein (dict): From above function. GO Term mapping to proteins with that GO attribute.
    '''
    organelle_proteins = list()

    for go_term, proteins in go_term2protein.items():
        if organelle_name in go_term:
            organelle_proteins += proteins
            
    organelle_proteins = list(set(organelle_proteins))

    return organelle_proteins    
