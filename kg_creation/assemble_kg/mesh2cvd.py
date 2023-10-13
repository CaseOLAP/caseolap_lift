#data handling
import pandas as pd
import json

from tqdm import tqdm

"""
NOTE: all functions are in order of calling
"""

def parse_mesh_tree(mesh_tree_file):
    lines = [l.strip("\n") for l in open(mesh_tree_file,"r").readlines()]
    mesh_term_to_code = {a:b for a,b in [l.split(";") for l in lines]}
    return mesh_term_to_code
    
def mesh2triples(mesh_tree_to_id_file, mesh_tree_file, disease_to_mesh, include_MeSH: bool):
    # reformat data into the following
    data = {h:[] for h in ['head','relation','tail','edge_type','weight']}

    # load the files
    mesh_tree_to_id_df = pd.read_csv(mesh_tree_to_id_file)
    mesh_tree_df = pd.read_csv(mesh_tree_file)

    #head and tail for the first df
    head_mesh_tree_df = [str(i) for i in mesh_tree_df['Disease (MeSH Tree)'].to_list()]
    tail_mesh_tree_df = [str(i) for i in mesh_tree_df['Disease (MeSH Tree).1'].to_list()]
    relation_mesh_tree_df = ['MeSH_hierarchy' for i in tail_mesh_tree_df]

    #head and tail for the second df
    head_mesh_tree_to_id_df = [str(i) for i in mesh_tree_to_id_df['Disease (MeSH Tree)'].to_list()]
    tail_mesh_tree_to_id_df = [str(i) for i in mesh_tree_to_id_df['Disease (MeSH)'].to_list()]
    relation_mesh_tree_to_id_df = ['MeSH_is' for i in tail_mesh_tree_to_id_df]

    head_cvd_mesh = []
    tail_cvd_mesh = []
    #iterate over dict
    for cvd, mesh_terms in disease_to_mesh.items():
        for m in mesh_terms:
            mesh = "MeSH_Tree_Disease:" + str(m)
            head_cvd_mesh.append(cvd)
            tail_cvd_mesh.append(mesh)
    relation_cvd_mesh = ["MeSH_CVD" for i in head_cvd_mesh]

    #combine final head column
    head = head_mesh_tree_df + head_mesh_tree_to_id_df + head_cvd_mesh
    tail = tail_mesh_tree_df + tail_mesh_tree_to_id_df + tail_cvd_mesh
    relation = relation_mesh_tree_df + relation_mesh_tree_to_id_df + relation_cvd_mesh
    weight = [1 for i in relation]

    mesh_kg = pd.DataFrame({"head" : head, "relation" : relation, "tail" : tail, "weight" : weight})

    #choose whether to include mesh tree or not
    if include_MeSH == False:
        mesh_kg = mesh_kg[mesh_kg["relation"] != "MeSH_hierarchy"]
    
    return mesh_kg
