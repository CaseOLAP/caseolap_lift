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
    
def mesh2triples(mesh_tree_to_id_file, mesh_tree_file, cvd_to_mesh_term_file):
    # reformat data into the following
    data = {h:[] for h in ['head','relation','tail','edge_type','weight']}
    
    # load the files
    mesh_tree_to_id_df = pd.read_csv(mesh_tree_to_id_file)
    mesh_tree_df = pd.read_csv(mesh_tree_file)    
    
    # make the mesh tree hierarchy
    for from_node, to_node in zip(mesh_tree_df['Disease (MeSH Tree)'], mesh_tree_df['Disease (MeSH Tree).1']):
        from_node_ = "MeSH_Tree_Disease:"+str(from_node)
        to_node_ = "MeSH_Tree_Disease:"+str(to_node)
        data['head'] += [from_node_]
        data['relation'] += [len(data['relation'])]
        data['tail'] += [to_node_]
        data['edge_type'] += ['MeSH_hierarchy']
        data['weight'] += [1]
    
    # the mesh tree node to mesh term mapping
    for from_node, to_node in zip(mesh_tree_to_id_df['Disease (MeSH Tree)'], mesh_tree_to_id_df['Disease (MeSH)']):
        data['head'] += [from_node]
        data['relation'] += [len(data['relation'])]
        data['tail'] += [to_node]
        data['edge_type'] += ['MeSH_is']
        data['weight'] += [1]
        
    # mapping 8CVDs to MeSH terms
    cvd_to_mesh_lines = [l.strip("\n") for l in open(cvd_to_mesh_term_file,"r").readlines()]
    cvd_to_mesh = {c:m.split(" ") for c,m in zip(['CM','ARR','CHD','VD','IHD','CCD','VOO','OTH'],cvd_to_mesh_lines)} #TODO 8CVDs need to be generalized
    for cvd, mesh_terms in cvd_to_mesh.items():
        for m in mesh_terms:
            mesh = "MeSH_Tree_Disease:"+str(m)
            data['head'] += [cvd]
            data['relation'] += [len(data['relation'])]
            data['tail'] += [mesh]
            data['edge_type'] += ['MeSH_CVD']
            data['weight'] += [1]

    mesh_kg = pd.DataFrame(data)
    mesh_kg["relation"] = mesh_kg["edge_type"]
    mesh_kg = mesh_kg.drop(columns = ["edge_type"])
    return mesh_kg
