import pandas as pd
import numpy as np

def mesh_nodes(mesh_kg: pd.DataFrame) -> pd.DataFrame:
    """
    node list for mesh
    [head, relation, tail]
    """
    all_nodes = list(set(mesh_kg["head"].to_list() + mesh_kg["tail"].to_list()))
    cvd_types = ["IHD", "CM", "ARR", "VD", "CHD", "CCD", "VOO", "OTH"]
    node_types = ["CVD" if i in cvd_types else "MeSH_Tree_Disease" for i in all_nodes]
    node_list = pd.DataFrame({"node": all_nodes, "node_type": node_types})
    return node_list

def caseolap_nodes(caseolap_kg: pd.DataFrame) -> pd.DataFrame:
    all_nodes = list(set(caseolap_kg["head"].to_list() + caseolap_kg["tail"].to_list()))
    cvd_types = ["IHD", "CM", "ARR", "VD", "CHD", "CCD", "VOO", "OTH"]
    node_types = ["CVD" if i in cvd_types else "Protein" for i in all_nodes]
    node_list = pd.DataFrame({"node": all_nodes, "node_type": node_types})
    return node_list

    
def reactome_nodes(reactome_kg: pd.DataFrame) -> pd.DataFrame:
    protein_nodes = list(set(reactome_kg["head"]))
    protein_type = ["Protein" for i in protein_nodes]

    pathway_nodes = list(set(reactome_kg["tail"]))
    pathway_type = ["Reactome_Pathway" for i in pathway_nodes]

    all_nodes = protein_nodes + pathway_nodes
    node_types = protein_type + pathway_type
    node_list = pd.DataFrame({"node": all_nodes, "node_type": node_types})
    return node_list


def ppi_nodes(ppi_kg: pd.DataFrame) -> pd.DataFrame:
    all_nodes = list(set(ppi_kg["head"].to_list() + ppi_kg["tail"].to_list()))
    node_types = ["Protein" for i in all_nodes]
    node_list = pd.DataFrame({"node": all_nodes, "node_type": node_types})
    return node_list


def graph_create(mesh_kg: pd.DataFrame, caseolap_kg: pd.DataFrame, \
                 reactome_kg: pd.DataFrame, ppi_kg: pd.DataFrame) -> pd.DataFrame:
    """
    given the edge list of all compiled datasets, generate merged edge list and merged node list
    """
    #filter out the stuff
    df_list = [mesh_kg, caseolap_kg, reactome_kg, ppi_kg]
    df_list = [i for i in df_list if i is not None]

    merged_edges = pd.concat(df_list)
    merged_edges.to_csv("graph_data/merged_edge_list.tsv", sep = "\t", index = False)

    node_list = pd.concat([mesh_nodes(mesh_kg), caseolap_nodes(caseolap_kg), \
                          reactome_nodes(reactome_kg), ppi_nodes(ppi_kg)])
    node_list = node_list.drop_duplicates(subset = ["node"])
    node_list.to_csv("graph_data/merged_node_list.tsv", sep = "\t", index = False)