#stl
import os
import shutil
import argparse
import time
import warnings
warnings.filterwarnings("ignore")

#data handling
import pandas as pd
import numpy as np
from tqdm import tqdm 

#grape
from grape import Graph

#assemble mesh tree/mesh--cvd
from assemble_kg.mesh2cvd import parse_mesh_tree
from assemble_kg.mesh2cvd import mesh2triples

#assemble caseolap (protein--cvd)
from assemble_kg.caseolap_cvd2protein import caseolap2triples

#assemble protein--pathway
from assemble_kg.protein2pathway import pathway2triples

#assemble reactome hierarchy
from assemble_kg.reactome_hierarchy import reactome2reactome

#assemble PPI protein--protein
from assemble_kg.protein2protein import protein2triples

#assemble TF protein--protein
from assemble_kg.transcription_protein2protein import transcription_protein2triples

#final graph data
from assemble_kg.final_graph_creation import graph_create

class caseolapLIFT_knowledge_graph:
    def __init__(self, include_STRING: bool, \
                 include_REACTOME: bool, \
                 include_MeSH: bool, \
                 include_TF: bool, \
                 caseolap_SCALING = "none"):
        """
        ARGS:
            include_STRING: bool whether to include STRING Protein Protein interactions
            include_REACTOME: bool whether to include REACTOME database
            include_REACTOME_HIERARCHY: bool whether to include all human pathways
            include_MeSH: bool whether to include MeSH hierarchy
            include_TF: bool whether to include transcription factor dependecies
            caseolap_SCALING: scaling type in caseolap: z-score, z-score+1, none

        RETURNS:
            NONE, but node list and edge list are saved in ../output/graph_data

        """

        #REQUIRED FILES
        ROOT = ".."

        #mesh (downloaded from irsyad-setup.py)
        mesh_tree_file = ROOT + "/data/MeSH/mtrees2021.bin"
        # categories_file = ROOT + "/caseolap/input/categories.txt"
        categories_file = ROOT + "/scratch/categories.txt"

        mesh_tree_to_mesh_id_file = ROOT + "/parsed_mappings/MeSH/edges_meshtree-IS-meshid_disease.csv"
        edges_meshtree2meshtree_hierarchy = ROOT + "/parsed_mappings/MeSH/edges_meshtree_to_meshtree.csv"

        #caseolap
        # caseolap_csv = ROOT + "/caseolap/result/all_proteins/all_caseolap.csv"
        caseolap_csv = ROOT + '/scratch/all_caseolap.csv'

        #reactome
        protein2pathway = ROOT + "/output/kg/reactome_edges.csv"

        #reactome hierarchy
        reactome_hierarchy = ROOT + "/data/Reactome/ReactomePathwaysRelation.txt"

        #string ppi
        string_edge = ROOT + "/output/kg/string_ppi.csv"

        #transcription factor dependence
        tf_edge = ROOT + "/output/kg/target_protein_id_2_tf_protein_id.json"

        print("\ncaseolapLIFT KG creation")

        print("\n-----------------------------------------------")

        print("\nfinding required files for assembly...")

        if os.path.exists(mesh_tree_file):
            print("\tmtrees2021.bin found")
        else:
            print("\tERROR: mtrees2021.bin not found")
            exit()

        if os.path.exists(mesh_tree_to_mesh_id_file):
            print("\tedges_meshtree-IS-meshid.csv found")
        else:
            print("\tERROR: edges_meshtree-IS-meshid.csv not found")
            exit()


        if os.path.exists(edges_meshtree2meshtree_hierarchy):
            print("\tedges_meshtree2meshtree_hierarchy.csv found")
        else:
            print("\tERROR: edges_meshtree2meshtree_hierarchy.csv not found")
            exit()

        if os.path.exists(categories_file):
            print("\tcategories.txt found")
        else:
            print("\tERROR: categories.txt not found")
            exit()

        if os.path.exists(caseolap_csv):
            print("\tcaseolap.csv found")
        else:
            print("\tERROR: caseolap.csv not found")
            exit()

        if os.path.exists(protein2pathway):
            print("\treactome_edges.csv found")
        else:
            print("\tERROR: reactome_edges.csv not found")
            exit()

        if os.path.exists(reactome_hierarchy):
            print("\tReactomePathwaysRelation.csv found")
        else:
            print("\tERROR: ReactomePathwaysRelation.csv not found")
            exit()
        
        if os.path.exists(string_edge):
            print("\tstring_ppi.csv found")
        else:
            print("\tERROR: string_ppi.csv not found")
            exit()

        if os.path.exists(tf_edge):
            print("\ttarget_protein_id_2_tf_protein_id.json found")
        else:
            print("\tERROR: target_protein_id_2_tf_protein_id.json not found")
            exit()


        print("\nassembling mesh--cvd relation...", end = " ")
        mesh_term_to_code = parse_mesh_tree(mesh_tree_file)
        

        #includes mesh tree, also conencted with the 8 cvd categories
        mesh_kg = mesh2triples(mesh_tree_to_mesh_id_file, edges_meshtree2meshtree_hierarchy,categories_file, include_MeSH)
        if include_MeSH:
            print("success (MeSH hierarchy included)")
        else:
            print("success")

        #caseolap
        print("assembling caseolap--protein relation...", end = " ")
        if caseolap_SCALING == "z-score":
            caseolap_kg = caseolap2triples(caseolap_csv, SCALING = "z-score")
        elif caseolap_SCALING == "z-score+1":
            caseolap_kg = caseolap2triples(caseolap_csv, SCALING = "z-score+1") 
        elif caseolap_SCALING == "none":
            caseolap_kg = caseolap2triples(caseolap_csv) 
        else:
            caseolap_kg = caseolap2triples(caseolap_csv) 
            caseolap_SCALING = "none"
        print("success (scaling = %s)" % caseolap_SCALING)


        reactome_kg = None
        reactome_hierarchy_kg = None
        if include_REACTOME:
            print("assembling protein--pathway relation...", end = " ")
            reactome_kg = pathway2triples(protein2pathway) #line that does all the work
            reactome_hierarchy_kg = reactome2reactome(reactome_hierarchy)
            print("success")

        ppi_kg = None
        if include_STRING:
            print("assembling STRING protein--protein relation...", end = " ")
            ppi_kg = protein2triples(string_edge)
            print("success")

        tf_kg = None
        if include_TF:
            print("assembling TRANSCRIPTION FACTOR protein--protein relation...", end = " ")
            tf_kg = transcription_protein2triples(tf_edge)
            print("success")

        if not os.path.exists("../output/graph_data"):
            os.mkdir("../output/graph_data")
        print("\nedge list and node list graph data in: ../output/graph_data/")

        
        graph_create(mesh_kg, caseolap_kg, reactome_kg, ppi_kg, tf_kg, reactome_hierarchy_kg)
        print("done")

        print("\n-----------------------------------------------\n")
        print("ready to initialize GRAPE knowledge graph\n")

    @staticmethod
    def knowledge_graph():
        return Graph.from_csv(node_path ="../output/graph_data/merged_node_list.tsv",
                              node_list_separator = "\t",
                              node_list_header = True,
                              nodes_column = "node",
                              node_list_node_types_column = "node_type",
                              edge_path ="../output/graph_data/merged_edge_list.tsv",
                              edge_list_separator = "\t",
                              edge_list_header = True,
                              sources_column = "head",
                              destinations_column = "tail",
                              edge_list_numeric_node_ids = False,
                              weights_column = "weight",
                              edge_list_edge_types_column = "relation",
                              directed = False,
                              verbose = True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="KG creation using CaseolapLIFT")
    parser.add_argument("--include_STRING", default=True, action="store_true", help = "include STRING protein to protein interactions")
    parser.add_argument("--include_REACTOME", default=True, action="store_true", help = "include REACTOME pathways to connect with proteins")
    parser.add_argument("--include_MESH", default=True, action="store_true", help = "include MeSH tree hierarchy")
    parser.add_argument("--caseolap_SCALING", default="none", type=str, help="caseolap scaling type")
    parser.add_argument("--include_TF", default=True, action="store_true", help="include transcription factor dependency")

    args = parser.parse_args()

    caseolapLIFT = caseolapLIFT_knowledge_graph(args.include_STRING, args.include_REACTOME, args.include_MESH, args.include_TF, args.caseolap_SCALING)
    graph = caseolapLIFT.knowledge_graph()
