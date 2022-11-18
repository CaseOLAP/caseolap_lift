import pandas as pd
import os
from grape import Graph
from kg_analysis.kg_analysis import get_edge_type_to_node_types_mapping, independent_edge_evaluation,add_predictions_to_kg

examine_core_proteins_only = True

core_proteins_file = '../output/core_proteins.txt'
node_path = '../output/graph_data/merged_node_list.tsv'
edge_path = '../output/graph_data/merged_edge_list.tsv'
output_folder = '../output/kg_analysis'
eval_output = os.path.join(output_folder,'eval_results.csv')
predictions_output = os.path.join(output_folder,'predictions.csv')

if not os.path.exists(output_folder):
   # Create a new directory because it does not exist
   os.makedirs(output_folder)

core_proteins = [l.strip("\n") for l in open(core_proteins_file,'r').readlines()]

g = Graph.from_csv(
  directed=False,
  node_path=node_path,
  edge_path=edge_path,
  verbose=True,
  nodes_column='node',
  node_list_node_types_column='node_type',
  default_node_type='None',
  sources_column='head',
  destinations_column='tail',
  edge_list_edge_types_column='relation',
  name="Mito KG"
)
g = g.remove_disconnected_nodes()

# how many of each edge type?
print(g.get_edge_type_names_counts_hashmap())

edge_type_to_node_types_mapping = get_edge_type_to_node_types_mapping(g, directed=False)

edge_pair = edge_type_to_node_types_mapping['CaseOLAP_score'][0]
eval_df, filtered_pred_df = independent_edge_evaluation(g, edge_pair)

# filter results
if examine_core_proteins_only:
    filtered_pred_df = filtered_pred_df[filtered_pred_df['destinations'].isin(core_proteins)]

# output
print(eval_df)
eval_df.to_csv(eval_output)
filtered_pred_df.to_csv(predictions_output,)
print("Number of predicted edges: %f"%filtered_pred_df.shape[0])
print("Number of unique proteins w/ predictions: %f"%len(set(filtered_pred_df['destinations'])))

## we were missing some protein-pathway relationships. Adding them now #TODO move this to kg creation
input_edges = pd.read_csv(edge_path, sep="\t")
added_pw_relations = pd.read_csv('../scratch/added_pathway_relations.csv')
added_pw_relations['weight'] = 1.0
merged_input_edges = input_edges.merge(added_pw_relations,on=['head','tail'],how='outer',indicator=True)
missing_pw_df = merged_input_edges[merged_input_edges['_merge']== 'right_only']
rels_to_add = missing_pw_df[['head','relation_y','tail','weight_y']]
rels_to_add.columns = ['head','relation','tail','weight']
input_edges = pd.concat([input_edges, rels_to_add])

filtered_pred_df = pd.read_csv(predictions_output)
print(filtered_pred_df.shape)
add_predictions_to_kg(input_edges,filtered_pred_df, output_folder=output_folder,debug=True)
