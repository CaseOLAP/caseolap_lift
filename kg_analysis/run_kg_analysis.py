
from grape import Graph
from kg_analysis.kg_analysis import get_edge_type_to_node_types_mapping, independent_edge_evaluation

examine_core_proteins_only = True

core_proteins_file = '../output/core_proteins.txt'
core_proteins = [l.strip("\n") for l in open(core_proteins_file,'r').readlines()]


g = Graph.from_csv(
  directed=False,
  node_path='./input/merged_node_list.tsv',
  edge_path='./input/merged_edge_list.tsv',
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

print(eval_df)
eval_df.to_csv("1_eval_results.csv")


if examine_core_proteins_only:
  filtered_pred_df = filtered_pred_df[filtered_pred_df['destinations'].isin(core_proteins)]
filtered_pred_df.to_csv("1_predictions.csv")
