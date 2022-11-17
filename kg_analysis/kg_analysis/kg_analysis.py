# imports
from grape import Graph
import pandas as pd

from grape.edge_prediction import PerceptronEdgePrediction
from grape.embedders import FirstOrderLINEEnsmallen

import os
import yaml
import pandas as pd
from grape import Graph

# node embedding imports
from embiggen.embedders.ensmallen_embedders.degree_spine import DegreeSPINE
from embiggen.embedders import GLEEEnsmallen
from embiggen.embedders import HOPEEnsmallen
from embiggen.embedders.pykeen_embedders.distmult import DistMultPyKEEN
from embiggen.embedders.pykeen_embedders.hole import HolEPyKEEN

# edge prediction imports
from embiggen.edge_prediction.edge_prediction_model import AbstractEdgePredictionModel
from grape.edge_prediction import PerceptronEdgePrediction
from grape.edge_prediction import MLPEdgePrediction
from grape.edge_prediction import GNNEdgePrediction


def add_predictions_to_kg(input_df, pred_df, output_folder='./results/kg_analysis/',debug=False):
    # format pred_df same as input_df
    pred_df.columns = ['idx','weight', 'head', 'tail', 'x']

    caseolap_edges = input_df[input_df['relation'] == 'CaseOLAP_score']
    caseolap_pairs = set(zip(caseolap_edges['head'], caseolap_edges['tail']))
    mask = [False if pair in caseolap_pairs else True for pair in
            zip(pred_df['head'], pred_df['tail'])]
    pred_df['relation'] = 'predicted_association'
    non_overlapping_filtered_pred_df = pred_df[mask]
    # non_overlapping_filtered_pred_df['relation'] = 'predicted_association'

    # append predictions to kg
    data_to_append = non_overlapping_filtered_pred_df[['head', 'relation', 'tail', 'weight']]
    out_df = pd.concat([input_df, data_to_append])

    output_file = os.path.join(output_folder,'kg_edges_with_predictions.csv')
    out_df.to_csv(output_file, index=False)

    if debug:
        print("%d predictions added"%(out_df.shape[0]-input_df.shape[0]))
        print("Predictions saved to %s"%(output_file))

def label_predictions_with_ground_truth(pred_df, test_graph, return_bool=True, directed=False):
    test_edges = set(test_graph.get_edge_node_names(directed=False))
    # consider reverse edges if directed is False
    if not directed:
        test_edges = test_edges.union([(t, f) for (f, t) in test_edges])

    # convert to dataframe to label edges
    test_edges_df = pd.DataFrame(test_edges)
    test_edges_df.columns = ['sources', 'destinations']

    labeled_pred_df = pred_df.merge(test_edges_df, how='left', indicator=True)
    test_truth = []
    for b in labeled_pred_df['_merge'] == 'both':
        if return_bool:
            y_ = b
        else:
            # return 1 or 0
            y_ = 1 if b else 0
        test_truth += [y_]
    labeled_pred_df['ground_truth'] = test_truth
    labeled_pred_df = labeled_pred_df.drop('_merge', axis=1)
    return labeled_pred_df


def label_negative_sample_pred(pred_df):
    pred_df['ground_truth'] = False
    return pred_df


def evaluate_predictions(m, labeled_pred_df):
    # convert DataFrame to numpy used for evaluation
    y_label = labeled_pred_df['ground_truth'].to_numpy()
    y_score = labeled_pred_df['predictions'].to_numpy()

    # handle case where no ground truth (for negative sampled edges)
    all_false = len(set(y_label)) == 1 and (~y_label[0])
    if all_false:
        e1 = {'auroc': float('NaN'), 'auprc': float('NaN')}
    else:
        e1 = m.evaluate_prediction_probabilities(y_label, y_score)
    e2 = m.evaluate_predictions(y_label, y_score)
    # return e1 | e2
    return  {**e1, **e2}


def get_edge_type_to_node_types_mapping(g, directed=False):
    '''
    This function returns a mapping from edge_type -> (from_types,to_types)
    '''
    # gather node and edge types for every edge as a DataFrame
    edge_to_type_dict = {h: [] for h in ["from", "to", "from_type", "to_type", "edge_type"]}
    for from_node_id, to_node_id in g.get_edge_node_ids(directed=False):
        # get node types
        from_node_type_ids = g.get_node_type_ids_from_node_id(from_node_id)
        from_node_type = [g.get_node_type_name_from_node_type_id(i) for i in from_node_type_ids]
        to_node_type_ids = g.get_node_type_ids_from_node_id(to_node_id)
        to_node_type = [g.get_node_type_name_from_node_type_id(i) for i in to_node_type_ids]

        # get edge type
        edge_id = g.get_edge_id_from_node_ids(from_node_id, to_node_id)
        edge_type = g.get_edge_type_name_from_edge_id(edge_id)

        # append to dict
        edge_to_type_dict['from'] += [from_node_id]
        edge_to_type_dict['to'] += [to_node_id]
        edge_to_type_dict['from_type'] += [from_node_type]
        edge_to_type_dict['to_type'] += [to_node_type]
        edge_to_type_dict['edge_type'] += [edge_type]
    edge_to_type_df = pd.DataFrame(edge_to_type_dict)

    # take unique node types for each edge type
    edge_type_to_node_types = {}
    for edge_type in set(edge_to_type_df['edge_type']):
        # only rows with the specified edge_type
        sub_df = edge_to_type_df[edge_to_type_df['edge_type'] == edge_type]

        # get unique pairs of node types
        unique_node_type_pairs = set()
        for from_node_types, to_node_types in zip(sub_df['from_type'], sub_df['to_type']):
            # enumerate all pairs, since these are lists of node types
            pairs = [(f, t) for f in from_node_types for t in to_node_types]

            # if undirected, do not include reverse node type
            # i.e. include (type_1,type_2) but not (type_2,type_1)
            if not directed:
                pairs_sorted = set()
                for f, t in pairs:
                    # sort pair alphabetically
                    pair = (f, t) if (f < t) else (t, f)
                    # keep unique pairs
                    pairs_sorted.add(pair)
                pairs = pairs_sorted

            # add unique pairs
            unique_node_type_pairs = unique_node_type_pairs.union(pairs)

        edge_type_to_node_types[edge_type] = list(unique_node_type_pairs)
    return edge_type_to_node_types


def filter_predictions(pred_df, threshold=0.9):
  filtered_pred_df = pred_df[pred_df['predictions'] > threshold]
  return filtered_pred_df


def independent_edge_evaluation(g, node_types):
    # get the node types we will be predicting edge between
    source_node_type_list = [node_types[0]]
    destination_node_type_list = [node_types[1]]
    print(source_node_type_list, destination_node_type_list)

    # split graph into train/test
    train, test = g.connected_holdout(train_size=0.7)

    # train embedding on train graph
    embedding = DistMultPyKEEN().fit_transform(train)


    # train model on train graph
    # model = PerceptronEdgePrediction(
    #     edge_features=None,
    #     number_of_edges_per_mini_batch=32,
    #     edge_embeddings="CosineSimilarity"
    # )

    #model = GNNEdgePrediction()
    model = MLPEdgePrediction()
    model.fit(
        graph=train,
        node_features=embedding
    )

    # predictions for train, test, and negative sampled graph
    train_pred_df = model.predict_proba_bipartite_graph_from_edge_node_types(graph=train,
                 node_features=embedding,
                 source_node_types=source_node_type_list,
                 destination_node_types=destination_node_type_list,
                 return_predictions_dataframe=True
                 )
    test_pred_df = model.predict_proba_bipartite_graph_from_edge_node_types(graph=test,
                node_features=embedding,
                source_node_types=source_node_type_list,
                destination_node_types=destination_node_type_list,
                return_predictions_dataframe=True
                )

    #TODO make negative_sampled edges be all edges NOT in train/test
    negative_sampled_pred_df = model.predict_proba_bipartite_graph_from_edge_node_types(
                    graph=g.sample_negative_graph(number_of_negative_samples=test.get_number_of_edges()),
                    node_features=embedding,
                    source_node_types=source_node_type_list,
                    destination_node_types=destination_node_type_list,
                    return_predictions_dataframe=True
                )

    #TODO do AUPRC curve and determine a good threshold
    thresh = 0.9

    # label predictions
    labeled_train_pred_df = label_predictions_with_ground_truth(train_pred_df, train)
    labeled_test_pred_df = label_predictions_with_ground_truth(test_pred_df, test)
    labeled_negative_sample_pred_df = label_negative_sample_pred(negative_sampled_pred_df)

    # evaluation
    data = [evaluate_predictions(model, labeled_train_pred_df),
            evaluate_predictions(model, labeled_test_pred_df),
            evaluate_predictions(model, labeled_negative_sample_pred_df)]
    eval_df = pd.DataFrame(data)

    filtered_pred_df = filter_predictions(labeled_negative_sample_pred_df, threshold=thresh)

    labeled_negative_sample_pred_df.to_csv("predictions.csv")

    return eval_df, filtered_pred_df

