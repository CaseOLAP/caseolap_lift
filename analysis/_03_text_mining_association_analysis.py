# general
import os
import pandas as pd
import numpy as np
from itertools import combinations
import networkx as nx
import json
import csv
import requests
import urllib

# stats
from sklearn import datasets # is this necessary?
from scipy.stats import zscore
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram

# plotting
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from matplotlib.lines import Line2D
from matplotlib.pyplot import figure
from matplotlib.tri import Triangulation
import seaborn as sns
from upsetplot import plot as upplot
 
'''
merge unirefs =============================
'''
def parse_uniref_to_uniprot_list(input_file, debug=False):
    # map uniref IDs back to proteins
    map_df = pd.read_csv(input_file, sep = "\t")
    map_df.columns = ['Mapped Proteins','Protein Groups']
    
    # make a reverse mapping from UniRef90 to UniProt list
    uniref_to_uniprot_list = {}
    for upid, uniref_id in zip(map_df['Mapped Proteins'], map_df['Protein Groups']):
        if uniref_id not in uniref_to_uniprot_list:
            uniref_to_uniprot_list[uniref_id] = []
        uniref_to_uniprot_list[uniref_id] += [upid]
    
    # return reverse mapping as a dataframe as well to add to summary table
    ret_df = pd.DataFrame([uniref_to_uniprot_list.keys(), ["; ".join(v) for v in uniref_to_uniprot_list.values()]]).T
    ret_df.columns = [ "Protein Groups", "Mapped Proteins"]
    ret_df = ret_df.set_index("Protein Groups")
    
    if debug:
        all_uniprot_ids = set()
        for ulist in uniref_to_uniprot_list.values():
            all_uniprot_ids = all_uniprot_ids.union(set(ulist))
        print("%d uniref ids and %d uniprot ids parsed"%(len(uniref_to_uniprot_list),
                                                        len(all_uniprot_ids)))
        
    return uniref_to_uniprot_list, ret_df


def prepare_ranked_ent(file_name):
    with open(file_name, 'r') as j:
        contents = json.loads(j.read())
        df = pd.DataFrame(contents)
        return df
    
    
def prepare_synonyms(ranked_ent_df):
    '''
    This function prepares the synonym column to be added to the summary table
    UniRef ID -> list of synonyms
    '''
    uniref_to_synonyms = {}
    for uniref in ranked_ent_df.index:
        all_synonyms_dict = ranked_ent_df.loc[uniref]['TOTAL']
        uniref_to_synonyms[uniref] = list(all_synonyms_dict.keys())
    return uniref_to_synonyms

def merge_redundant_ids(prot_to_syn, scores):
    # protein groups to synonym
    temp_sets = {k: set(v) for k, v in prot_to_syn.items()}  # duplicate synonyms will be in same order
    temp = {k: ", ".join(v) for k, v in temp_sets.items()}  # joining the synonyms
    prot_to_syn_df = pd.DataFrame.from_dict(temp, orient='index')  # creating dataframe out of dict
    prot_to_syn_df.columns = ["Synonyms"]
    prot_to_syn_df.index.name = "entity"

    # creating the data frame with "Mapped Proteins" as index
    scores = scores.set_index("entity")
    merged_df = scores.merge(prot_to_syn_df, right_index=True, left_index=True, how='left')

    # obtaining list of all synonyms
    all_syns = prot_to_syn_df[['Synonyms']].values

    # no_repeat_syns = set(all_syns)
    no_repeat_syns = []
    for i in all_syns:
        if i not in no_repeat_syns:
            no_repeat_syns.append(i)

    # for each synonym
    for syn in no_repeat_syns:

        for i in syn:
            syn = i

        # take a subset of that data frame - only those synonyms
        syn_df = merged_df.loc[merged_df['Synonyms'] == syn]
        # no merging if it's the only entry with those synonyms
        if syn_df.shape[0] == 1:
            continue

        # obtaining list of all mapped proteins in this subset
        all_mappedProteins = syn_df.index.to_list()

        # obtaining list of all mapped proteins in reverse
        reverse_all_mappedProteins = []
        for i in all_mappedProteins:
            reverse_all_mappedProteins.insert(0, i)

        # these are the proteins that will be marked 'deleted' later on
        deleted_prots = []

        # starting at the beginning of the protein list
        for protein in all_mappedProteins:

            # delete protein from the reversed list so we don't match it against itself
            del reverse_all_mappedProteins[-1]

            # if this protein was not marked 'deleted', obtain its scores and synonyms
            if protein not in deleted_prots:

                # protein of focus' score
                prot_score = scores.loc[protein].to_numpy()

                # protein of focus' synonym
                prot_synonym = prot_to_syn[protein]

                # these will stores the UniRefs and mapped proteins to be appended
                mappedProts_to_append = []

                # now to compared against another protein
                for other_protein in reverse_all_mappedProteins:

                    # if this protein was not marked 'deleted', obtain its scores and synonyms as well
                    if other_protein != 'deleted':

                        # score of protein being compared
                        other_prot_score = scores.loc[other_protein].to_numpy()

                        # synonym of protein being compared
                        other_prot_syns = prot_to_syn[other_protein]

                        # if protein of focus and the one being compared match in score and synonym:
                        if (prot_synonym == other_prot_syns) and (
                        np.array_equal(prot_score, other_prot_score, equal_nan=True)):

                            # other_protein and its UniRef need to be appended
                            mappedProts_to_append.append(other_protein)

                            # marked other_protein from original and reversed list as 'deleted'
                            deleted_prots.append(other_protein)
                            reverse_all_mappedProteins = [i.replace(other_protein, 'deleted') for i in
                                                          reverse_all_mappedProteins]

                            # removing other_protein entry in df
                            merged_df = merged_df.drop(index=other_protein)

                # append other proteins to same entry as protein of focus
                mappedProts_to_append.append(protein)

                mappedProts_to_append = ";".join(mappedProts_to_append)
                protein = str(protein)

                merged_df = merged_df.rename(index={protein: mappedProts_to_append})

    # setting the index back to protein groups
    merged_df = merged_df.reset_index()

    syns_column = merged_df.pop('Synonyms')
    merged_df.insert(1, 'Synonyms', syns_column)

    return merged_df


'''
idk =============================
'''

# heatmap visualizing the caseolap scores

def make_caseolap_score_heatmap(data, score_cutoff=0.05,
                    out_file = "./output/CaseOLAP_score_heatmap.pdf"):
    '''
    Makes a heatmap with caseolap scores. Only shows rows with at least one 
    score above score_cutoff
    '''    

    # remove protein column
    data = data.copy()
    data = data.drop(['entity'],axis=1)
    # only take those with at least one score above score_cutoff
    data = data[( data > score_cutoff).sum(axis=1) > 1]
    max_score = data.to_numpy().max()
    print(data.head())
    print(data.shape)
    print(max_score)
    # plot the heamap
    chart = sns.clustermap(data, cmap='YlGnBu', xticklabels=True, 
                                            yticklabels=False,col_cluster=False,row_cluster=True,
                                            vmin = 0.0, vmax=max_score, method='average',
                                            cbar_kws={"shrink": .5})
    #chart_axis = chart.ax_row_dendrogram.axes
    chart.ax_row_dendrogram.set_visible(False)

    # make it save/display in high-resolution
    fig = chart.fig
    fig.set_size_inches(10,10)

    # colorbar
    # fig.colorbar(ax.get_children()[0], orientation="horizontal")

    # # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,box_inches = 'tight')
    
    return chart

# pca for the caseolap scores
def make_caseolap_pca(data, out_file='./output/CaseOLAP_pca.pdf'):

    # remove protein column
    data = data.copy()
    data = data.drop(['entity'],axis=1)
    # scale data
    scaler = StandardScaler()
    scaler.fit(data)
    X=scaler.transform(data)

    pca = PCA()
    pca.fit(X)
    
    print("explained variance")
    print(pca.explained_variance_ratio_)
    print("singular values")
    print(pca.singular_values_)
    
    X_transformed = pca.transform(X)     

    legend_elements = [Line2D([0], [0], color='r', label='IHD'),
                                     Line2D([0], [0], color='g', label='CM'),
                                     Line2D([0], [0], color='b', label='ARR'),
                                     Line2D([0], [0], color='c', label='VD'),
                                     Line2D([0], [0], color='burlywood', label='CHD'),
                                     Line2D([0], [0], color='purple', label='CCD'),
                                     Line2D([0], [0], color='orange', label='VOO'),
                                     Line2D([0], [0], color='chartreuse', label='OTH'),
                                     ]
    fig, ax = plt.subplots()

    ax.legend(handles=legend_elements, loc='upper right')
    #Call the function. Use only the 2 PCs.
    score = X_transformed[:,0:2]
    coeff = np.transpose(pca.components_[0:2, :])
    
    xs = score[:,0]
    ys = score[:,1]
    n = coeff.shape[0]
    scalex = 1.0/(xs.max() - xs.min())
    scaley = 1.0/(ys.max() - ys.min())

    plt.scatter(xs*scalex ,ys*scaley, marker = '.', alpha = 0.5)
    plt.grid()
    
    #plot all the arrows along with catagory names
    IHD_arrow = plt.arrow(0, 0, coeff[0,0], coeff[0,1], color = 'r', alpha = 0.7,linestyle = '-',linewidth = 2, overhang=0.2, label = 'IHD')
    plt.text(coeff[0,0]* 1.15, coeff[0,1] * 1.15, "IHD", color = 'k', ha = 'center', va = 'bottom')

    plt.arrow(0, 0, coeff[1,0], coeff[1,1], color = 'g', alpha = 0.7,linestyle = '-',linewidth = 2, overhang=0.2)
    plt.text(coeff[1,0]* 1.15, coeff[1,1] * 1.15, "CM", color = 'k', ha = 'right', va = 'center')

    plt.arrow(0, 0, coeff[2,0], coeff[2,1], color = 'b', alpha = 0.7,linestyle = '-',linewidth = 2, overhang=0.2)
    plt.text(coeff[2,0]* 1.15, coeff[2,1] * 1.15, "ARR", color = 'k', ha = 'center', va = 'center')

    plt.arrow(0, 0, coeff[3,0], coeff[3,1], color = 'c', alpha = 0.7,linestyle = '-',linewidth = 2, overhang=0.2)
    plt.text(coeff[3,0]* 1.15, coeff[3,1] * 1.15, "VD", color = 'k', ha = 'right', va = 'bottom')

    plt.arrow(0, 0, coeff[4,0], coeff[4,1], color = 'burlywood', alpha = 0.7,linestyle = '-',linewidth = 2, overhang=0.2)
    plt.text(coeff[4,0]* 1.15, coeff[4,1] * 1.15, "CHD", color = 'k', ha = 'center', va = 'center')

    plt.arrow(0, 0, coeff[5,0], coeff[5,1], color = 'purple', alpha = 0.7,linestyle = '-',linewidth = 2, overhang=0.2)
    plt.text(coeff[5,0]* 1.15, coeff[5,1] * 1.15, "CCD", color = 'k', ha = 'center', va = 'center')

    plt.arrow(0, 0, coeff[6,0], coeff[6,1], color = 'orange', alpha = 0.7,linestyle = '-',linewidth = 2, overhang=0.2)
    plt.text(coeff[6,0]* 1.15, coeff[6,1] * 1.15, "VOO", color = 'k', ha = 'right', va = 'top')

    plt.arrow(0, 0, coeff[7,0], coeff[7,1], color = 'chartreuse', alpha = 0.7,linestyle = '-',linewidth = 2, overhang=0.2)
    plt.text(coeff[7,0]* 1.15, coeff[7,1] * 1.15, "OTH", color = 'k', ha = 'center', va = 'bottom')

    plt.xlim(-0.15,1.25)
    plt.ylim(-0.8,0.8)
    plt.xlabel("PC{}".format(1))
    plt.ylabel("PC{}".format(2))
    
    plt.grid()
    plt.axhline(y=0, color="black", linestyle="-")
    plt.axvline(x=0, color="black")
    plt.show()


    # # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,box_inches = 'tight')
            
            
'''

'''

def convert_to_zscore(original_df, include_zeros = False, columns_to_ignore=['Protein Groups']):

    # make a copy of original, drop unnecessary columns
    new_df = original_df.copy().drop(columns_to_ignore, axis=1)
    temp_headers = new_df.columns

    # mask zeros if needed
    if not include_zeros:
        new_df[new_df == 0] = np.nan

    # convert to z-score
    zscore_values = zscore(new_df, axis = 0, ddof = 6, nan_policy = 'omit')
    
    # convert back to dataframe
    new_df = pd.DataFrame(zscore_values)
    new_df.columns = temp_headers
    # add dropped columns and return in same column ordering
    new_df = new_df.join(original_df[columns_to_ignore])
    new_df = new_df[original_df.columns]

    return new_df

def cvd_proteomes_violin_plot(data, subproteome_plot_values, show_figure=True,
                              out_file = './caseolap_zscore_violin_plot.pdf'):

    sns.set_theme(style="darkgrid")
    figure(figsize=(12, 6), dpi=80)

    # plot our data
    ax = sns.violinplot(data=data)
    sns.swarmplot( data=subproteome_plot_values, color="gray", edgecolor="white", size=8)
    plt.plot([7.5, -0.5], [3, 3], linewidth=2)
    
                
    # gather counts for figure
    cvd_to_sig_count = {}
    cvd_to_x = {}
    for cvd in subproteome_plot_values.columns:
        sig_ps = list(subproteome_plot_values[~subproteome_plot_values[cvd].isna()][cvd])
        cvd_to_sig_count[cvd] = len(sig_ps)
        cvd_to_x[cvd] = len(cvd_to_x)

    # plot count above violin plot
    if cvd_to_x and cvd_to_sig_count:
        max_y = max(list(subproteome_plot_values.max()))
        text_height = max_y+0.75
        for cvd, x_ind in cvd_to_x.items():
            pcount = cvd_to_sig_count[cvd]
            # plot text
            plt.text(x_ind,text_height,str(pcount),
                    {'weight':'bold', 'fontsize':14, 'color':'royalblue'},
                    ha='center', va='center')


    # formatting
    ax.set_facecolor('white')
    # ax.set_title(plot_title, fontsize = 20)
    ax.tick_params(axis='both', labelsize=14)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_color('black')
    ax.spines['bottom'].set_color('black')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    ax.grid(axis='y', alpha=0.75)
    # ax.set_xlabel(x_axis_text, fontsize = 20)
    y_axis_text = "z-score"
    ax.set_ylabel(y_axis_text, fontsize = 20, weight='bold')
    ax.set_xticklabels(ax.get_xmajorticklabels(), weight='bold')
#     ax.set_yticklabels(ax.get_ymajorticklabels(), weight='bold')

    plt.savefig(out_file)
    if show_figure:
        plt.show()
    

'''

'''
def identify_cvd_unique_uniref_ids(zscores_df, score_threshold = 3.0, debug=False):
    
    # get list of cvds from header
    cvds = list(zscores_df.columns)

    # get unirefs which pass score threshold for that cvd
    cvd_to_unirefs = {}
    for cvd in cvds:
        # for a uniref id to be unique, it needs to pass score threshold ONLY for that CVD
        sub_df = zscores_df[ (zscores_df[cvd] >= score_threshold) ]
        cvd_to_unirefs[cvd] = set(sub_df.index)

    # find uniref ids that only show up in one cvd
    cvd_to_unique_unirefs = {}
    for this_cvd in cvds:
        unirefs = cvd_to_unirefs[this_cvd]
        for other_cvd in cvds:
            if this_cvd != other_cvd:
                unirefs = unirefs.difference(cvd_to_unirefs[other_cvd])
        cvd_to_unique_unirefs[this_cvd] = unirefs

    if debug:
        for cvd, unirefs in cvd_to_unirefs.items():
            print("%s has %d uniref ids with scores above %f"%(cvd,len(unirefs),score_threshold))
        count = 0
        for cvd, unirefs in cvd_to_unique_unirefs.items():
            print("%s has %d unique uniref ids"%(cvd,len(unirefs)))
            count += len(unirefs)
        print("In total, %d unique uniref ids"%count)
    
    return cvd_to_unique_unirefs

# Make table to show how many proteins in each group, with different threhsolds of z-scores
def zscore_cutoff_table(zscores_df):
    result_dict = {'z-score cutoff':[]}
    for thresh in [1.0,2.0,3.0,5.0,10.0]:
        # perform filtering
        cvd_to_unique_unirefs = identify_cvd_unique_uniref_ids(zscores_df,
                                                            score_threshold = thresh, 
                                                            debug = False)
        
        #initialize dict to hold values
        for cvd in cvd_to_unique_unirefs.keys():
            if cvd not in result_dict:
                result_dict[cvd] = []

        # add values to table
        result_dict['z-score cutoff'] += [thresh]
        for cvd in cvd_to_unique_unirefs.keys():
            # reports the cardinality of number of unirefs
            result_dict[cvd] += [len(cvd_to_unique_unirefs[cvd])]

    result_df = pd.DataFrame(result_dict)
    return result_df


def make_category_i_heatmap(category_i_scores, sort_by='max', n=5, 
                            out_file = './figures/Category_I_Heatmap.pdf'):
    # make a copy of original data
    scores = category_i_scores.copy(deep=True)
    
    # determine how to sort
    if sort_by == 'max':
        scores['max'] = scores.max(axis=1)
        scores = scores.sort_values(['max'], ascending = False).drop_duplicates()
        del scores['max']
    elif sort_by == 'mean':
        scores['mean'] = scores.mean(axis=1)
        scores = scores.sort_values(['mean'], ascending = False).drop_duplicates()
        del scores['mean']
    
    # take the first n rows
    scores = scores.head(n)
    
    # rename the y-labels with '*' indicating multiple UniRef90 IDs
    new_y_labels = []
    for label in scores.index:
        label_ = label
        if ";" in label:
            label_ = label.split(";")[0]+"*"
        new_y_labels += [label_]
    scores.index = new_y_labels
    
    # plot
    plt.figure(figsize = [8,4])
    s = sns.heatmap(scores,\
                cmap="YlGnBu",\
                yticklabels=True,\
                vmin = 4,vmax = 8, \
                annot=True)
    s.set_ylabel('Protein Groups')
    
    # make it save/display in high-resolution
    fig = plt.gcf()
    fig.set_size_inches(10,8)

    # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,box_inches = 'tight')
    

# heatmap for the proteins unique to each cvd
def make_category_ii_heatmap(summary_table, n=3,
                        out_file='./output/unique_unirefs_heatmap.pdf'):
    
    category_headers = [c for c in summary_table.columns if "Category II" in c]
    cvds = [c.split(" ")[0] for c in category_headers]
    
    # extract heatmap values 
    heatmap_rows = []
    for col in category_headers:
        cvd = col.split(' ')[0]
        
        # get the rows where it is true in Category II
        sub_df = summary_table[summary_table[col] == True]
        # keep the top n row names
        sub_df = sub_df.sort_values([cvd], ascending = False)
        heatmap_rows += list(sub_df.head(n).index)
#         print(list(sub_df.head(n).index))
#         print(sub_df.shape)
    heatmap_df = summary_table[summary_table.index.isin(heatmap_rows)][cvds]
    heatmap_df = heatmap_df.reindex(heatmap_rows)
    
    # rename the y-labels with '*' indicating multiple UniRef90 IDs
    new_y_labels = []
    for label in heatmap_df.index:
        label_ = label
        if ";" in label:
            label_ = label.split(";")[0]+"*"
        new_y_labels += [label_]
    heatmap_df.index = new_y_labels
    heatmap_df.index.name='entity'
    
    max_score = np.nan_to_num(heatmap_df.to_numpy()).max()
    
    # plot the heamap
    chart = sns.clustermap(heatmap_df, cmap='YlGnBu', xticklabels=True, 
                                            yticklabels=True,col_cluster=False,row_cluster=False,
                                            vmin = -1.0, vmax=max_score, annot=True, fmt='.2f',
                                            cbar_kws={"shrink": .5}, figsize=(10,20))

    chart.ax_heatmap.xaxis.tick_top() # x axis on top
    chart.ax_heatmap.xaxis.set_label_position('top')
    chart.ax_heatmap.yaxis.tick_left() # y axis on top
    chart.ax_heatmap.yaxis.set_label_position('left')

    # make it save/display in high-resolution
    fig = plt.gcf()
    fig.set_size_inches(10,8)

    # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,box_inches = 'tight')
            
            
#obtaining all the different combinations of CVDs
def get_combinations(my_list):
    '''
    This function gets every combination (in order) from my_list
    my_list is size n, this function returns combinations of sizes [2 ... n]
    '''
    combs = []
    for size in (range(2,len(my_list)+1)):
        combs += list(combinations(my_list,size))
    return combs

def get_bool_table(categories, ordering=None):
    label_ordering = []
    if ordering:
        label_ordering = ordering
    else:
        # get unique labels in alphabetical order
        unique_labels = set()
        for c in categories:
            for c_ in c:
                unique_labels.add(c)
        label_ordering = list(c)
        label_ordering.sort()

    # populate boolean arrays                                                                                                
    bool_arrays = []
    for l in label_ordering:
        b_a = []
        for c in categories:
            b = l in c
            b_a += [b]
        bool_arrays += [b_a]

    # multi index from labels and bool arrays
    multi_index = pd.MultiIndex.from_arrays(bool_arrays, names = label_ordering)
    return multi_index

'''
Pathway analysis
'''


def download_reactome_pathway_analysis_results(token, out_file = "results.csv"):
    '''
    This downloads the results from a Reactome Pathway Analysis via its token as out_file
    '''
    # prepare parameters for curl command
    headers = {
        'accept': '*/*',
    }
    url = 'https://reactome.org/AnalysisService/download/%s/pathways/TOTAL/result.csv'%(token)
    
    # post request
    response = requests.get(url, headers=headers)

    if "error" in response.text:
        print("Error occured")
        print(response.text)
        return
    
    # parse csv files, ignore comma within quotation marks
    data = [l for l in csv.reader(response.text.splitlines(), quotechar='"', delimiter=',',
                     quoting=csv.QUOTE_ALL, skipinitialspace=True)]

    # convert to pandas DataFrame, rename and drop empty columns
    df = pd.DataFrame(data)
    df = df.rename(columns=df.iloc[0]).drop(df.index[0])
    df = df.loc[:, df.columns.notna()]
    
    # save a csv
    if out_file:
        df.to_csv(out_file,index=False)
        print("Saved to file %s"%out_file)
    
    return df

def submit_reactome_pathway_analysis(protein_list, out_file='results.csv', debug=False):
    '''
    This is Reactome Pathway Analysis using project to human.
    Source: https://reactome.org/AnalysisService/#/identifiers/getPostTextToHuman
    '''
    # prepare parameters for curl command
    headers = {
        'accept': '*/*',
        'Content-Type': 'text/plain',
    }
    params = {
        'interactors': 'false',
        'pageSize': '20',
        'page': '1',
        'sortBy': 'ENTITIES_PVALUE',
        'order': 'ASC',
        'resource': 'TOTAL',
        'pValue': '1',
        'includeDisease': 'true',
    }

    # account for proteins separated by ';'
    proteins = []
    for p in protein_list:
        if ";" in p:
            proteins += p.split(";")
        else:
            proteins += [p]

    # prepare input data
    data = "\n".join(proteins)
    # post request
    response = requests.post('https://reactome.org/AnalysisService/identifiers/projection', 
                             params=params, headers=headers, data=data)
    # extract token
    result_dict = json.loads(response.text)
    token = result_dict['summary']['token']
    token = urllib.parse.quote(token) # format for url
    
    if debug:
        print("Sucessfully posted pathway analysis with token %s"%(result_dict['summary']['token']))
    
    # download result
    result_df = download_reactome_pathway_analysis_results(token, out_file)
    
    return result_df

def convert_uniref_ids_to_proteins(uniref_id_list, uniref_to_uniprot_list):
    plist = set()
#     for uniref_ids in list_uniref_ids.split("; "):
    for l in uniref_id_list:
        uniref_ids = l.split("; ")
        for uniref_id in uniref_ids:
#             print(uniref_id)
            proteins = uniref_to_uniprot_list[uniref_id]
            plist = plist.union(set(proteins))
    return plist


def extract_coverage(df, pathway_list):
    '''
    This returns the 'coverage' of a pathway based on proteins from a condition
    based on the pathways in pathway_list
    Calculate coverage using #Entities found / #Entities total
    '''

    # use the pathway ID's and get the rows for each cvd from this list
    sub_df = df[df['Pathway identifier'].isin(pathway_list)].copy()
    # Calculate coverage using #Entities found / #Entities total
    sub_df['coverage'] = sub_df['#Entities found'] / sub_df['#Entities total']
    
    # put as a dict, to extract in the correct order
    coverages_dict = dict(zip(sub_df['Pathway identifier'],sub_df['coverage']))
    # fill in missing data with 0
    for pid in pathway_list:
        if pid not in coverages_dict:
            coverages_dict[pid] = 0
    coverages = [coverages_dict[pid] for pid in pathway_list]
    
    return coverages

def convert_to_pathway_name(df, pathways_to_include):
    '''
    Convert pathway IDs to pathway name.
    Borrows a Reactome PA dataframe and returns the appropriate names in order
    '''
    id_to_name = dict(zip(df['Pathway identifier'], df['Pathway name']))
    return [id_to_name[p] for p in pathways_to_include]


def extract_heatmap_data(title_to_df, pval_cutoff = 0.05, max_num = 15):
    '''
    Calculates the coverage for each CVD, based on the top num_rows pathways
    found in the 'all' analysis. 
    num_rows: The number of Reactome Pathways we want to keep as top n from 'all'
    '''    

    # extract the pathway ID's for the first n rows from all protein analysis
    all_df = title_to_df['category_i']
#     print(all_df.head())
    pathways_to_include = list(all_df[all_df['Entities pValue'] < pval_cutoff]['Pathway identifier'].head(n=max_num))
    print(len(pathways_to_include))
    # assemble data for heatmap
    headers = [] # stores cvd names
    coverages = [] # stores the coverage values
    for df_title in title_to_df.keys():
        if "nonzero" in df_title: # only take those with nonzero value
            cvd = df_title.split("_")[0]
            coverage = extract_coverage(title_to_df[df_title], pathways_to_include)
            coverages += [coverage]
            headers += [cvd]
    # make into dataframe
    df = pd.DataFrame(coverages)
    df.columns = convert_to_pathway_name(all_df, pathways_to_include) # set row names
    df = df.T # invert
    df.columns = [h.split("_")[0] for h in headers] # set column names
    return df, pathways_to_include

def extract_pathway_to_scores(pathways_to_extract, uniref_to_score_df, 
                                                            reactome_results_df, uniprot_to_uniref, debug=False):
    '''
    '''
#     cvds = list(uniref_to_score_df.copy().drop('protein', axis=1).columns) 
    cvds = list(uniref_to_score_df.copy())
                
    # step 1: get matching proteins in each pathway
    pathway_to_proteins = {}
    for pathway in pathways_to_extract:
        sub_df = reactome_results_df[reactome_results_df['Pathway identifier'] == pathway]
        if sub_df.shape[0] > 0:
            temp1 = list(sub_df['Submitted entities found'])[0].split(";")
            pathway_to_proteins[pathway] = temp1
        else:
            pathway_to_proteins[pathway] = []
    # step 2: get matching unirefs in each pathway
    pathway_to_unirefs = {}
    for pathway in pathways_to_extract:
        proteins = pathway_to_proteins[pathway]
        unirefs = set()
        for p in proteins:
            if p in uniprot_to_uniref:
                uniref = uniprot_to_uniref[p]
                unirefs.add(uniref)
                # print("%s maps to %s"%(p,uniref))
            # else:
                # print("%s not in uniprot to uniref mapping"%p)
        pathway_to_unirefs[pathway] = list(unirefs)
    # print(pathway_to_unirefs)

    # # step 3: convert uniref_to_score_df to protein level df
    # protein_to_score_df = {}

    # step 4: extract relevant uniref_to_score_df subset for each pathway
    pathway_to_uniref_scores = {}
    for pathway, unirefs in pathway_to_unirefs.items():
        sub_df = uniref_to_score_df[uniref_to_score_df.index.isin(unirefs)]
        pathway_to_uniref_scores[pathway] = sub_df

    # # step 5: extract relevant protein_to_score_df subset for each pathway
    # pathway_to_protein_scores = {}

    # step 6: debug and return
    if debug:
        summary_df = {'Pathway ID':[],'Pathway size':[]}
        for cvd in cvds:
            summary_df[cvd] = [] 
        for pathway in pathways_to_extract:
            sub_df = reactome_results_df[reactome_results_df['Pathway identifier'] == pathway]
            pathway_size = list(sub_df['#Entities total'])[0]
            uniref_score_df = pathway_to_uniref_scores[pathway]
            # protein_score_df = pathway_to_protein_scores[pathway]
            summary_df['Pathway ID'] += [pathway]
            summary_df['Pathway size'] += [pathway_size]
            # print(uniref_score_df)
            for cvd in cvds:
                num_unirefs = uniref_score_df.shape[0] - (uniref_score_df[cvd].isna()).sum()
                temp_text = str(num_unirefs)
                # num_proteins = ...
                # temp_text = "(%d, %d)"%(num_unirefs,num_proteins)
                summary_df[cvd] += [temp_text]
        summary_df = pd.DataFrame(summary_df)
        print(summary_df)

    return pathway_to_uniref_scores


def get_uniref_to_uniprot_list(prot_to_uniref_file, reverse_mapping=False):
    # map uniref IDs back to proteins
    map_df = pd.read_csv(prot_to_uniref_file, sep = "\t")
    map_df.columns = ['UniProt list','UniRef90']
    uniref_to_uniprot_list = {}
    for upid_list, uniref_id in zip(map_df['UniProt list'], map_df['UniRef90']):
        uniref_to_uniprot_list[uniref_id] = upid_list.split(",")
    if not reverse_mapping:
        return uniref_to_uniprot_list
    else:
        uniprot_to_uniref_list = {}
        for uniref, uniprots in uniref_to_uniprot_list.items():
            for uniprot in uniprots:
                uniprot_to_uniref_list[uniprot] = uniref
        return uniprot_to_uniref_list


def extract_pathway_avg_zscore_matrix(pathway_to_uniref_scores):
    score_dict = {'Pathway':[]}
    for pathway, uniref_scores_df in pathway_to_uniref_scores.items():
        cvds = list(uniref_scores_df.copy().columns) 
#         cvds = list(uniref_scores_df.copy().drop('protein', axis=1).columns) 
        for cvd in cvds:
            if cvd not in score_dict:
                score_dict[cvd] = []
            avg_zscore = np.mean(uniref_scores_df[cvd])
            score_dict[cvd] += [avg_zscore]
        score_dict['Pathway'] += [pathway]
    z_data = pd.DataFrame(score_dict).set_index('Pathway')
    return z_data

def make_heatmap(data, v_lim = (0.0,1.0), ytick_max=60,
                                 out_file = "./output/CVD_Reactome_coverage_heatmap_no_dendrogram.pdf"):
    '''
    Makes a heatmap with dendrogram using seaborne and the prev function as data
    '''    

    # truncate the pathway names if its too long
    df = data.copy(deep=True)
    df.index = [s[:ytick_max] if (len(s) > ytick_max) else s for s in df.index]
    
    # plot the heamap
    chart = sns.clustermap(df, cmap='YlGnBu', annot=True, xticklabels=True, 
                                            yticklabels=True,col_cluster=False,row_cluster=False,
                                            vmin = v_lim[0], vmax = v_lim[1], fmt='.2f',cbar_kws={"shrink": .5})
    chart_axis = chart.ax_row_dendrogram.axes
    # dendrogram(linkage_matrix, ax=chart_axis, orientation='left', 
    #                         color_threshold=0, above_threshold_color='black')
    
    chart.ax_heatmap.xaxis.tick_top() # x axis on top
    chart.ax_heatmap.xaxis.set_label_position('top')

    # make it save/display in high-resolution
    fig = chart.fig
    fig.set_size_inches(15,5)

    # colorbar
    # fig.colorbar(ax.get_children()[0], orientation="horizontal")

    # # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,box_inches = 'tight')
    
    return chart

# reorder the heatmap rows to match the dendrogram
def reorder_heatmap(heatmap_df, new_heatmap_order_pid, pid_to_pname, pid_to_proteins):

    # stores only reactome pathway names, needed for re-ordering the hmap
    new_index_order = []
    # stores reactome pathway name and size, for renaming the hmap
    new_index_names = []
    for pid in new_heatmap_order_pid:
        pathway_id = pid.split("_")[0]
        pathway_name = pid_to_pname[pathway_id]
        proteins = pid_to_proteins[pathway_id]
        new_index_order += [pathway_name]
        new_index_names += ["%s (%d proteins)" % (pathway_name, len(proteins))]

    # re-order and re-name
    new_df = heatmap_df.copy()
    new_df = heatmap_df.reindex(new_index_order)
    new_df.index = new_index_names

    return new_df

def make_heatmap(data, linkage_matrix, v_lim = (0.0,1.0),
                                 out_file = "./output/dendrogram_CVD_Reactome_coverage_heatmap.pdf"):
    '''
    Makes a heatmap with dendrogram using seaborne and the prev function as data
    '''    

    # plot the heamap
    chart = sns.clustermap(data, cmap='YlGnBu', annot=True, xticklabels=True, 
                                            yticklabels=True,col_cluster=False,row_cluster=False,
                                            vmin = v_lim[0], vmax = v_lim[1], fmt='.2f',cbar_kws={"shrink": .5})
    chart_axis = chart.ax_row_dendrogram.axes
    dendrogram(linkage_matrix, ax=chart_axis, orientation='left', 
                            color_threshold=0, above_threshold_color='black')
    
    chart.ax_heatmap.xaxis.tick_top() # x axis on top
    chart.ax_heatmap.xaxis.set_label_position('top')

    # make it save/display in high-resolution
    fig = chart.fig
    fig.set_size_inches(15,10)

    # colorbar
    # fig.colorbar(ax.get_children()[0], orientation="horizontal")

    # # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,box_inches = 'tight')
    
    return chart



# source: https://stackoverflow.com/questions/63530701/python-package-to-plot-two-heatmaps-in-one-split-each-square-into-two-triangles/63531813#63531813
def double_heatmap(c_data,z_data,dend_data):
        # sanity check
        if c_data.shape != z_data.shape:
                print("Data frames are not the same shape!")
                return
        n_pathways = c_data.shape[0]
        n_cvds = c_data.shape[1]
        
        # make coordinates
        M = n_cvds
        N = n_pathways
        x = np.arange(M + 1)
        y = np.arange(N + 1)
        xs, ys = np.meshgrid(x, y)
        # triangles1 = [(i + j*(M+1), i+1 + j*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
        # triangles2 = [(i+1 + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
        # I flipped the triangles the other way
        triangles1 = [(i + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
        triangles2 = [(i+1 + j*(M+1), i+1 + (j+1)*(M+1), i + (j)*(M+1)) for j in range(N) for i in range(M)]
        triang1 = Triangulation(xs.ravel(), ys.ravel(), triangles1)
        triang2 = Triangulation(xs.ravel(), ys.ravel(), triangles2)

        # extract intensities
        c_zs = []
        z_zs = []
        for i in range(N):
            idx = N-i-1
            c_zs += list(c_data.iloc[idx,:])
            z_zs += list(z_data.iloc[idx,:])

        # plot
        img1 = plt.tripcolor(triang1, c_zs, cmap='Reds', vmax=1)
        img2 = plt.tripcolor(triang2, z_zs, cmap='Blues')

        fig, ax = plt.subplots()
        plt.colorbar(img2, ticks=range(10), pad=-0.05)
        plt.colorbar(img1, ticks=range(10))
        plt.xlim(x[0], x[-1])
        plt.ylim(y[0], y[-1])
        plt.xticks(x, rotation=90)
        plt.yticks(y)
        plt.show()

        fig.savefig("./output/temp.pdf",box_inches = 'tight')
    
        # print(z_zs)
        return
    
    

# source: https://stackoverflow.com/questions/63530701/python-package-to-plot-two-heatmaps-in-one-split-each-square-into-two-triangles/63531813#63531813
def double_heatmap(c_data,z_data,dend_data, ytick_max=60,
                                 out_file="./output/double_heatmap_figure.pdf"):
    # sanity check
    if c_data.shape != z_data.shape:
            print("Data frames are not the same shape!")
            return
    n_pathways = c_data.shape[0]
    n_cvds = c_data.shape[1]

    # make coordinates
    M = n_cvds
    N = n_pathways
    x = np.arange(M + 1)
    y = np.arange(N + 1)
    xs, ys = np.meshgrid(x, y)
    # triangles1 = [(i + j*(M+1), i+1 + j*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
    # triangles2 = [(i+1 + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
    # I flipped the triangles the other way
    triangles1 = [(i + j*(M+1), i+1 + (j+1)*(M+1), i + (j+1)*(M+1)) for j in range(N) for i in range(M)]
    triangles2 = [(i+1 + j*(M+1), i+1 + (j+1)*(M+1), i + (j)*(M+1)) for j in range(N) for i in range(M)]
    triang1 = Triangulation(xs.ravel(), ys.ravel(), triangles1)
    triang2 = Triangulation(xs.ravel(), ys.ravel(), triangles2)

    # extract intensities
    c_zs = []
    z_zs = []
    for i in range(N):
        idx = N-i-1
        c_zs += list(c_data.iloc[idx,:])
        z_zs += list(z_data.iloc[idx,:])

    # make two subplots
    f, axes = plt.subplots(2, 2, gridspec_kw={'width_ratios': [1, 3],
                                               'height_ratios':[20,1]},
                                               figsize=(15, 9))
    a0 = axes[0][0]
    a1 = axes[0][1]
    a2 = axes[1][0]
    a3 = axes[1][1]   
    
    ### dendrogram ###
    dendrogram(dend_data, ax=a0, orientation='left',
                        color_threshold=0, above_threshold_color='black') 

    # plot
    img1 = a1.tripcolor(triang1, c_zs, cmap='Reds', vmax=1)
    img2 = a1.tripcolor(triang2, z_zs, cmap='Blues')
    a1.set_xlim(x[0], x[-1])
    a1.set_ylim(y[0], y[-1])

#     plt.colorbar(img2, ticks=range(10), orientation='horizontal', pad=-2.0)
#     plt.colorbar(img1, ticks=range(10), orientation='horizontal', pad=-3.0)
#     plt.colorbar(img2, ticks=range(10), orientation='horizontal', ax=a3, aspect=20,location='bottom')
#     plt.colorbar(img1, ticks=range(10), orientation='horizontal', ax=a3, aspect=20,location='top')
    cbar1 = plt.colorbar(img1, ticks=range(10), orientation='horizontal', ax=a3, aspect=40, fraction=0.7, pad=0.04,location='top', panchor=(1.0,0.5))
    cbar2 = plt.colorbar(img2, ticks=range(10), orientation='horizontal', ax=a3, aspect=40, fraction=2.1, pad=0.04,location='bottom', panchor=(1.0,0))
    cbar1.set_label("c-score")
    cbar2.set_label("Average z-score")
    ### formatting ###
    # truncate columns
    truncated_ylabels = [s[:ytick_max] if (len(s) > ytick_max) else s for s in list(c_data.index)[::-1]]
    
    # remove figure borders
    a0.spines['top'].set_visible(False)
    a0.spines['right'].set_visible(False)
    a0.spines['bottom'].set_visible(False)
    a0.spines['left'].set_visible(False)
    a1.spines['top'].set_visible(False)
    a1.spines['right'].set_visible(False)
    a1.spines['bottom'].set_visible(False)
    a1.spines['left'].set_visible(False)    
    a0.patch.set_facecolor('white')
    a1.patch.set_facecolor('white')

    # disable a2 and a3
    a2.axis('off')
    a3.axis('off')
    
    # Turn off tick labels in dendrogram
    a0.set_yticks([])
    a0.set_xticks([])
    # re-label and move axes for bubble plot
    a1.yaxis.tick_right()
    a1.xaxis.tick_top()
    a1.set_xticks([x+0.5 for x in list(range(n_cvds))])
    
    a1.set_xticklabels(c_data.columns)
    a1.set_yticks([y+0.5 for y in list(range(n_pathways))])
    a1.set_yticklabels(truncated_ylabels)

    if len(out_file) > 0:
        plt.tight_layout()
        f.savefig(out_file,box_inches = 'tight')

    return
    
def add_proteins_to_uniref_scores(pathway_to_uniref_scores, 
                                reactome_results_df, uniref_to_uniprot_list,
                                pathway_id_to_pathway_name,
                                reactome_pathway_to_unique_proteins):
    '''
    This function adds two column to each dataframe in pathway_to_uniref_scores
    pathway_to_uniref_scores is a dict mapping pathway -> df showing uniref CaseOLAP scores
    The new function columns will rename 'protein' to 'ProteinGroup', and add
    two new columns (1) the list of mammalian proteins in that protein group 
    according to uniref_to_uniprot_list, (2) the human proteins that were 
    identified to be in the pathway based on Reactome Pathway Analysis
    '''
    # identify which proteins were mapped by reactome
    pathway_to_submitted_entities_found = {}
    pathway_to_mapped_entities = {}
    for pathway,sub_ent,mapped_ent in zip(reactome_results_df['Pathway identifier'], 
                                    reactome_results_df['Submitted entities found'],
                                    reactome_results_df['Mapped entities']):
        # pathway_to_mapped_proteins[pathway] = mapped_ent.split(";")
        pathway_to_submitted_entities_found[pathway] = sub_ent
        pathway_to_mapped_entities[pathway] = mapped_ent
        
    # add the columns to each df
    pathway_to_uniref_proteins_scores = {}
    for pathway in pathway_to_uniref_scores.keys():
        df = pathway_to_uniref_scores[pathway].copy()
        # rename 'protein' to 'ProteinGroup'
        df.columns = ['ProteinGroup'] + list(df.columns[1:])

        # "What proteins are in this protein group?"
        # add mammalian proteins to df
        mammalian_proteins = []
        for uniref_id in df['ProteinGroup']:
            proteins_in_protein_group = uniref_to_uniprot_list[uniref_id]
            mammalian_proteins += [";".join(proteins_in_protein_group)]
        df['Proteins'] = mammalian_proteins
        
        file_output_name = "(%s) %s.csv"%(pathway,pathway_id_to_pathway_name[pathway])
        # reorder columns
        df = df[[df.columns[0]] + [df.columns[-1]] + list(df.columns[1:-1])]

        pathway_to_uniref_proteins_scores[file_output_name] = df

    return pathway_to_uniref_proteins_scores

def save_pathway_to_proteins(pathways, pathway_id_to_pathway_name,
                             reactome_pathway_to_unique_proteins,
                             out_file_name = "./output/heatmap_data/pathway_to_proteins.csv"):
    with open(out_file_name,"w") as out_file:
        header = ['PathwayID','PathwayName','Proteins']
        out_file.write(",".join(header)+"\n")

        for pathway_id in pathways:
            pathway_name = pathway_id_to_pathway_name[pathway_id]
            proteins = reactome_pathway_to_unique_proteins[pathway_id]
            proteins_text = ";".join(proteins)

            out_file.write(",".join([pathway_id,pathway_name,proteins_text])+"\n")
    
# save heatmap data
def save_heatmap_data(c_data,z_data,pathway_to_uniref_proteins_scores, 
                      out_folder="./output/heatmap_data/"):
    c_data.to_csv(out_folder+"c_data.csv",index=False)
    z_data.to_csv(out_folder+"z_data.csv",index=False)

    for out_name, df in pathway_to_uniref_proteins_scores.items():
        df.to_csv(out_folder + out_name,index=False)
        

def make_heatmap_unique_to_cvd(unique_to_cvd_pathway_list, zscores_df, 
                                                     all_union_pa_df, uniprot_to_uniref, cvds,
                                                     pathway_id_to_pathway_name, pid_to_proteins,
                                                     max_pathway_name_length=60,
                                 out_file = "./pathways_unique_to_cvd_heatmap.pdf"):
    '''
    Makes a heatmap with dendrogram using seaborne and the prev function as data
    '''    

    # extract the data
    pathways = []
    # cvds = unique_to_cvd_pathway_list.keys()
    for cvd in cvds:
        if cvd in unique_to_cvd_pathway_list:
            pathway_list = unique_to_cvd_pathway_list[cvd]
            for pathway in pathway_list:
                pathways+=[pathway]
        
    pathway_to_scores = extract_pathway_to_scores(pathways, zscores_df, 
                                                        all_union_pa_df, uniprot_to_uniref, debug=False)
    pathway_names = []

    for pathway_id in pathways:
        # truncate the pathway name if too long
        pathway_name = pathway_id_to_pathway_name[pathway_id]
        proteins = pid_to_proteins[pathway_id]
        p_name = pathway_name[:max_pathway_name_length]
        # add elipses
        if len(pathway_name) > max_pathway_name_length:
            p_name = p_name[:-4]+" ..."
        # add protein size of pathway
        p_name = "%s (%d proteins)"%(p_name,len(proteins))
        pathway_names += [p_name]

    data = extract_pathway_avg_zscore_matrix(pathway_to_scores)
    data.index = pathway_names
    
    # plot the heamap
    chart = sns.clustermap(data, cmap='YlGnBu', annot=True, xticklabels=True, 
                                            yticklabels=True,col_cluster=False,row_cluster=False,
                                            fmt='.2f',cbar_kws={"shrink": .5}, 
                                            cbar_pos = (0.045, 0.04, 0.02, 0.76))
    
    chart.ax_heatmap.xaxis.tick_top() # x axis on top
    chart.ax_heatmap.xaxis.set_label_position('top')

    chart.ax_heatmap.set_xticklabels(chart.ax_heatmap.get_xmajorticklabels(), weight='bold')
    chart.ax_heatmap.set_yticklabels(chart.ax_heatmap.get_ymajorticklabels(), weight='bold')



    # make it save/display in high-resolution
    fig = chart.fig
    fig.set_size_inches(30,10)

    # colorbar
    # fig.colorbar(ax.get_children()[0], orientation="horizontal")

    # # save file if file name is given
    if len(out_file) > 0:
            fig.savefig(out_file,box_inches = 'tight')
    
    return chart
        
'''
Reactome stuff
'''
# Imported functions from Alex's github repo
# source: https://github.com/arpelletier/dglke_workspace/blob/main/bin/kg_utils.py

'''
This function reads the hierarchy information from Reactome and only returns the pathway ID's
which are for humans (HSA). These relationships are represented as an edge between a pair of nodes,
reported as (ancestor_pathway, decendent_pathway). Although this information isn't typically used,
these information is important to identify root pathways (pathways with no ancestors) and leaf
pathways (pathways with no decendents) by looking at the in/out degree with a directed graph.
This function also returns the set of unique reactome ids found in these relationships.
'''


def extract_human_hierarchical_information(reactome_hierarchy_to_pathway_file):
    # Read hierarchical information
    reactome_hierarchical_information = [f.strip("\n").split("\t")
                                         for f in open(reactome_hierarchy_to_pathway_file, 'r').readlines()]

    human_reactome_hierarchical_info = []
    for v1, v2 in reactome_hierarchical_information:
        if "HSA" in v1 and "HSA" in v2:
            human_reactome_hierarchical_info += [(v1, v2)]

    # Describe the data
    unique_human_reactome_ids = set()
    for v1, v2 in human_reactome_hierarchical_info:
        unique_human_reactome_ids.add(v1)
        unique_human_reactome_ids.add(v2)

    print("Number of edges: " + str(len(human_reactome_hierarchical_info)))
    print("Number of nodes: " + str(len(unique_human_reactome_ids)))
    return human_reactome_hierarchical_info, unique_human_reactome_ids


'''
This function reports the set of proteins corresponding to a pathway_id found in the input file. This
information is found in the reactome_uniprot_to_pathway_file table downloaded from the Reactome
website. Exclude isoforms removes protein accession which correspond to a protein isoform, reporting
only the cannonical protein accession instead.
'''


def extract_pathway_to_proteins(reactome_uniprot_to_pathway_file,
                                exclude_isoforms=True):
    # Load file for mapping reactome id -> proteins
    reactome_uniprot_to_pathway_table = pd.read_csv(reactome_uniprot_to_pathway_file, sep='\t', header=None)
    reactome_uniprot_to_pathway_table.columns = ["Uniprot", "Pathway identifier", "URL", "Pathway name", "X",
                                                 "Organism"]

    # add all proteins to each pathway set
    pathway_to_unique_proteins = {}
    for protein, pathway in zip(reactome_uniprot_to_pathway_table['Uniprot'],
                                reactome_uniprot_to_pathway_table['Pathway identifier']):
        # initialize set for each pathway
        if pathway not in pathway_to_unique_proteins:
            pathway_to_unique_proteins[pathway] = set()
        pathway_to_unique_proteins[pathway].add(protein)

    # remove protein isoforms
    if exclude_isoforms:
        for pathway in pathway_to_unique_proteins.keys():

            # identify isoforms
            isoforms = set()
            for protein in pathway_to_unique_proteins[pathway]:
                if "-" in protein:
                    isoforms.add(protein)

            # remove isoforms and add base protein
            for protein in isoforms:
                base_protein = protein.split("-")[0]
                pathway_to_unique_proteins[pathway].remove(protein)
                pathway_to_unique_proteins[pathway].add(base_protein)

    # Information about what we gathered
    num_proteins_per_pathway = []
    for pathway in pathway_to_unique_proteins.keys():
        proteins = pathway_to_unique_proteins[pathway]
        num_proteins_per_pathway += [len(proteins)]
    num_proteins_frequency = {x: num_proteins_per_pathway.count(x) for x in num_proteins_per_pathway}
    num_pathways_with_non_zero_proteins = 0
    for i in num_proteins_frequency.keys():
        freq = num_proteins_frequency[i]
        if i > 0:
            num_pathways_with_non_zero_proteins += freq
    arr = np.array(num_proteins_per_pathway)
    mean_num_proteins = np.mean(arr)
    median_num_proteins = np.median(arr)
    stdev_num_proteins = np.std(arr)

    # Describe the data
    print("Number of pathways extracted: " + str(len(pathway_to_unique_proteins)))
    print("Number of pathways extracted with > 0 proteins: " + str(num_pathways_with_non_zero_proteins))
    print("Mean pathway size: " + str(mean_num_proteins))
    print("Median pathway size: " + str(median_num_proteins))
    print("Standard deviation pathway size: " + str(stdev_num_proteins))

    return pathway_to_unique_proteins


'''
Returns a dict mapping pathway id -> pathway name
'''


def extract_pathway_id_to_pathway_name(file_name):
    reactome_id_to_pathway_name_lines = [f.strip("\n").split("\t") for f in
                                         open(file_name, 'r').readlines()]
    pathway_id_to_name = {}
    for line in reactome_id_to_pathway_name_lines:
        pathway_id_to_name[line[0]] = line[1]
    return pathway_id_to_name


'''
This function returns a list of pathways which are the descendants of each pathway in ancestor_list. sort_by_size will sort
the list based on the number of descendants it has. If largest_first, the sorted list will be in the order of largest number
of descendants first, otherwise it will be smallest first. Each pathway is added to the list of its ancestor pathways in the
specified order. If remove_duplicates is true, each descendant will only show up in the first ancestor list it has appeared.
After removing duplicates, the lists will not necessarily be sorted, as the size of each pathway has changed upon removing
duplicate nodes.
'''


def get_descendant_pathways(ancestor_list, G, sort_by_size=True,
                            largest_first=True,
                            remove_duplicates=True):
    pair_nodes = []
    for root_node in ancestor_list:
        bfs_traversal = ([root_node] +
                         [v1 for v1, v2 in
                          list(nx.algorithms.traversal.breadth_first_search.bfs_predecessors(G, root_node))])
        # print(bfs_traversal)
        # print(root_node+" has "+str(len(bfs_traversal))+" nodes in its pathway")
        pair_nodes += [(bfs_traversal, len(bfs_traversal))]

    if sort_by_size:
        pair_nodes.sort(key=lambda x: x[1])

        if largest_first:
            pair_nodes = pair_nodes[::-1]

    filtered_list = []
    if remove_duplicates:
        used_pathways = set()
        for nodes, size in pair_nodes:
            node_list = []
            for node in nodes:
                if node not in used_pathways:
                    node_list += [node]
                    used_pathways.add(node)
            if len(node_list) > 0:
                filtered_list += [(node_list, len(node_list))]
            else:
                # it seems there is a parent pathway which is labeled by another parent pathway, so...
                # this is redundant...
                filtered_list += [([nodes[0]], 0)]
    else:
        filtered_list = pair_nodes

    # want to have root path in position 0 and list in position 1
    return_list = []
    for nodes, size in filtered_list:
        # print(pathway_id_to_pathway_name[nodes[0]] + " has " + str(size) + " nodes")
        return_list += [(nodes[0], nodes)]
    return return_list


'''
Perform this function after calling get_descendant_pathways. It creates a reverse mapping of descendant -> ancestor.
For this to work properly, you should remove the duplicate pathways such that one pathway only shows up in one ancestor list.
'''


def get_familiar_relations(descendants):
    relations = {}
    for family, members in descendants:
        for member in members:
            relations[member] = family
    return relations



'''
This function prunes edges of the DAG, breaking ties with a labeling function. Edges between a descendant and ancestor which
share the same label will be preserved. If no labeling function is provided or if there are multiple ancestors with the same 
label, then the ancestor which has a larger subtree (larger number of total descendants) is retained. All other descendant-
ancestor edges are pruned from the DAG. The resulting structure should be a tree.
'''


def prune_extra_edges(G, labeling_function=None, debug=False):
    nodes = G.nodes()
    leaves = set(n for n in nodes if G.out_degree(n) == 0)
    roots = set(n for n in nodes if G.in_degree(n) == 0)
    inner_nodes = [n for n in nodes if G.out_degree(n) > 0]

    # Compute the direct decendants and ancestors of each node
    descendants = dict()
    ancestors = dict()
    for ancestor in dict(G.adjacency()).keys():
        children_dict = dict(G.adjacency())[ancestor]
        children = list(children_dict.keys())
        descendants[ancestor] = children
        for child in children:
            if child not in ancestors:
                ancestors[child] = []
            if ancestor not in ancestors[child]:
                ancestors[child] += [ancestor]

    if (debug):
        print("Descendants: ")
        print(descendants)
        print("Ancestors: ")
        print(ancestors)

    # Compute the size of each subtree number of nodes
    subtree = dict((n, [n]) for n in leaves)
    for u in inner_nodes:
        children = set()
        nodes = list(descendants[u])
        while len(nodes) > 0:
            v = nodes.pop(0)
            children.add(v)
            nodes += descendants[v]
        subtree[u] = sorted(children)

    if (debug):
        print("Subtree: ")
        print(subtree)

    # identify those in the tree with extra edges
    extra_edges_nodes = [x for x in G.nodes() if G.in_degree(x) > 1]

    if (debug):
        print("Number of nodes with extra edges: " + str(len(extra_edges_nodes)))

    # make a copy of graph, add root node
    pruned_G = G.copy()
    for node in roots:
        pruned_G.add_edge("Root", node)

    # prune edges from extra_edges
    num_edges_pruned = 0
    for node in extra_edges_nodes:

        # there are necessarily more than one ancestor
        parents = ancestors[node]
        # store the node name and subtree size of the ancestor you want to keep
        keep_this_parent = ("", 0)

        # extract node label if labeling function exists
        node_label = None
        if labeling_function != None and node in labeling_function:
            node_label = labeling_function[node]

        if (debug):
            print(node)
            print(node_label)

        # determine which parent to keep
        for parent in parents:

            # extract parent label if labeling function exists
            parent_label = None
            if labeling_function != None and parent in labeling_function:
                parent_label = labeling_function[parent]

            if (debug):
                print(parent)
                print(parent_label)

            # keep the parent if it matches the labeling function of node
            if node_label == parent_label:
                subtree_size = len(subtree[parent])

                # keep the parent with largest subtree
                if subtree_size > keep_this_parent[1]:
                    keep_this_parent = (parent, subtree_size)

        # Rare case where parent label does not match node label
        if len(keep_this_parent[0]) < 1:
            # keep the parent pathway with the laregest subtree
            for parent in parents:
                subtree_size = len(subtree[parent])
                if subtree_size > keep_this_parent[1]:
                    keep_this_parent = (parent, subtree_size)

        # prune the edges to parents we aren't keeping
        if (debug):
            print("Decided to keep " + keep_this_parent[0] + ", with subtree size " + str(keep_this_parent[1]))

        for parent in parents:
            if parent != keep_this_parent[0]:
                # remove edge from graph
                pruned_G.remove_edge(parent, node)
                num_edges_pruned += 1

    print("Pruned " + str(num_edges_pruned) + " edges.")

    return pruned_G


'''
For every node in the node_list, this function computes the descendant->ancestor->root path on the graph.
Returned value is a dictionary of node -> path.
'''


def path_to_root(G, node_list, debug=False):
    # check to see if it is a tree
    if not nx.is_tree(G):
        print("Graph is not a tree! Please prune the unwanted edges first!")
        return
    elif debug:
        print("Graph is a tree.")

    rev_G = G.reverse()
    paths = dict()
    for node in node_list:
        bfs = list(nx.algorithms.traversal.breadth_first_search.bfs_successors(rev_G, node))
        path = [node]
        for descendant_node, parent_nodes in bfs:
            # take only the first parent node
            parent = parent_nodes[0]

            # ignore the paths corresponding to a parent node we did not select
            # this was originally to account for nodes with multiple parents, but if you pruned the tree
            # first, this is not an issue.
            if descendant_node in path:
                path += [parent]
        paths[node] = path
    return paths


'''
This function takes a sub-graph of the original graph returning only the nodes in the node_list and all of its
ancestor nodes, up to and including the root node. As well, any intermediate (non-leaf) nodes in the node list are
represented with a leaf node, flagged by "_leaf" suffix. The input graph should be a tree, not a DAG. 
'''


def construct_hierarchy_subgraph(G, node_list, debug=False):
    # node path to root for every node in node_list
    node_list_paths_to_root = path_to_root(G, node_list)

    # preserve the node_list but replace the "_leaf" labeled ones
    new_node_list = node_list.copy()

    # add every pair of edges to the new graph
    sub_G = nx.DiGraph()
    edges = []
    for node in node_list:
        path = node_list_paths_to_root[node]

        for i in range(len(path) - 1):
            to_node = path[i]
            from_node = path[i + 1]
            edges += [(from_node, to_node)]

        # add edges for daughter nodes which are parents of another daughter node
        for node in path[1:]:
            leaf_node = node + "_leaf"
            if node in node_list and leaf_node not in new_node_list:
                edges += [(node, leaf_node)]
                # replace name in new node list
                new_node_list = [n.replace(node, leaf_node) for n in new_node_list]
    sub_G.add_edges_from(edges)

    return sub_G, new_node_list


'''
Constructs a linkage matrix for a hierarchy tree. The input graph must be a DAG. The leaves of the dendrogram are specified
by the node_list. The graph is first pruned using a labeling function if provided, converting the graph into a tree if it
is not already. A sub-tree is then extracted, keeping only the necessary nodes in node_list and all of its ancestors. If a
node in the node_list is not a leaf node, then a new node will be added to make it so. This is necessary to construct the
dendrogram.
'''


def construct_dendrogram(G, node_list, labeling_function=None, debug=False):

    if not nx.is_tree(G):
        # prune the graph into a tree
        pruned_G = prune_extra_edges(G, labeling_function=labeling_function)

        # extract only relevant nodes from graph
        sub_G, leaves = construct_hierarchy_subgraph(pruned_G, node_list)
    else:
        sub_G, leaves = construct_hierarchy_subgraph(G, node_list)

    nodes = sub_G.nodes()
    inner_nodes = [n for n in nodes if sub_G.out_degree(n) > 0]

    if (debug):
        print("Number of nodes in node_list: " + str(len(node_list)))
        print("Number of leaf nodes in subgraph: " + str(len(leaves)))
        print("Number of nodes in node_list but not in leaves: " + str(len(set(node_list).difference(set(leaves)))))
        print("Number of nodes in leaves but not in node_list: " + str(len(set(leaves).difference(set(node_list)))))

    # Compute the direct decendants of each node
    descendants = dict()
    for key in dict(sub_G.adjacency()).keys():
        children_dict = dict(sub_G.adjacency())[key]
        children = list(children_dict.keys())
        descendants[key] = children
    if (debug):
        print(descendants)

    # Compute the size of each subtree number of nodes
    subtree = dict((n, [n]) for n in leaves)
    # Compute the size (i.e. how many leaf nodes) of each subtree
    leaf_subtree = dict((n, [n]) for n in leaves)
    for u in inner_nodes:
        children = set()
        nodes = list(descendants[u])
        while len(nodes) > 0:
            v = nodes.pop(0)
            children.add(v)
            nodes += descendants[v]
        subtree[u] = sorted(children)
        leaf_subtree[u] = sorted(children & set(leaves))
    if (debug):
        print("Subtree: ")
        print(subtree)
        print("Leaf subtree: ")
        print(leaf_subtree)

    # Compute the depth of each node in subtree
    depth_from_root = nx.shortest_path_length(sub_G, "Root")
    max_depth = max(depth_from_root.values())
    height_of_dendrogram = {node: (max_depth - depth) for node, depth in
                            zip(depth_from_root.keys(), depth_from_root.values())}

    # Compute the linkage matrix
    inner_nodes.sort(key=lambda n: len(subtree[n]))  # <-- order inner nodes ascending by subtree size, root is last
    if (debug):
        print("inner nodes:")
        print(inner_nodes)
    node_to_index = {}
    num_unique_indexes = 0
    for node in leaves:
        node_to_index[node] = num_unique_indexes
        num_unique_indexes += 1

    if (debug):
        print(node_to_index)

    linkage_matrix = []
    for ancestor_node in inner_nodes:
        if (debug):
            print("Ancestor node: " + ancestor_node)

        children = descendants[ancestor_node]
        if len(children) > 1:
            # merge the branches if there are more than 1 children
            for j in range(len(children) - 1):
                child_1 = children[j]
                child_2 = children[j + 1]
                child_1_index = node_to_index[child_1]
                child_2_index = node_to_index[child_2]
                node_to_index[ancestor_node] = num_unique_indexes
                linkage_matrix += [[float(child_1_index),
                                    float(child_2_index),
                                    float(len(subtree[ancestor_node])),
                                    float(len(leaf_subtree[ancestor_node]))]]
                if (debug):
                    print("Combine " + child_1 + " and " + child_2 + " into ancestor node " + ancestor_node +
                          ", index: " + str(node_to_index[ancestor_node]))
                node_to_index[child_1] = node_to_index[ancestor_node]
                node_to_index[child_2] = node_to_index[ancestor_node]
                num_unique_indexes += 1
        else:
            # otherwise, the ancestor inherits the index of its single child
            node_to_index[ancestor_node] = node_to_index[children[0]]
            if (debug):
                print("Ancestor node " + ancestor_node + " inherits the index of its single child "
                      + children[0] + ", index: " + str(node_to_index[ancestor_node]))

    # Visualize
    dend = dendrogram(linkage_matrix, labels=leaves, leaf_rotation=90)
    dendrogram_leaf_order = dend['ivl'][::-1]
    plt.show()

    return linkage_matrix, dendrogram_leaf_order


