
import os, sys
from analysis._03_text_mining_association_analysis import *

def parse_id_to_synonym(entity_count_files):
    # prepare ranked_ent_df
    ranked_ent ={}
    for f in entity_count_files:
        tag = f.split("_")[-1].strip(".json")
        ent_dict = json.load(open(f,'r'))
        ranked_ent[tag] = ent_dict
    ranked_ent_df = pd.DataFrame(ranked_ent)
    return prepare_synonyms(ranked_ent_df)


def get_required_files(root_directory,merge_proteins, use_core_proteins):
    '''
    Required files: caseolap.csv, 3 reactome files, ent2syn2count.json (some sort of it)
    :param root_directory:
    :return:
    '''
    data_folder = os.path.join(root_directory,'data')

    ### input files ###
    caseolap_result_folder = os.path.join(root_directory,"output/text_mining_results/all_proteins/")
    if use_core_proteins:
        caseolap_result_folder = os.path.join(root_directory, "output/text_mining_results/core_proteins/")

    caseolap_results_file = os.path.join(caseolap_result_folder,"caseolap.csv")

    reactome_uniprot_to_pathway_file = os.path.join(data_folder,"Reactome/UniProt2Reactome_All_Levels.txt")
    reactome_hierarchy_to_pathway_file = os.path.join(data_folder,"Reactome/ReactomePathwaysRelation.txt")
    reactome_id_to_pathway_name_file = os.path.join(data_folder,"Reactome/ReactomePathways.txt")

    # synonym related files for redundancy removal
    entity_count_files = []
    if merge_proteins:
        entity_count_folder = os.path.join(caseolap_result_folder,'ranked_proteins/ranked_caseolap_score')

        for f in os.listdir(entity_count_folder):
            if 'dictionary_ranked_proteins' in f:
                entity_count_files += [os.path.join(entity_count_folder,f)]

    ret = {'caseolap_results':caseolap_results_file,
           'reactome_data':[reactome_uniprot_to_pathway_file,
                            reactome_hierarchy_to_pathway_file,
                            reactome_id_to_pathway_name_file],
           'entity_count_files':entity_count_files}

    # check files
    for k,v in ret.items():
        if type(v) == list:
            files_to_check = v
        else:
            files_to_check = [v]
        for f in files_to_check:
            if not os.path.isfile(f):
                print("Input file missing! %s missing at path %s"%(k,f))
                return

    return ret



def load_files(root_directory, merge_proteins, use_core_proteins):

    input_files = get_required_files(root_directory, merge_proteins, use_core_proteins)
    if not input_files:
        print("Error, missing file! Make sure the previous steps ran properly")
        sys.exit()

    # load CaseOLAP scores
    raw_caseolap_scores = pd.read_csv(input_files['caseolap_results'])

    # load Reactome data
    reactome_uniprot_to_pathway, reactome_hierarchy_to_pathway, reactome_id_to_pathway_name = input_files['reactome_data']
    # Load in the reactome hierarchial information for human. See function definition in first cell.
    hierarchical_relationships, unique_reactome_ids = \
        extract_human_hierarchical_information(reactome_hierarchy_to_pathway)

    # Load in pathway id to proteins relationships. See function definition in first cell.
    reactome_pathway_to_unique_proteins = extract_pathway_to_proteins(reactome_uniprot_to_pathway,
                                                                      exclude_isoforms=False)

    # Load in mapping of reactome pathway id to pathway name
    pathway_id_to_pathway_name = extract_pathway_id_to_pathway_name(reactome_id_to_pathway_name)

    reactome_data = [hierarchical_relationships, unique_reactome_ids,reactome_pathway_to_unique_proteins, pathway_id_to_pathway_name]

    # load id_to_synonym
    id_to_synonym = parse_id_to_synonym(input_files['entity_count_files'])

    return raw_caseolap_scores, reactome_data, id_to_synonym


def run_pathway_analysis(summary_table):
    ### cvd subproteome protein lists ###
    print("### CVD Sub-proteomes ###")
    cvd_subproteome_unirefs = {}
    cvd_subproteome_proteins = {}
    for cvd in cvds:
        tag = "%s Sub-proteome" % (cvd)
        proteins = list(summary_table[summary_table[tag]].index)
        #     unirefs = list(summary_table[summary_table[tag]].index)
        #     proteins = convert_uniref_ids_to_proteins(unirefs,uniref_to_uniprot_list)
        print("%s: %d proteins" % (tag, len(proteins)))

        cvd_subproteome_unirefs[cvd] = unirefs
        cvd_subproteome_proteins[cvd] = proteins

    ### category i protein lists ###
    print("### Category I ###")
    category_i_proteins = list(summary_table[summary_table['Category I']].index)
    # category_i_unirefs = list(summary_table[summary_table['Category I']].index)
    # category_i_proteins = convert_uniref_ids_to_proteins(category_i_unirefs,uniref_to_uniprot_list)
    print("%d proteins in Category I" % len(category_i_proteins))

    ### category ii protein lists ###
    print("### Category II ###")
    category_ii_unirefs = {}
    category_ii_proteins = {}
    for cvd in cvds:
        tag = "%s Category II" % (cvd)
        proteins = list(summary_table[summary_table[tag]].index)
        #     unirefs = list(summary_table[summary_table[tag]].index)
        #     proteins = convert_uniref_ids_to_proteins(unirefs,uniref_to_uniprot_list)
        print("%s: %d proteins" % (tag, len(proteins)))

        category_ii_unirefs[cvd] = unirefs
        category_ii_proteins[cvd] = proteins

    for cvd in cvds:
        proteins = cvd_subproteome_proteins[cvd]

        if len(proteins) > 0:
            out_name = "./reactome_results/%s_subproteome_pathway_analysis_result.csv" % cvd
            pwa_df = submit_reactome_pathway_analysis(proteins, out_file=out_name)
            num_pathways = pwa_df.shape[0]
            print("%d pathways identified. Saved to %s" % (num_pathways, out_name))
        else:
            print("Insufficient number of proteins for %s" % cvd)

    # run pathway analysis for category i proteins
    cat_i_pa_result_df = submit_reactome_pathway_analysis(category_i_proteins, debug=True,
                                                          out_file="./reactome_results/category_i_proteins_pathway_analysis_result.csv")
    cat_i_pa_result_df.head()

    category_ii_pa_df = {}
    for cvd in cvds:
        proteins = category_ii_proteins[cvd]

        if len(proteins) > 0:
            out_name = "./reactome_results/%s_category_ii_pathway_analysis_result.csv" % cvd
            pwa_df = submit_reactome_pathway_analysis(proteins, out_file=out_name)
            num_pathways = pwa_df.shape[0]
            print("%d pathways identified. Saved to %s" % (num_pathways, out_name))
        else:
            print("Insufficient number of proteins for %s" % cvd)

    # all proteins?
    # all_unirefs = list(uniref_to_uniprot_list.keys())
    # all_proteins = set()
    # for proteins in uniref_to_uniprot_list.values():
    #     all_proteins = all_proteins.union(set(proteins))
    # print("%d unirefs and %d proteins"%(len(all_unirefs),len(all_proteins)))
    all_proteins = set(raw_caseolap_scores['entity'])
    print("%d proteins" % (len(all_proteins)))

    out_name = "./reactome_results/all_proteins_pathway_analysis_result.csv"
    pwa_df = submit_reactome_pathway_analysis(all_proteins, out_file=out_name)
    num_pathways = pwa_df.shape[0]
    print("%d pathways identified. Saved to %s" % (num_pathways, out_name))

    # proteins that had a score in each CVD at ALL
    nonzero_cvd_uniref = {}
    nonzero_cvd_proteins = {}
    for cvd in cvds:

        # extract proteins
        proteins = list(summary_table[summary_table[cvd] > 0].index)
        #     unirefs = list(summary_table[summary_table[cvd] > 0].index)
        #     proteins = convert_uniref_ids_to_proteins(unirefs,uniref_to_uniprot_list)
        nonzero_cvd_uniref[cvd] = unirefs
        nonzero_cvd_proteins[cvd] = proteins
        print("%s: %d proteins" % (cvd, len(proteins)))

        # pathway analysis
        if len(proteins) > 0:
            out_name = "./reactome_results/%s_nonzero_proteins_pathway_analysis_result.csv" % cvd
            pwa_df = submit_reactome_pathway_analysis(proteins, out_file=out_name)
            num_pathways = pwa_df.shape[0]
            print("%d pathways identified. Saved to %s" % (num_pathways, out_name))
        else:
            print("Insufficient number of proteins for %s" % cvd)


    load_reactome_data() #TODO

    # read all pathway analysis results into dataframes

    pathway_analysis_directory = './reactome_results'
    title_to_df = {}
    for filename in os.listdir(pathway_analysis_directory):
        if filename.endswith(".csv"):
            pwa_file = os.path.join(pathway_analysis_directory, filename)

            # rename some keys for formatting
            tag = filename.split("_")[0]  # cvd abbreviation
            if "category_ii" in filename:
                tag += "_unique"
            elif "subproteome" in filename:
                tag += '_subproteome'
            elif "nonzero" in filename:
                tag += "_nonzero"
            elif "category_i" in filename:
                tag = "category_i"
            elif "all" in filename:
                tag = "all_proteins"

            title_to_df[tag] = pd.read_csv(pwa_file)

        else:
            continue
    reactome_results = title_to_df['category_i']
    print(title_to_df.keys())

    # make heatmap using these values
    heatmap_data, pathways_in_heatmap = extract_heatmap_data(title_to_df)
    heatmap_data = heatmap_data[cvds]
    heatmap_data

    # the results of putting all the proteins into reactome pathway analysis
    all_union_pa_df = title_to_df['all_proteins']

    # get reverse and forward mapping of uniprot proteins to uniref90 groups
    # uniref_to_uniprot_list = get_uniref_to_uniprot_list(prot_to_uniref_file)
    # uniprot_to_uniref = get_uniref_to_uniprot_list(prot_to_uniref_file, reverse_mapping=True)
    uniprot_to_uniprot = {p: p for p in all_proteins}

    # extract pathway to its set of corresponding CaseOLAP scores for each protein in that pathway
    # pathway_to_uniref_scores = extract_pathway_to_scores(pathways_in_heatmap, zscores_df,
    #                                                     all_union_pa_df, uniprot_to_uniref, debug=True)
    pathway_to_uniref_scores = extract_pathway_to_scores(pathways_in_heatmap, zscores_df,
                                                         all_union_pa_df, uniprot_to_uniprot, debug=True)

    c_data = heatmap_data[cvds]
    z_data = extract_pathway_avg_zscore_matrix(pathway_to_uniref_scores)
    # relabel rows and reorder columns
    z_data.index = [pathway_id_to_pathway_name[p] for p in z_data.index]
    z_data = z_data[cvds]
    # dend_data = linkage_matrix

    make_heatmap(reordered_z_data, linkage_matrix, v_lim=(None, None))



    ### Heatmap unique to CVDs

    unique_to_cvd_pathway_list = {}
    for title, reactome_results_df in title_to_df.items():
        if "unique" in title:
            pathways = list(reactome_results_df['Pathway identifier'])
            temp_df = reactome_results_df.sort_values('Entities pValue',
                                                      ascending=True)  # make sure lowest pvalue is first
            cvd = title.split("_")[0]
            unique_to_cvd_pathway_list[cvd] = pathways[:3]
    for cvd, p in unique_to_cvd_pathway_list.items():
        print("%s: %s" % (cvd, p))
        pathway_names = []
        for pp in p:
            pathway_names += [pathway_id_to_pathway_name[pp]]
        print(pathway_names)

    make_heatmap_unique_to_cvd(unique_to_cvd_pathway_list, zscores_df,
                               all_union_pa_df, uniprot_to_uniprot, cvds,
                               pathway_id_to_pathway_name, reactome_pathway_to_unique_proteins)

    # make_heatmap_unique_to_cvd(unique_to_cvd_pathway_list, zscores_df,
    #                                                      all_union_pa_df, uniprot_to_uniref, cvds,
    #                                              pathway_id_to_pathway_name, reactome_pathway_to_unique_proteins)

def load_reactome_data():
    # Load in the reactome hierarchial information for human. See function definition in first cell.
    hierarchical_relationships, unique_reactome_ids = \
        extract_human_hierarchical_information(reactome_hierarchy_to_pathway)

    # Load in pathway id to proteins relationships. See function definition in first cell.
    reactome_pathway_to_unique_proteins = extract_pathway_to_proteins(reactome_uniprot_to_pathway,
                                                                      exclude_isoforms=False)

    # Load in mapping of reactome pathway id to pathway name
    pathway_id_to_pathway_name = extract_pathway_id_to_pathway_name(reactome_id_to_pathway_name)

    # Create the graph
    G = nx.DiGraph()
    G.name = 'Reactome hierarchical tree'
    G.add_edges_from(hierarchical_relationships)
    print(nx.info(G))

    # Calculate the familiar relations of each node
    # That is, label a node based on its membership as a descendent of a specific set of pathways (root and parent pathways)

    root_nodes = [x for x in G.nodes() if G.out_degree(x) > 0 and G.in_degree(x) == 0]
    print("Number of root nodes: " + str(len(root_nodes)))

    print("### Root Pathway Descendants ###")
    descendants = get_descendant_pathways(root_nodes, G)
    root_pathway_familiar_relations = get_familiar_relations(descendants)

    # Prune graph based on membership to root pathway
    pruned_G_root_pathways = prune_extra_edges(G,
                                               labeling_function=root_pathway_familiar_relations,
                                               debug=False)
    pruned_G_root_pathways.name = 'Reactome hierarchical tree pruned with root pathway membership'
    print(nx.info(pruned_G_root_pathways))

def generate_category_table(zscore_caseolap_file, z_score_threshold=3.0):

    # reading the csv file
    summary_table = pd.read_csv(zscore_caseolap_file)
    summary_table = summary_table.set_index('entity')
    #TODO zscores_df = subset of summary_table or copy

    # obtaining list of CVDs
    CVDs = list(summary_table.columns)
    print(CVDs)
    summary_table.head()

    # Sub-proteomes: Proteins with z-score above the threshold
    sub_proteome_headers = []
    for cvd in CVDs:
        cvd_sub = cvd + " Sub-proteome"
        summary_table[cvd_sub] = summary_table[cvd] >= z_score_threshold
        sub_proteome_headers += [cvd_sub]

    summary_table[sub_proteome_headers].sort_values(sub_proteome_headers[0], ascending=False).head()


    # extract the z-scores of protein groups in the subproteome for each CVD
    cvd_subproteome_to_scores = {}
    for cvd in CVDs:
        # extract protein groups in subproteome
        tag = cvd + ' Sub-proteome'
        proteins_in_subproteome = set(summary_table[summary_table[tag]].index)
        #     proteins_in_subproteome = uniref_list_to_protein_list(unirefs_in_subproteome,uniref_to_uniprot_list)
        print("%s %d proteins" % (tag, len(proteins_in_subproteome)))

        # extract z-scores
        scores = list(summary_table[summary_table.index.isin(proteins_in_subproteome)][cvd])
        cvd_subproteome_to_scores[cvd] = scores

    # prepare plot values above threshold line
    subproteome_plot_values = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in cvd_subproteome_to_scores.items()]))

    cvd_proteomes_violin_plot(zscores_df, subproteome_plot_values)

    # Category I: Protein groups with score > 0 for all CVDs
    summary_table['Category I'] = (summary_table.isnull().sum(axis=1) == 0)
    temp = CVDs + ['Category I']
    summary_table[summary_table['Category I'] == True][temp].head()

    # print out summary statistics
    category_i_proteins = set(summary_table[summary_table['Category I'] == True].index)
    # category_i_unirefs = set(summary_table[summary_table['Category I']==True].index)
    # category_i_proteins = uniref_list_to_protein_list(category_i_unirefs,uniref_to_uniprot_list)

    # print("%d unique protein groups and %d proteins in Category I"%(len(category_i_unirefs), len(category_i_proteins)))
    print("%d proteins in Category I" % (len(category_i_proteins)))

    make_category_i_heatmap(summary_table[summary_table['Category I'] == True][CVDs], sort_by='mean')

    zscore_cutoff_table(zscores_df)
    print("z-score threshold used for this analysis: %f" % (z_score_threshold))

    # Category II: explaionaation here

    category_ii_headers = []
    for cvd in CVDs:
        cvd_sub = cvd + " Sub-proteome"
        cvd_catII = cvd + " Category II"

        # gets the cvds that are not the ones of current interest
        cvd_other = [other for other in CVDs if (other != cvd)]

        # accesses the sub-proteome columns of the "other" cvds
        cvd_other = [x + " Sub-proteome" for x in cvd_other]

        # the new category column is if the sub-proteome of interest is >= 3 and the other sub-proteomes are <= 3
        summary_table[cvd_catII] = summary_table[cvd_sub] & ~(
            (summary_table[cvd_other]).any(bool_only=True, axis='columns'))

        category_ii_headers += [cvd_catII]

    # summary_table[category_ii_headers].head()

    # print out summary statistics
    for col_name in summary_table.columns:
        if col_name in category_ii_headers:
            num_true = summary_table[col_name].sum()
            print("%s\t%d" % (col_name, num_true))

    make_category_ii_heatmap(summary_table)

    combs = get_combinations(CVDs)

    # calculate new category III columns
    category_iii_headers = []
    for combination in combs:
        # gets the cvds that are not the ones of current interest
        cvd_other = [other for other in CVDs if (other not in list(combination))]
        # print(combination, cvd_other)

        # accesses the sub-proteome columns of the "other" cvds
        cvd_other = [x + " Sub-proteome" for x in cvd_other]

        # gets the cvds that are the ones of current interest
        # cvd_interest = [interest for interest in CVDs if (interest == cvd)]

        # accesses the sub-proteome columns of the cvds of interest
        cvd_interest = [x + " Sub-proteome" for x in combination]

        # the new category column is if the sub-proteomes of interest are >= 3 and the other sub-proteomes are <= 3
        summary_table[combination] = ((summary_table[cvd_interest]).all(bool_only=True, axis='columns')) & ~(
            (summary_table[cvd_other]).any(bool_only=True, axis='columns'))

        category_iii_headers += [combination]

    summary_table[category_iii_headers].sort_values(category_iii_headers[0], ascending=False).head()

    # print out summary statistics
    counter = 0
    upset_plot_dict = {}
    for col_name in summary_table.columns[25:]:
        num_true = summary_table[col_name].sum()
        if num_true != 0:
            counter = counter + 1
            #         print("%s\t%d"%(col_name,num_true))
            upset_plot_dict[col_name] = num_true

    # print(counter)
    for l, c in upset_plot_dict.items():
        print(l, c)
    print(len(upset_plot_dict))

    multi_index_categories = get_bool_table(list(upset_plot_dict.keys()), ordering=CVDs)
    # multi_index_categories
    data_series = pd.Series(data=upset_plot_dict.values(), index=multi_index_categories)
    data_series

    upplot(data_series)
    plt.show()

    summary_table['Mapped Proteins'] = merged['Mapped Proteins']
    summary_table['Synonyms'] = merged['Synonyms']
    new_order = ['Mapped Proteins', 'Synonyms'] + list(summary_table.columns[:-2])
    summary_table = summary_table[new_order]
    summary_table.head()

    # write table
    summary_table.to_csv('./results/SupplementaryData2_Summary_Table.csv') #TODO


def analyze_results(root_directory, z_score_thresh=3.0, merge_proteins=True, use_core_proteins=True, debug=True):

    output_directory = os.path.join(root_directory,'output/analyze_results_output')
    figure_direcotry = os.path.join(output_directory,'figures')
    # make data folders if doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    if not os.path.exists(figure_direcotry):
        os.makedirs(figure_direcotry)

    # error checking and load files
    raw_caseolap_scores, reactome_data, id_to_synonym = load_files(root_directory, merge_proteins, use_core_proteins)
    df =  raw_caseolap_scores
    hierarchical_relationships, unique_reactome_ids, reactome_pathway_to_unique_proteins, pathway_id_to_pathway_name = reactome_data

    # merge based on shared synonyms
    if merge_proteins:

        merged = merge_redundant_ids(id_to_synonym, df)

        if debug:
            print("Number of merged proteins before merging: %d"%df.shape[0])
            print("Number of merged proteins after merging by synonyms: %d"%merged.shape[0])
            print("%f%% of original"%(100*merged.shape[0]/df.shape[0]))

        df = merged.copy(deep=True).reset_index().drop(["Synonyms"], axis=1)

    # print(df.shape)
    if debug:
        print("%d unique UniRef IDs associated with at least one CVD" % (df.shape[0]))
        print("%d protein group-disease relationships (including 0's)" % (df.shape[0] * df.shape[1]))
        print("%d nonzero protein group-disease relationships" % (np.count_nonzero(df, axis=None)))
        print("%d zero protein group-disease relationships" % (df.shape[0] * df.shape[1] - np.count_nonzero(df, axis=None)))

    # all cardiovascular disease types
    cvds = list(df.columns)[2:]

    # make a boolean (True/False) table from original data
    df_bool = df.copy()
    for col in cvds:
        df_bool[col] = df[col] > 0

    # get only rows with all true

    all_true_df_bool = df_bool[df_bool[cvds].T.all()]
    unirefs_in_all = set(all_true_df_bool['entity'])
    print("%d UniRef90 IDs were found in all %s CVDs" % (len(unirefs_in_all), len(cvds)))

    # print top scoring protein in each category
    for cvd in cvds:
        max_val = np.max(df[cvd])
        max_idx = df[cvd].idxmax()
        prot_group = df.iloc[max_idx]['entity']
        print(cvd, max_val, prot_group)

    # df is the original CaseOLAP scores matrix.
    zscores_df = convert_to_zscore(df, include_zeros=False, columns_to_ignore=['entity'])
    zscores_df = zscores_df.set_index('entity')
    zscores_df.head()
    zscores_df.to_csv(os.path.join(output_directory,"merged_caseolap_zscores.csv"), index=True)



