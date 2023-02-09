
import os, sys
from analysis._03_text_mining_association_analysis import *

def parse_ranked_ent_df(entity_count_folder):

    # prepare ranked_ent_df
    ranked_ent ={}
    for filename in os.listdir(entity_count_folder):
        f = os.path.join(entity_count_folder, filename)
        # checking if it is a file
        if os.path.isfile(f) and "json" in f:
    #         print(tag)
            tag = f.split("_")[-1].strip(".json")
            ent_dict = json.load(open(f,'r'))
    #         print(ent_dict)
            ranked_ent[tag] = ent_dict
    ranked_ent_df = pd.DataFrame(ranked_ent)

    print(ranked_ent_df)
    return ranked_ent_df


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


def get_required_files(root_directory,merge_proteins, use_core_proteins):
    '''
    Required files: caseolap.csv, 3 reactome files, ent2syn2count.json (some sort of it)
    :param root_directory:
    :return:
    '''
    data_folder = os.path.join(root_directory,'data')

    ### input files ###

    # caseolap_result_folder = os.path.join(root_directory,"output/text_mining_results/all_proteins/")#TODO
    caseolap_result_folder = os.path.join(root_directory,"output/all_proteins/")
    caseolap_results_file = os.path.join(caseolap_result_folder, "all_caseolap.csv")
    if use_core_proteins:
        # caseolap_result_folder = os.path.join(root_directory, "output/text_mining_results/core_proteins/")#
        caseolap_result_folder = os.path.join(root_directory, "output/core_proteins/")
        caseolap_results_file = os.path.join(caseolap_result_folder, "core_caseolap.csv")

    reactome_uniprot_to_pathway_file = os.path.join(data_folder,"Reactome/UniProt2Reactome_All_Levels.txt")
    reactome_hierarchy_to_pathway_file = os.path.join(data_folder,"Reactome/ReactomePathwaysRelation.txt")
    reactome_id_to_pathway_name_file = os.path.join(data_folder,"Reactome/ReactomePathways.txt")

    # synonym related files for redundancy removal
    entity_count_file = os.path.join(caseolap_result_folder,"ranked_proteins/ranked_caseolap_score")
    # entity_count_files = []
    # if merge_proteins:
    #     entity_count_folder = os.path.join(caseolap_result_folder,'ranked_proteins/ranked_caseolap_score')
    #
    #     for f in os.listdir(entity_count_folder):
    #         if 'dictionary_ranked_proteins' in f:
    #             entity_count_files += [os.path.join(entity_count_folder,f)]

    ret = {'caseolap_results':caseolap_results_file,
           'reactome_data':[reactome_uniprot_to_pathway_file,
                            reactome_hierarchy_to_pathway_file,
                            reactome_id_to_pathway_name_file],
           'entity_count_files':entity_count_file}

    # check files
    for k,v in ret.items():
        if type(v) == list:
            files_to_check = v
        else:
            files_to_check = [v]
        for f in files_to_check:
            if not os.path.exists(f):
            # if (not os.path.isfile(f):
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


    hierarchical_relationships, unique_reactome_ids, reactome_pathway_to_unique_proteins, pathway_id_to_pathway_name, pruned_G_root_pathways,G, root_pathway_familiar_relations = load_reactome_data(input_files['reactome_data'])
    reactome_data = [hierarchical_relationships, unique_reactome_ids,reactome_pathway_to_unique_proteins, pathway_id_to_pathway_name,pruned_G_root_pathways,root_pathway_familiar_relations]

    # load id_to_synonym
    id_to_synonym = prepare_synonyms(parse_ranked_ent_df(input_files['entity_count_files']))

    return raw_caseolap_scores, reactome_data, id_to_synonym


def run_pathway_analysis(summary_table, cvds, 
                        hierarchical_relationships, unique_reactome_ids, reactome_pathway_to_unique_proteins, pathway_id_to_pathway_name, G, root_pathway_familiar_relations,
                        output_directory='.'):
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

        #cvd_subproteome_unirefs[cvd] = unirefs
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

        #category_ii_unirefs[cvd] = unirefs
        category_ii_proteins[cvd] = proteins
    
    reactome_results_dir = os.path.join(output_directory,"reactome_results")
    # make reactome result folder. TODO move to analyze_results function
    if not os.path.exists(reactome_results_dir):
        os.makedirs(reactome_results_dir)
    for cvd in cvds:
        proteins = cvd_subproteome_proteins[cvd]

        if len(proteins) > 0:
            out_name = os.path.join(reactome_results_dir, "%s_subproteome_pathway_analysis_result.csv" % cvd)
            pwa_df = submit_reactome_pathway_analysis(proteins, out_file=out_name)
            num_pathways = pwa_df.shape[0]
            print("%d pathways identified. Saved to %s" % (num_pathways, out_name))
        else:
            print("Insufficient number of proteins for %s" % cvd)

    # run pathway analysis for category i proteins
    cat_i_outfile = os.path.join(reactome_results_dir,"category_i_proteins_pathway_analysis_result.csv")
    cat_i_pa_result_df = submit_reactome_pathway_analysis(category_i_proteins, debug=True,
                                                          out_file=cat_i_outfile)
    cat_i_pa_result_df.head()

    category_ii_pa_df = {}
    for cvd in cvds:
        proteins = category_ii_proteins[cvd]

        if len(proteins) > 0:
            out_name = os.path.join(reactome_results_dir,"%s_category_ii_pathway_analysis_result.csv" % cvd)
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
    all_proteins = set(summary_table.index)
    print("%d proteins" % (len(all_proteins)))
    
    out_name = os.path.join(reactome_results_dir,"all_proteins_pathway_analysis_result.csv")
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
        #nonzero_cvd_uniref[cvd] = unirefs
        nonzero_cvd_proteins[cvd] = proteins
        print("%s: %d proteins" % (cvd, len(proteins)))

        # pathway analysis
        if len(proteins) > 0:
            out_name = os.path.join(reactome_results_dir,"%s_nonzero_proteins_pathway_analysis_result.csv" % cvd)
            pwa_df = submit_reactome_pathway_analysis(proteins, out_file=out_name)
            num_pathways = pwa_df.shape[0]
            print("%d pathways identified. Saved to %s" % (num_pathways, out_name))
        else:
            print("Insufficient number of proteins for %s" % cvd)

    # read all pathway analysis results into dataframes

    pathway_analysis_directory = reactome_results_dir #TODO fix
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

    zscores_df = summary_table[cvds]


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
    linkage_matrix, dendrogram_order = construct_dendrogram(G,pathways_in_heatmap, labeling_function=root_pathway_familiar_relations)
    heatmap_outfile = os.path.join(output_directory,"output/figures/CVD_Reactome_coverage_heatmap_no_dendrogram.pdf")
    reordered_c_data = reorder_heatmap(c_data,dendrogram_order,pathway_id_to_pathway_name,reactome_pathway_to_unique_proteins)
    reordered_z_data = reorder_heatmap(z_data,dendrogram_order,pathway_id_to_pathway_name,reactome_pathway_to_unique_proteins)
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

    pathway_unique_to_cvd_heatmap_outfile = os.path.join(output_directory,"output/figures/pathways_unique_to_cvd_heatmap.pdf")
    make_heatmap_unique_to_cvd(unique_to_cvd_pathway_list, zscores_df,
                               all_union_pa_df, uniprot_to_uniprot, cvds,
                               pathway_id_to_pathway_name, reactome_pathway_to_unique_proteins)

    # make_heatmap_unique_to_cvd(unique_to_cvd_pathway_list, zscores_df,
    #                                                      all_union_pa_df, uniprot_to_uniref, cvds,
    #                                              pathway_id_to_pathway_name, reactome_pathway_to_unique_proteins)


def load_reactome_data(input_files):
    reactome_uniprot_to_pathway, reactome_hierarchy_to_pathway, reactome_id_to_pathway_name = input_files
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
   
    return hierarchical_relationships, unique_reactome_ids, reactome_pathway_to_unique_proteins, pathway_id_to_pathway_name, pruned_G_root_pathways,G,root_pathway_familiar_relations


def generate_category_table(zscore_caseolap_file, merged, output_directory=".", z_score_threshold=3.0):

    # reading the csv file
    summary_table = pd.read_csv(zscore_caseolap_file)
    summary_table = summary_table.set_index('entity')
    summary_table = summary_table.drop('index',axis=1)
    zscores_df = summary_table.copy(deep=True)
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

    violin_plot_output_file = os.path.join(output_directory,"figures/caseolap_zscore_violin_plot.pdf")
    cvd_proteomes_violin_plot(zscores_df, subproteome_plot_values, show_figure=False, out_file=violin_plot_output_file)

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
    cat_i_heatmap_outfile = os.path.join(output_directory,'figures/Category_I_Heatmap.pdf')

    make_category_i_heatmap(summary_table[summary_table['Category I'] == True][CVDs], sort_by='mean',
                            out_file = cat_i_heatmap_outfile)

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

    category_ii_heatmap_outfile = os.path.join(output_directory,"figures/unique_unirefs_heatmap.pdf")
    make_category_ii_heatmap(summary_table, out_file=category_ii_heatmap_outfile)

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
#TODO make opplot output to table
    upplot(data_series)
#    plt.show()

    summary_table['Synonyms'] = list(merged['Synonyms'])
    new_order = ['Synonyms'] + list(summary_table.columns[:-1])
    summary_table = summary_table[new_order]
    summary_table.head()

    # write table
    out_file = os.path.join(output_directory,"SupplementaryData2_Summary_Table.csv")
    summary_table.to_csv(out_file) #TODO
    return summary_table

def analyze_results(root_directory, z_score_thresh=3.0, merge_proteins=True, use_core_proteins=True, debug=True):

    output_directory = os.path.join(root_directory,'output/analyze_results_output')
    figure_directory = os.path.join(output_directory,'figures')
    # make data folders if doesn't exist
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    if not os.path.exists(figure_directory):
        os.makedirs(figure_directory)

    # error checking and load files
    raw_caseolap_scores, reactome_data, id_to_synonym = load_files(root_directory, merge_proteins, use_core_proteins)
    df =  raw_caseolap_scores
    hierarchical_relationships, unique_reactome_ids, reactome_pathway_to_unique_proteins, pathway_id_to_pathway_name, G, root_pathway_familiar_relations = reactome_data
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
    zscore_table_outfile = os.path.join(output_directory,"merged_caseolap_zscores.csv")
    zscores_df.to_csv(zscore_table_outfile, index=True)

    summary_table = generate_category_table(zscore_table_outfile,merged,output_directory=output_directory,z_score_threshold=z_score_thresh)
    run_pathway_analysis(summary_table,cvds, hierarchical_relationships, unique_reactome_ids, reactome_pathway_to_unique_proteins, pathway_id_to_pathway_name, G, root_pathway_familiar_relations, output_directory=output_directory)
