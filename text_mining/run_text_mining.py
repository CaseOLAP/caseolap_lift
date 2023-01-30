
import json, sys, os, time, traceback, sys
from collections import defaultdict
import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns
from elasticsearch_dsl import Search, Q
from multiprocessing import cpu_count, Process
from elasticsearch import Elasticsearch

# setting path
sys.path.append('..')
from text_mining.caseolap._01_download import *
from text_mining.caseolap._02_parsing import *
from text_mining.caseolap._03_mesh2pmid import *
from text_mining.caseolap._05_index_populate import *
from text_mining.caseolap._06_textcube import *
from text_mining.caseolap._07_vary_synonyms_cases import *
from text_mining.caseolap._08_count_synonyms import *
from text_mining.caseolap._09_screen_synonyms import *
from text_mining.caseolap._10_make_entity_counts import *
from text_mining.caseolap._11_metadata_update import *
from text_mining.caseolap._12_caseolap_score import *
from text_mining.caseolap._13_inspect_entity_scores import *


def run_text_mining(root_dir, data_folder, mapping_folder, analysis_output_folder,
                    date_range=None,
                    include_full_text=False,
                    include_label_imputation=False,
                    check_synonyms=False,
                    rerun_scoring=False):

    # directories
    # data_dir = os.path.join(root_dir,'data')
    # result_dir = os.path.join(root_dir,'result') # Main folder where the results from this section will be stored
    log_dir = os.path.join(root_dir,'log')
    config_dir = os.path.join(root_dir,'config')
    input_dir = os.path.join(root_dir,'input') #TODO how is input folder different?


    # Input 01
    logFilePath = os.path.join(log_dir,'download_log.txt')
    baseline_dir = os.path.join(data_folder, 'ftp.ncbi.nlm.nih.gov/pubmed/baseline/') # Documents prior to this year
    update_files_dir = os.path.join(data_folder,'ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/') # Documents prior to this year

    download_config_file_path = os.path.join(config_dir,'download_config.json')
    ftp_config_file_path = os.path.join(config_dir,'ftp_config.json')

    # Input 02
    parsing_config_file = os.path.join(config_dir,'parsing_config.json')

    # Output 02
    pubmed_path = os.path.join(data_folder,'pubmed.json')     # The parsed PubMed documents (dictionary)
    filestat_path = os.path.join(data_folder,'filestat.json') # Unzipped PubMed download file name : #PMIDs
    logfile_path = os.path.join(log_dir,'parsing_log.txt') # Logs progress on parsing


    ### Input 03

    ### Output 03
    mesh2pmid_outputfile = os.path.join(data_folder,"mesh2pmid.json" )   # {"MeSH Term":[PMID1,...], ...}
    mesh2pmid_statfile = os.path.join(data_folder,"mesh2pmid_stat.json") # {"MeSH Term": #PMIDs, ...}
    logFilePath = os.path.join(log_dir,'mesh2pmid_log.txt')           # Logs mapping progress


    # Index parameters 04
    index_name = 'pubmed_lift'      # Index name (match with 05 index populate file)
    type_name = 'pubmed_meta_lift'  # Index type name (match with 05 index populate file)
    number_shards = 1          # Set to 1 if no cluster
    number_replicas = 0
    case_sensitive = True   # Index the text as case sensitive (True) or lower case (False)

    # Input file 04
    index_init_config_file = os.path.join(config_dir,'index_init_config.json')


    # Input 05
    index_populate_config_file = os.path.join(config_dir,'index_populate_config.json')

    # Output 05
    logfile_path = os.path.join(log_dir,'indexing_log.txt')            # Reports progress on indexing

    # Names of the index you want to create 05

    # Input data 06
    meshtree = os.path.join(data_folder,'MeSH/mtrees2021.bin')                 # MeSH Tree ontology
    root_cat = os.path.join(data_folder,'categories.txt')                 # Categories' root MeSH Tree nums
    textcube_config = os.path.join(config_dir,'textcube_config.json')   # ['CategoryName1',...]
    mesh2pmid = mesh2pmid_outputfile # from step 2

    # Output data directories 06
    textcube_pmid2category = os.path.join(data_folder,'textcube_pmid2category.json') # Map PMID to category
    textcube_category2pmid = os.path.join(data_folder,'textcube_category2pmid.json') # Map category to PMID
    textcube_stat = os.path.join(data_folder,'textcube_stat.txt')                    # Num. documents per category
    MeSHterms_percat = os.path.join(data_folder,'meshterms_per_cat.json')            # MeSH *Terms* per category
    logfile_path = os.path.join(log_dir,'textcube_log.txt')                       # Logs file's messages


    # Input 07
    #TODO move this one folder above, not in TFD
    entity_dict_path_no_cs = os.path.join(mapping_folder,'Transcription_Factor_Dependence/id2synonyms_not_case_varied.json') # entity dict
    species = ['human', 'pig', 'mouse', 'rat']  # Species studied. Used for permuting syns.

    # Output 07
    case_varied_entites_outpath = os.path.join(data_folder,'casesensitive_entities.txt') # Case sensitive entity dict
    case_varied_entity_dict_path = os.path.join(input_dir,'id2syns.json') # entity dict
    core_id2syns=os.path.join(input_dir,'core_id2syns.json')
    core_proteins_file=os.path.join(root_dir,'output/core_proteins.txt')

    # Input 08
    entity_dict_path = case_varied_entity_dict_path
    #entity_dict_path = 'input/id2syns.json'
    textcube_pmid2category = os.path.join(data_folder,'textcube_pmid2category.json')
    if date_range != None:
        start_year = int(date_range[0].split("-")[0])
        end_year = int(date_range[1].split("-")[0])
    else:
        start_year = 1800
        end_year = 2022
        # TODO
    # start_year = 2012    # <-- Change as you see fit for your publications of interest
    # end_year = 2022      # <-- Same as above comment

    # Intermediary file (produced as    c
    syn_pmid_count = os.path.join(data_folder,'syn_pmid_count.txt')

    # Output 08
    pmid_syn_count_out = os.path.join(data_folder,'pmid_synonym_counts_2012-2022.json')   # PMID Syn|Count...Syn|Cnt
    synfound_pmid2cat = os.path.join(data_folder,'synfound_pmid2category_2012-2022.txt')  # PMID--->CategoryNumber
    logfile = os.path.join(log_dir,'synonymcount_log_2012-2022.txt')                   # #hits:Synonym


    # Input 09
    id2syns = os.path.join(input_dir,'id2syns.json')

    # Output 09
    eng_path = os.path.join(data_folder,'suspect_synonyms_english_words.txt')  # Ambiguous syns: English word
    short_path = os.path.join(data_folder,'suspect_synonyms_too_short.txt')  # Ambiguous syns: Short words
    rem_path = os.path.join(data_folder,'remove_these_synonyms.txt')         # Synonyms to remove


    # Other parameters 09
    index_name = 'pubmed_lift' # Index name
    key = 'abstract'        # Choose if searching the abstracts and titles
    key = 'full_text'      # Choose if searchines the abstracts, titles, and full text
    key = 'full_text_no_methods' # Same as above, excluding abstracts



    # Input 10
    remove_syns_infile = os.path.join(data_folder,'remove_these_synonyms.txt')
    all_id2syns_path = os.path.join(input_dir,'id2syns.json')
    core_id2syns_path = os.path.join(input_dir,'core_id2syns.json')
    pmid_syn_count_in = os.path.join(data_folder,'pmid_synonym_counts_2012-2022.json')

    # Output 10
    all_entitycount_outfile = os.path.join(data_folder,'all_entitycount_2012-2022.txt')
    core_entitycount_outfile = os.path.join(data_folder,'core_entitycount_2012-2022.txt')



    # Input file paths 11
    all_entitycount_path = os.path.join(data_folder,'all_entitycount_2012-2022.txt')              # PMID Entity|Count ...
    core_entitycount_path = os.path.join(data_folder,'core_entitycount_2012-2022.txt')              # PMID Entity|Count ...
    pmid2category_path = os.path.join(data_folder,'textcube_pmid2category.json')# PMIDs of interest to category
    category_names_file = os.path.join(config_dir,'textcube_config.json')  # Category names


    # Output file paths 11
    all_outfile_pmid2entity2count = os.path.join(data_folder,'all_metadata_pmid2entity2count_2012-2022-2.json') # {PMID:{Entity:Count,...},...}
    core_outfile_pmid2entity2count = os.path.join(data_folder,'core_metadata_pmid2entity2count_2012-2022-2.json') # {PMID:{Entity:Count,...},...}
    all_cat2pmids_path = os.path.join(data_folder,'all_metadata_category2pmids_2012-2022-2.json')  # {CatName:[PMID,...], ...}
    core_cat2pmids_path = os.path.join(data_folder,'core_metadata_category2pmids_2012-2022-2.json')  # {CatName:[PMID,...], ...}
    all_logfile_path = os.path.join(log_dir,'all_metadata_update_log_2012-2022.txt')         # Similar to pmid2pcount
    core_logfile_path = os.path.join(log_dir,'core_metadata_update_log_2012-2022.txt')         # Similar to pmid2pcount

    # Input data directories 12
    all_cat2pmids_path =  os.path.join(data_folder,'all_metadata_category2pmids_2012-2022-2.json') # {CatName:[PMID,...], ...}
    core_cat2pmids_path = os.path.join(data_folder,'core_metadata_category2pmids_2012-2022-2.json')  # {CatName:[PMID,...], ...}
    all_pmid2entity2count_path = os.path.join(data_folder,'all_metadata_pmid2entity2count_2012-2022-2.json') # {PMID:{Entity:Count,...},...}
    core_pmid2entity2count_path = os.path.join(data_folder,'core_metadata_pmid2entity2count_2012-2022-2.json') # {PMID:{Entity:Count,...},...}
    category_names_path = os.path.join(config_dir,'textcube_config.json')    # ['CategoryName1',...]
    # config_dir
    # Output data path 12
    all_logFilePath = os.path.join(log_dir,'all_caseolap_score_log.txt') # Logs #PMIDs for each category
    core_logFilePath = os.path.join(log_dir,'core_caseolap_score_log.txt') # Logs #PMIDs for each category
    all_caseolap_name = 'all_caseolap'  # Name of dataframe/spreadsheet for the caseolap scores
    core_caseolap_name = 'core_caseolap'

    # Input path 13
    id2syns_in = os.path.join(input_dir,'id2syns.json')                # The case-varied entity dict
    pmid_syn_count_in = os.path.join(data_folder,'pmid_synonym_counts_2012-2022.json')  # Counts of each synonym
    remove_syns_in = os.path.join(data_folder,'remove_these_synonyms.txt')    # Syns that were not used

    all_caseolap_scores_in = os.path.join(analysis_output_folder,'all_proteins/all_caseolap.csv')  # The CaseOLAP scores
    all_popular_scores_in =  os.path.join(analysis_output_folder,'all_proteins/all_popularity_score.csv')
    all_distinct_scores_in =  os.path.join(analysis_output_folder,'all_proteins/all_distinctiveness_score.csv')
    all_cat2pmids_in = os.path.join(data_folder,'all_metadata_category2pmids_2012-2022-2.json')  # {CatName:[PMID,...], ...}

    core_caseolap_scores_in =  os.path.join(analysis_output_folder,'core_proteins/core_caseolap.csv') # The CaseOLAP scores
    core_popular_scores_in =  os.path.join(analysis_output_folder,'core_proteins/core_popularity_score.csv')
    core_distinct_scores_in =  os.path.join(analysis_output_folder,'core_proteins/core_distinctiveness_score.csv')
    core_cat2pmids_in = os.path.join(data_folder,'core_metadata_category2pmids_2012-2022-2.json')  # {CatName:[PMID,...], ...}

    # Output paths 13
    all_ranked_syns_out = os.path.join(analysis_output_folder,'all_proteins/ranked_synonyms/ranked_synonyms.txt')  # Syns ranked by counts
    all_direct = os.path.join(analysis_output_folder,'all_proteins/ranked_proteins/')
    all_ranked_caseolap_out = all_direct + 'ranked_caseolap_score/ranked_proteins.txt'  # Proteins ranked by score
    all_ranked_popular_out = all_direct + 'ranked_popularity_score/ranked_proteins.txt'
    all_ranked_distinct_out = all_direct + 'ranked_distinctiveness_score/ranked_proteins.txt'

    core_ranked_syns_out = os.path.join(analysis_output_folder,'core_proteins/ranked_synonyms/ranked_synonyms.txt') # Syns ranked by counts
    core_direct = os.path.join(analysis_output_folder,'core_proteins/ranked_proteins/')
    core_ranked_caseolap_out = core_direct + 'ranked_caseolap_score/ranked_proteins.txt'  # Proteins ranked by score
    core_ranked_popular_out = core_direct + 'ranked_popularity_score/ranked_proteins.txt'
    core_ranked_distinct_out = core_direct + 'ranked_distinctiveness_score/ranked_proteins.txt'


    # make directories
    for dir in [analysis_output_folder, log_dir,input_dir,baseline_dir,update_files_dir]:
        isExist = os.path.exists(dir)
        if not isExist:
            # Create a new directory because it does not exist
            os.makedirs(dir)
#    print("01_run_download")
#    text_mining_01_run_download(data_folder, logFilePath, download_config_file_path, ftp_config_file_path,
#                               baseline_dir,update_files_dir)
#
#    print("02_run_parsing")
#    text_mining_02_run_parsing(baseline_dir, update_files_dir,
#                              parsing_config_file, pubmed_path, filestat_path,
#                              logfile_path)
#
#    print("03_run_mesh2pmid")
#    text_mining_03_run_mesh2pmid(pubmed_path,
#                                mesh2pmid_outputfile, mesh2pmid_statfile, logFilePath)
#    print("04_run_index_init")
#    text_mining_04_run_index_init(index_name, type_name, index_init_config_file,
#                                 number_shards=1, number_replicas=0, case_sensitive=True)
#    print("05_run_index_populate")
#    text_mining_05_run_index_populate(pubmed_path, index_populate_config_file,
#                                  logfile_path, index_name, type_name)
#    print("06_run_textcube")
#    text_mining_06_run_textcube(textcube_category2pmid, meshtree, mesh2pmid, root_cat,
#                               textcube_config, textcube_pmid2category,
#                               textcube_stat, MeSHterms_percat,
#                               logfile_path)
#    print("07_run_vary_synonyms_cases")
#    text_mining_07_run_vary_synonyms_cases(entity_dict_path_no_cs, species,
#                                          case_varied_entites_outpath,
#                                          case_varied_entity_dict_path,
#                                          core_proteins_file,
#                                          core_id2syns)
#    print("08_run_count_synonyms")
#    text_mining_08_run_count_synonyms(entity_dict_path, textcube_pmid2category,
#                                  start_year, end_year, syn_pmid_count, pmid_syn_count_out,
#                                  synfound_pmid2cat, logfile, index_name, key)
#    print("09_run_screen_synonyms")
#    text_mining_09_run_screen_synonyms(data_folder, id2syns, eng_path, short_path, rem_path)

#    print("10_run_make_entity_counts")
#    text_mining_10_run_make_entity_counts(remove_syns_infile, all_id2syns_path, core_id2syns_path, pmid_syn_count_in,
#                                        all_entitycount_outfile, core_entitycount_outfile)
#    print("11_run_metadata_update")
#    text_mining_11_run_metadata_update(all_entitycount_path, core_entitycount_path, pmid2category_path,
#                                      category_names_file, all_outfile_pmid2entity2count,
#                                      core_outfile_pmid2entity2count, all_cat2pmids_path, core_cat2pmids_path,
#                                      all_logfile_path, core_logfile_path)
#    print("12_run_caseolap_score")
#    text_mining_12_run_caseolap_score(all_cat2pmids_path, core_cat2pmids_path, all_pmid2entity2count_path,
#                                      core_pmid2entity2count_path, category_names_path, all_logFilePath,
#                                      core_logFilePath, all_caseolap_name, core_caseolap_name, analysis_output_folder)

    print("13_run_inspect_entity_scores")
    text_mining_13_run_inspect_entity_scores(id2syns_in, pmid_syn_count_in, remove_syns_in, all_caseolap_scores_in,
                                             all_popular_scores_in, all_distinct_scores_in, all_cat2pmids_in,
                                             core_caseolap_scores_in, core_popular_scores_in, core_distinct_scores_in,
                                             core_cat2pmids_in, all_ranked_syns_out, all_ranked_caseolap_out,
                                             all_ranked_popular_out, all_ranked_distinct_out, core_ranked_syns_out,
                                             core_ranked_caseolap_out, core_ranked_popular_out,
                                             core_ranked_distinct_out)
    print("Finished with Text Mining Module")
# pubmed_path!!pubmed_path = os.path.join(data_dir,"pubmed.json")
# pubmed_path!!pubmed_path = os.path.join(data_dir,'pubmed.json') #'data/pubmed.json' # Parsed publications
# baseline_dir!!baseline_dir = os.path.join(data_dir,'ftp.ncbi.nlm.nih.gov/pubmed/baseline' ) #
# update_files_dir!!update_files_dir = os.path.join(data_dir,'ftp.ncbi.nlm.nih.gov/pubmed/updatefiles') # Documents from this year
# remove the load!!index_populate_config = json.load(open(index_populate_config_file))
# mesh2pmid_outputfile!!mesh2pmid = os.path.join(data_dir,'mesh2pmid.json') # MeSH to PMID mapping



def text_mining_01_run_download(data_dir, logFilePath, download_config_file_path, ftp_config_file_path,
                               baseline_dir,update_files_dir):
    '''
    The purpose of this file is to download the zipped files containing
    the PubMed publications (i.e. documents). These will be mined later.
    '''
    # Start the download, verification, and extraction process

    # Open log file
    logfile = open(logFilePath, 'w') 
    
    # Load config files
    download_config = json.load(open(download_config_file_path))
    ftp_config = json.load(open(ftp_config_file_path))
    
    # Check main directory
    if not os.path.isdir(data_dir):
        print('Directory not found:', data_dir) 
    
    # Start download 
    download_pubmed(data_dir, download_config, ftp_config, logfile)
    
    # Verify download: 'baseline files' & 'update files'
    check_all_md5_in_dir(baseline_dir, logfile, linux = True, mac = False)
    check_all_md5_in_dir(update_files_dir, logfile, linux = True, mac = False)

    # Extract downloaded files: 'baseline files' & 'update files'
    extract_all_gz_in_dir(baseline_dir, logfile)
    extract_all_gz_in_dir(update_files_dir, logfile)
    
    logfile.close()


def text_mining_02_run_parsing(baseline_dir, update_files_dir,
                              parsing_config_file, pubmed_path, filestat_path,
                              logfile_path):
    '''
    The purpose of this file is to parse the downloaded PubMed documents,
    saving their information into a dictionary.
    '''
    # Start time
    t1 = time.time()

    # Open the files
    logfile = open(logfile_path, 'w')
    pubmed = open(pubmed_path, 'w')
    filestat = open(filestat_path, 'w')
    parsing_config = json.load(open(parsing_config_file, 'r'))
    print(parsing_config)

    # Parse the files (baseline and updatefiles)
    parse_dir(baseline_dir, pubmed, filestat, 'baseline', parsing_config, logfile)
    parse_dir(update_files_dir, pubmed, filestat, 'updatefiles', parsing_config, logfile)

    # Close the files
    pubmed.close()
    filestat.close()
    logfile.close()
    t2 = time.time()

    # Print info
    print('Parsing finished, results dumped to', pubmed_path)
    print('Total Time: ', str((t2 - t1) / 360), 'hours')


def text_mining_03_run_mesh2pmid(pubmed_path,
                                mesh2pmid_outputfile, mesh2pmid_statfile, logFilePath):
    '''
    The purpose of this file is to map MeSH Terms (metadata) to PMIDs
    '''
    # Open log file
    logfile = open(logFilePath, "w")

    # Initialize class
    mesh2PMID = MeSH2PMID()

    # Map MeSH to PMID
    mesh2PMID.mesh2pmid_mapping(pubmed_path, mesh2pmid_outputfile, logfile)

    # Map MeSH to #PMIDs
    mesh2PMID.mesh2pmid_mapping_stat(mesh2pmid_statfile, logfile)

    # Close log file
    logfile.close()


def text_mining_04_run_index_init(index_name, type_name, index_init_config_file,
                                 number_shards=1, number_replicas=0, case_sensitive=True):
    '''
    The purpose of this file is to initialize the ElasticSearch index
    which is where the indexed PubMed text data will be.
    This version allows the option to preserve case-sensitivity in the
    indexed text.
    '''
    # Load the indexing config file
    index_init_config = json.load(open(index_init_config_file, 'r'))

    # Start elasticsearch
    # es = Elasticsearch('https://localhost:9200')
    es = Elasticsearch()

    # Delete the old index if it exists
    if es.indices.exists(index=index_name):
        res = es.indices.delete(index=index_name)
        print('Deleted index:', index_name, '\nResponse:', res, '\n')

    # Request Body Parameters
    mappings = {type_name: {'properties': index_init_config}}
    settings = {'number_of_shards': number_shards, 'number_of_replicas': number_replicas}

    if case_sensitive == True:
        settings['analysis'] = {'analyzer': {'casesensitive_text': {'type': 'custom',
                                                                    'tokenizer': 'standard',
                                                                    'filter': ['stop']}}}
    # Create an index
    res = es.indices.create(index=index_name, settings=settings, mappings=mappings)
    print('Created index:', index_name, '\nResponse:', res)


def text_mining_05_run_index_populate(pubmed_path,index_populate_config_file,
                                     logfile_path,index_name, type_name):
    '''
    The purpose of this file is to populate the ElasticSearch index.
    Make sure ElasticSearch is running.
    '''
    # Open the log file
    logfile = open(logfile_path, 'w')
    index_populate_config = json.load(open(index_populate_config_file))
    # Populate the index
    populate_index(pubmed_path, logfile, index_name, type_name, index_populate_config)

    # Close the log file
    logfile.close()


def text_mining_06_run_textcube(textcube_category2pmid, meshtree, mesh2pmid, root_cat,
                               textcube_config, textcube_pmid2category,
                               textcube_stat, MeSHterms_percat,
                               logfile_path):
    '''
    The purpose of this file is to create mappings for the PMIDs to their categories.
    Categories (e.g., diseases) are identified via MeSH term metadata.
    The TextCube is mappings between PMID to category
    '''
    logfile = open(logfile_path, 'w')

    # Get category names
    category_names = json.load(open(textcube_config, 'r'))
    print(str(len(category_names)), 'categories: ', category_names)

    # Initialize TextCube class
    TC = TextCube(category_names)

    # Get categories' descendant MeSH terms
    TC.descendant_mesh(root_cat, meshtree, MeSHterms_percat, logfile)

    # Get categories' PMIDs
    TC.category2pmids_mapping(mesh2pmid, textcube_category2pmid, logfile)

    # Get PMIDs' categories
    TC.pmid2category_mapping(textcube_pmid2category, logfile)

    # Display the number of documents per category
    TC.category_statistics(textcube_stat, logfile)

    logfile.close()


def text_mining_07_run_vary_synonyms_cases(entity_dict_path, species,
                                          case_varied_entites_outpath,
                                          case_varied_entity_dict_path,
                                          core_proteins,
                                          core_id2syns_outfile):
    '''
    The purpose of this file is to create case-sensitive variations
    of the synonyms. The synonyms will then be queried.
    '''
    # Instantiates the class, loads entity dictionary mapping ID to synonyms
    VSC = VarySynonymsCases(entity_dict_path, case_varied_entites_outpath)
    # Adds some more synonyms
    VSC.add_species_syns(species)

    # Makes case-sensitive variations of the synonyms
    VSC.gets_final_syns()

    # Make a dictionary of the ID to synonym data
    make_id2syns_dict(case_sensitive_entities_file=case_varied_entites_outpath,
                id2syns_outfile=case_varied_entity_dict_path,
                core_proteins_file=core_proteins,
                core_id2syns_outfile=core_id2syns_outfile)


def text_mining_08_run_count_synonyms(entity_dict_path, textcube_pmid2category,
                                     start_year, end_year, syn_pmid_count, pmid_syn_count_out,
                                     synfound_pmid2cat, logfile, index_name, key):
    # Instantiate the object
    CS = CountSynonyms(entity_dict_path, textcube_pmid2category)
    # Search for the synonyms in the indexed text
    CS.synonym_search(key, logfile, syn_pmid_count, index_name, start_year, end_year)

    # Finalize the output files
    CS.finish_synonym_search(logfile, syn_pmid_count, \
                             pmid_syn_count_out, synfound_pmid2cat)


def text_mining_09_run_screen_synonyms(data_dir, id2syns, eng_path, short_path, rem_path):
    '''
    The purpose of this file is to screen and modify entity synonyms.
    - It screens the entity synonyms (checks for potentially ambiguous synonyms).
    - It saves the potentially ambiguous synonyms (rem_path). You will look at these
    synonyms and determine if they are too ambiguous for you.
    - It also may also add some synonyms (add_species_syns).

    Suspect = potentially bad = ambiguous
    '''
    # Instantiate the class (Load the entity dictionary, synonyms)
    SE = ScreenSynonyms(id2syns)

    # Identifies synonyms to possible remove
    SE.find_suspect_synonyms(data_dir=data_dir)

    # Merge the temp files of suspect synonyms
    SE.merge_suspect_syns(data_dir=data_dir)

    # Save the synonyms you might want to remove. Inspect and modify this file
    SE.export_suspect_synonyms(eng_path, short_path, rem_path)


def text_mining_10_run_make_entity_counts(remove_syns_infile, all_id2syns_path, core_id2syns_path, pmid_syn_count_in,
                                         all_entitycount_outfile, core_entitycount_outfile):
    '''
    The purpose of this class is to make the entitycount file that maps
    PMIDs to entities to entity counts. This uses a "PMIDs to synonyms to synonym counts"
    mapping and an "entity to synonym" mapping to do this. No ElasticSearch querying
    is needed in this step.
    '''
    #### Expanded set proteins
    # Instantiate and initialize class
    print("MEC for all entities: ",all_id2syns_path)
    MEC = MakeEntityCounts(all_id2syns_path, remove_syns_infile, pmid_syn_count_in)

    # Makes "pmid->entity->count" from "pmid->syn->count" & "entity->syn" mappings
    MEC.entitycount(all_entitycount_outfile)

    #### Core Proteins
    # Instantiate and initialize class
    print("MEC for core entities: ",core_id2syns_path)
    MEC = MakeEntityCounts(core_id2syns_path, remove_syns_infile, pmid_syn_count_in)

    # Makes "pmid->entity->count" from "pmid->syn->count" & "entity->syn" mappings
    MEC.entitycount(core_entitycount_outfile)


def text_mining_11_run_metadata_update(all_entitycount_path, core_entitycount_path, pmid2category_path,
                                      category_names_file, all_outfile_pmid2entity2count,
                                      core_outfile_pmid2entity2count, all_cat2pmids_path, core_cat2pmids_path,
                                      all_logfile_path, core_logfile_path):
    '''
    The purpose of this file is to format mappings that can be used for the CaseOLAP scoring.
    The mappings are the outputs:
    - pmid2pcount_path: {'PMID':{'Entity':'Hits',...},...}
    - category2pmids_path: {'Category1': ['PMID_1',...,'PMID_n' ], ...}

    Starting with the textcube provided PMIDs in each category, this produces mappings of
    PMIDs in each category only if those PMIDs were found to contain entities of interest.
    '''

    #### Expanded protein set
    # Open log file
    logfile = open(all_logfile_path, 'w')

    # Get category names
    category_names = json.load(open(category_names_file, 'r'))

    # Initialize class
    MU = MetadataUpdate(category_names)

    # Rewrite PMID->Entity->Entity Count as a nested dictionary
    MU.update_pmid2entity2count(all_entitycount_path, all_outfile_pmid2entity2count, logfile)

    # Category->PMID (PMIDs in which queried entities were discovered)
    MU.map_category2pmid_pmids_with_entities(pmid2category_path, all_cat2pmids_path, logfile)

    # Close log file
    logfile.close()

    ### Core protein set
    # Open log file
    logfile = open(core_logfile_path, 'w')

    # Get category names
    category_names = json.load(open(category_names_file, 'r'))

    # Initialize class
    MU = MetadataUpdate(category_names)

    # Rewrite PMID->Entity->Entity Count as a nested dictionary
    MU.update_pmid2entity2count(core_entitycount_path, core_outfile_pmid2entity2count, logfile)

    # Category->PMID (PMIDs in which queried entities were discovered)
    MU.map_category2pmid_pmids_with_entities(pmid2category_path, core_cat2pmids_path, logfile)

    # Close log file
    logfile.close()

def text_mining_12_run_caseolap_score(all_cat2pmids_path, core_cat2pmids_path, all_pmid2entity2count_path,
                                      core_pmid2entity2count_path, category_names_path, all_logFilePath,
                                      core_logFilePath, all_caseolap_name, core_caseolap_name, result_dir):

    '''
    The purpose of this file is to produce CaseOLAP scores for the entities
    based on their hits in each document (pmid2pcount_path) and the documents'
    category (category2pmids_path).
    '''

    #### All proteins ####
    print('======All Proteins======')
    logfile = open(all_logFilePath, 'w')
    category2pmids = json.load(open(all_cat2pmids_path, 'r'))
    pmid2entity2count = json.load(open(all_pmid2entity2count_path, 'r'))

    ''' Initial Calculations'''
    # Initialize object with input data
    C = Caseolap(category2pmids, pmid2entity2count, result_dir, category_names_path, logfile)

    # Print info on categories and their number of publications
    C.print_categories2pmid(all_or_core='all', dump=True, verbose=True)

    # Map Category to its PMIDs to its Entities to the Entity's Counts
    C.map_category2pmid2entity2count()

    # Save all entities
    C.get_all_entities('all', dump=True, verbose=True)

    ''' Popularity Score (note: relies on some previous sections above)'''
    # Get the entity counts for each category
    C.get_entity_counts_per_category()

    # Maps category to its entities to their counts (includes zero count entities)
    C.category2entity2tf_finder()

    # Calculate the popularity scores for all entities
    C.calculate_all_popularity_scores(all_or_core='all', dump=True)

    ''' Distinctiveness Score (note: relies on some previous sections above)'''
    # Map entities to the count of their PMIDs
    C.category2entity2num_pmids_finder()

    # Calculate normalized term frequencies
    C.calculate_category2entity2ntf()

    # Calculate normalized document frequencies
    C.calculate_category2entity2ndf()

    # Calculate ratio of normalized term frequency over normalized document frequency
    C.calculate_entity2ntf_ndf_ratio()

    # Calculate distinctiveness score
    C.calculate_all_distinctiveness_scores(all_or_core='all', dump=True)

    '''Final Score'''
    # Calculate CaseOLAP Score (combine popularity & distinctiveness)
    C.calculate_caseolap_score(all_or_core='all', dump=True)

    # Close logfile
    logfile.close()

    #### Core proteins #####
    print('\n======Core Proteins======')
    logfile = open(core_logFilePath, 'w')
    category2pmids = json.load(open(core_cat2pmids_path, 'r'))
    pmid2entity2count = json.load(open(core_pmid2entity2count_path, 'r'))

    ''' Initial Calculations'''
    # Initialize object with input data
    C = Caseolap(category2pmids, pmid2entity2count, result_dir, category_names_path, logfile)

    # Print info on categories and their number of publications
    C.print_categories2pmid(all_or_core='core', dump=True, verbose=True)

    # Map Category to its PMIDs to its Entities to the Entity's Counts
    C.map_category2pmid2entity2count()

    # Save all entities
    C.get_all_entities('core', dump=True, verbose=True)

    ''' Popularity Score (note: relies on some previous sections above)'''
    # Get the entity counts for each category
    C.get_entity_counts_per_category()

    # Maps category to its entities to their counts (includes zero count entities)
    C.category2entity2tf_finder()

    # Calculate the popularity scores for all entities
    C.calculate_all_popularity_scores(all_or_core='core', dump=True)

    ''' Distinctiveness Score (note: relies on some previous sections above)'''
    # Map entities to the count of their PMIDs
    C.category2entity2num_pmids_finder()

    # Calculate normalized term frequencies
    C.calculate_category2entity2ntf()

    # Calculate normalized document frequencies
    C.calculate_category2entity2ndf()

    # Calculate ratio of normalized term frequency over normalized document frequency
    C.calculate_entity2ntf_ndf_ratio()

    # Calculate distinctiveness score
    C.calculate_all_distinctiveness_scores(all_or_core='core', dump=True)

    '''Final Score'''
    # Calculate CaseOLAP Score (combine popularity & distinctiveness)
    C.calculate_caseolap_score(all_or_core='core', dump=True)

    # Close logfile
    logfile.close()




def text_mining_13_run_inspect_entity_scores(id2syns_in, pmid_syn_count_in, remove_syns_in, all_caseolap_scores_in,
                                        all_popular_scores_in, all_distinct_scores_in, all_cat2pmids_in,
                                        core_caseolap_scores_in, core_popular_scores_in, core_distinct_scores_in,
                                        core_cat2pmids_in, all_ranked_syns_out, all_ranked_caseolap_out,
                                        all_ranked_popular_out, all_ranked_distinct_out, core_ranked_syns_out,
                                        core_ranked_caseolap_out, core_ranked_popular_out,
                                        core_ranked_distinct_out):
    '''
    The purpose of this file is so that the user can see why each entity scored
    highly; the user can see the ranked entities.
    Get the ranked entities.
    Get the ranked synonyms.
    '''
    # Initialize the class
    IES = InspectEntityScores(all_caseolap_scores_in, id2syns_in, remove_syns_in, \
                              pmid_syn_count_in, all_cat2pmids_in)

    '''ranked_synonyms'''
    # Export the ranked synonyms
    IES.get_ranked_synonyms_found(all_cat2pmids_in, pmid_syn_count_in, all_ranked_syns_out)

    '''Ranked Entities (CaseOLAP Score: Popular + Distinct)'''
    # Sort entities by CaseOLAP scores
    IES.sort_all_scores(all_caseolap_scores_in)

    # Display the proportion entities found / entities searched for
    IES.prop_entities_found()

    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()

    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(all_ranked_caseolap_out, score_component='CaseOLAP Score')

    '''Ranked Entities (Popular)'''
    # Initialize the class
    IES = InspectEntityScores(all_popular_scores_in, id2syns_in, remove_syns_in, \
                              pmid_syn_count_in, all_cat2pmids_in)

    IES.get_ranked_synonyms_found(all_cat2pmids_in, pmid_syn_count_in, all_ranked_syns_out,
                                  export=False)

    # Sort entities by CaseOLAP scores
    IES.sort_all_scores(all_popular_scores_in)

    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()

    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(all_ranked_popular_out, score_component='Popularity Score')

    '''Ranked Entities (Distinct)'''
    # Initialize the class
    IES = InspectEntityScores(all_distinct_scores_in, id2syns_in, remove_syns_in, \
                              pmid_syn_count_in, all_cat2pmids_in)

    # Get the ranked synonyms
    IES.get_ranked_synonyms_found(all_cat2pmids_in, pmid_syn_count_in, all_ranked_syns_out, \
                                  export=False)

    # Sort entities by CaseOLAP scores
    IES.sort_all_scores(all_distinct_scores_in)

    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()

    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(all_ranked_distinct_out, score_component='Distinctiveness Score')

    # Initialize the class
    IES = InspectEntityScores(core_caseolap_scores_in, id2syns_in, remove_syns_in, \
                              pmid_syn_count_in, core_cat2pmids_in)

    '''ranked_synonyms'''
    # Export the ranked synonyms
    IES.get_ranked_synonyms_found(core_cat2pmids_in, pmid_syn_count_in, core_ranked_syns_out)

    '''Ranked Entities (CaseOLAP Score: Popular + Distinct)'''
    # Sort entities by CaseOLAP scores
    IES.sort_core_scores(core_caseolap_scores_in)

    # Display the proportion entities found / entities searched for
    IES.prop_entities_found()

    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()

    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(core_ranked_caseolap_out, score_component='CaseOLAP Score')

    '''Ranked Entities (Popular)'''
    # Initialize the class
    IES = InspectEntityScores(core_popular_scores_in, id2syns_in, remove_syns_in, \
                              pmid_syn_count_in, core_cat2pmids_in)

    IES.get_ranked_synonyms_found(core_cat2pmids_in, pmid_syn_count_in, core_ranked_syns_out,
                                  export=False)

    # Sort entities by CaseOLAP scores
    IES.sort_core_scores(core_popular_scores_in)

    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()

    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(core_ranked_popular_out, score_component='Popularity Score')

    '''Ranked Entities (Distinct)'''
    # Initialize the class
    IES = InspectEntityScores(core_distinct_scores_in, id2syns_in, remove_syns_in, \
                              pmid_syn_count_in, core_cat2pmids_in)

    # Get the ranked synonyms
    IES.get_ranked_synonyms_found(core_cat2pmids_in, pmid_syn_count_in, core_ranked_syns_out, \
                                  export=False)

    # Sort entities by CaseOLAP scores
    IES.sort_core_scores(core_distinct_scores_in)

    # Save ranked entites, their synonyms, and their synonym counts
    IES.get_entity_syn_counts()

    # Export & display the ranked entities, their synonyms, and synonym counts
    IES.rank_each_category_score(core_ranked_distinct_out, score_component='Distinctiveness Score')


