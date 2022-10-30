'''
The purpose of this file is to produce CaseOLAP scores for the entities
based on their hits in each document (pmid2pcount_path) and the documents'
category (category2pmids_path).
'''
import pandas as pd, numpy as np, matplotlib.pyplot as plt, seaborn as sns, json
from caseolap._12_caseolap_score import *



'''
Parameters
'''
# Input data directories
cat2pmids_path = './data/metadata_category2pmids_2012-2022.json' # {CategoryName:[PMID,...],...}
all_pmid2entity2count_path = 'data/metadata_pmid2entity2count_2012-2022.json' # {PMID:{Entity:Count,...},...}
core_pmid2entity2count_path = 'data/metadata_pmid2entity2count_2012-2022.json' # {PMID:{Entity:Count,...},...}
category_names_path = 'config/textcube_config.json'    # ['CategoryName1',...]

# Output data path
result_dir = 'result/' # Main folder where the results from this section will be stored
all_logFilePath = './log/all_caseolap_score_log.txt' # Logs #PMIDs for each category
core_logFilePath = './log/core_caseolap_score_log.txt' # Logs #PMIDs for each category
all_caseolap_name = 'all_caseolap'  # Name of dataframe/spreadsheet for the caseolap scores
core_caseolap_name = 'core_caseolap'

'''
Main Code
'''
if __name__ == '__main__':

    #### All proteins ####
    logfile = open(all_logFilePath, 'w') 
    category2pmids = json.load(open(cat2pmids_path, 'r'))   
    pmid2entity2count = json.load(open(all_pmid2entity2count_path, 'r'))


    ''' Initial Calculations'''
    # Initialize object with input data
    C = Caseolap(category2pmids, pmid2entity2count, result_dir, category_names_path, logfile)
    
    # Print info on categories and their number of publications
    C.print_categories2pmid(dump = True, verbose =True)
    
    # Map Category to its PMIDs to its Entities to the Entity's Counts
    C.map_category2pmid2entity2count()
    
    # Save all entities 
    C.get_all_entities(dump = True, verbose = True)

    
    ''' Popularity Score (note: relies on some previous sections above)'''
    # Get the entity counts for each category   
    C.get_entity_counts_per_category()
    
    # Maps category to its entities to their counts (includes zero count entities)
    C.category2entity2tf_finder()

    # Calculate the popularity scores for all entities
    C.calculate_all_popularity_scores(dump = True)
    
    
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
    C.calculate_all_distinctiveness_scores(dump = True)
    
    
    '''Final Score'''
    # Calculate CaseOLAP Score (combine popularity & distinctiveness)
    C.calculate_caseolap_score(all_caseolap_name, dump = True)
    
    # Close logfile
    logfile.close()
    
    
    
    
    
    
    #### Core proteins ##### 
    logfile = open(core_logFilePath, 'w') 
    category2pmids = json.load(open(cat2pmids_path, 'r'))   
    pmid2entity2count = json.load(open(core_pmid2entity2count_path, 'r'))


    ''' Initial Calculations'''
    # Initialize object with input data
    C = Caseolap(category2pmids, core_pmid2entity2count_path, result_dir, category_names_path, logfile)
    
    # Print info on categories and their number of publications
    C.print_categories2pmid(dump = True, verbose =True)
    
    # Map Category to its PMIDs to its Entities to the Entity's Counts
    C.map_category2pmid2entity2count()
    
    # Save all entities 
    C.get_all_entities(dump = True, verbose = True)

    
    ''' Popularity Score (note: relies on some previous sections above)'''
    # Get the entity counts for each category   
    C.get_entity_counts_per_category()
    
    # Maps category to its entities to their counts (includes zero count entities)
    C.category2entity2tf_finder()

    # Calculate the popularity scores for all entities
    C.calculate_all_popularity_scores(dump = True)
    
    
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
    C.calculate_all_distinctiveness_scores(dump = True)
    
    
    '''Final Score'''
    # Calculate CaseOLAP Score (combine popularity & distinctiveness)
    C.calculate_caseolap_score(core_caseolap_name, dump = True)
    
    # Close logfile
    logfile.close()
    
    
    
    
    
    
    
    
    
    