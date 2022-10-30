from caseolap._08_count_synonyms import *
import sys, json, time, os
from elasticsearch import Elasticsearch
from elasticsearch_dsl import Search, Q
from multiprocessing import cpu_count, Process

'''
Parameters
'''
# Input
entity_dict_path = 'input/id2syns.json'  
textcube_pmid2category = 'data/textcube_pmid2category.json'
start_year = 2012    # <-- Change as you see fit for your publications of interest
end_year = 2022      # <-- Same as above comment 

# Intermediary file (produced as output, used as input)
syn_pmid_count = 'data/syn_pmid_count.txt'

# Output 
pmid_syn_count_out = 'data/pmid_synonym_counts_2012-2022.json'   # PMID Syn|Count...Syn|Cnt
synfound_pmid2cat = 'data/synfound_pmid2category_2012-2022.txt'  # PMID--->CategoryNumber
logfile = 'log/synonymcount_log_2012-2022.txt'                   # #hits:Synonym



# Other parameters
index_name = 'pubmed_lift' # Index name
key = 'abstract'        # Choose if searching the abstracts and titles
key = 'full_text'      # Choose if searchines the abstracts, titles, and full text
key = 'full_text_no_methods' # Same as above, excluding abstracts


'''
Main code
'''
if __name__ == '__main__':
    
    # Instantiate the object
    CS = CountSynonyms(entity_dict_path, textcube_pmid2category)
    
    # Search for the synonyms in the indexed text
    CS.synonym_search(key, logfile, syn_pmid_count, index_name, start_year, end_year) 
   
    # Finalize the output files
    CS.finish_synonym_search(logfile, syn_pmid_count,\
                             pmid_syn_count_out, synfound_pmid2cat)
