'''
The purpose of this file is to populate the ElasticSearch index.
Make sure ElasticSearch is running.
'''
import time, re, sys, os, json
from collections import defaultdict
from elasticsearch import Elasticsearch
# setting path
sys.path.append('..')

from text_mining.caseolap._05_index_populate import *



'''
Parameters
'''
# Input
root_dir = '/caseolap_lift_shared_folder/'
data_dir = os.path.join(root_dir,'data')
parsed_text = os.path.join(data_dir,'pubmed.json') #'data/pubmed.json' # Parsed publications
index_populate_config_file = os.path.join(root_dir,'config/index_populate_config.json')
index_populate_config = json.load(open(index_populate_config_file))

# Output
logfile_path = os.path.join(root_dir,'log/indexing_log.txt')            # Reports progress on indexing

# Names of the index you want to create
index_name = 'pubmed_lift'
type_name = 'pubmed_meta_lift'



'''
Main Code
'''
if __name__ == '__main__':
    # Open the log file
    logfile = open(logfile_path, 'w')

    # Populate the index
    populate_index(parsed_text, logfile, index_name, type_name, index_populate_config) 

    # Close the log file
    logfile.close()
