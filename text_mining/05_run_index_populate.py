'''
The purpose of this file is to populate the ElasticSearch index.
Make sure ElasticSearch is running.
'''
import time, re, sys, os, json
from collections import defaultdict
from elasticsearch import Elasticsearch
from caseolap._05_index_populate import *



'''
Parameters
'''
# Input
parsed_text = '../../../caseolap/data/pubmed.json'#'data/pubmed.json' # Parsed publications
index_populate_config = json.load(open('./config/index_populate_config.json')) 

# Output
logfile_path = './log/indexing_log.txt'            # Reports progress on indexing

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