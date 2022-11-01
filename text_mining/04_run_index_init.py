'''
The purpose of this file is to initialize the ElasticSearch index
which is where the indexed PubMed text data will be.
This version allows the option to preserve case-sensitivity in the
indexed text.
'''
import json
import os
from elasticsearch import Elasticsearch


'''
Parameters
'''
# Index parameters
root_dir = '/caseolap_lift_shared_folder/'
index_name = 'pubmed_lift'      # Index name (match with 05 index populate file)
type_name = 'pubmed_meta_lift'  # Index type name (match with 05 index populate file)
number_shards = 1          # Set to 1 if no cluster
number_replicas = 0    
case_sensitive = True   # Index the text as case sensitive (True) or lower case (False)

# Input file
index_init_config_file = os.path.join(root_dir,'config/index_init_config.json')


'''
Main Code
'''
if __name__ == '__main__':
    # Load the indexing config file
    index_init_config = json.load(open(index_init_config_file,'r')) 
    
    # Start elasticsearch 
    es = Elasticsearch('https://localhost:9200') 
        
    # Delete the old index if it exists
    if es.indices.exists(index = index_name):
        res = es.indices.delete(index = index_name)
        print('Deleted index:',index_name,'\nResponse:',res,'\n')
    
    
    # Request Body Parameters
    mappings = {type_name: {'properties':index_init_config}}
    settings = {'number_of_shards': number_shards, 'number_of_replicas': number_replicas}

    if case_sensitive == True:
        settings['analysis'] = {'analyzer':{'casesensitive_text':{'type':'custom',
                                                                  'tokenizer':'standard',
                                                                  'filter': ['stop']}}}
    # Create an index 
    res = es.indices.create(index = index_name, settings = settings, mappings = mappings)
    print('Created index:',index_name,'\nResponse:',res)
    
    
