import time, re, sys, os, json
from collections import defaultdict
from elasticsearch import Elasticsearch


def populate_index(parsed_text_infile, logfile, index_name, type_name, 
                   index_populate_config):
    '''
    FUNCTION:
    - This populates the ElasticSearch index with information from PubMed
      publications: PMID, title, abstract, MeSH terms, date, author, 
      location, journal, and full text if available.. 
    
    PARAMS:
    - parsed_text_infile (JSON): The input file of parsed publications
      (keys and values for
      title:'what the title actually is', etc.)
    - logfile: Output file, where the progress will be printed.
    - index_name (str): The name you chose for the index
    - type_name (str): The name you chose for the index type
    - index_populate_config (JSON): Indicates which fields (e.g., title) to 
      index.
    '''
    
    es = Elasticsearch()
    batch_num = 0
    
    with open(parsed_text_infile, 'r') as fin: 
        start = time.time()
        BATCH_SIZE = 500    # number of document processed in each batch index
        batch_data = []     # data in batch index

        
        '''Collect a batch of documents'''
        # Iterate through all documents
        for doc_count, document in enumerate(fin): 

            # Get information on each document (e.g., title, abstract, date)
            paper_info = json.loads(document.strip())
            data_dict = {}

            # PMID, title, abstract, full text
            data_dict['pmid'] = paper_info.get('PMID', '-1')                     
            data_dict['title'] = paper_info.get('ArticleTitle')                   
            data_dict['abstract'] = paper_info.get('Abstract', '')               
            data_dict['full_text'] = paper_info['full_text']    
            data_dict['introduction'] = ''
            data_dict['methods'] = ''
            data_dict['results'] = ''
            data_dict['discussion'] = ''
            
            # Date
            if index_populate_config['date']:                                   
                data_dict['date'] = str(paper_info['PubDate'])
                
                # Year
                year = paper_info['PubDate']['Year']
                try: 
                    year = int(year)
                except:
                    try:
                        year = int(hit.date['MedlineDate'].split(' ')[0])
                    except:
                        year = -1
                data_dict['year'] = str(year)

            # MeSH
            if index_populate_config['MeSH']:                                   
                data_dict['MeSH'] = paper_info['MeshHeadingList']

            # Location (where paper was published)
            if index_populate_config['location']:                               
                data_dict['location'] = paper_info['Country']

            # Authors
            if index_populate_config['author']:                                 
                data_dict['author'] = paper_info['AuthorList']

            # Journal
            if index_populate_config['journal']:                                 
                data_dict['journal'] = paper_info['Journal']

            op_dict = {'index': {
                        '_index': index_name,
                        '_type': type_name,
                        '_id': data_dict['pmid']}}
            
            # Put current data into the batch 
            batch_data.append(op_dict)
            batch_data.append(data_dict) 

            
            '''Index a batch of documents'''
            # Bulk indexing
            if doc_count % BATCH_SIZE == 0 and doc_count != 0:
                
                # Index
                es.bulk(index = index_name, body = batch_data, request_timeout = 500)
                batch_data = []
                
                # Print progress
                elapsed_time = round(time.time() - start, 3)
                msg = str(doc_count)+' documents indexed in '\
                     +str(round(elapsed_time/60,3))+' minutes'
                logfile.write(msg+'\n')
                batch_num += 1
                if batch_num % 100 == 0:
                    print(msg)
                
        '''Index the remainder in the last batch'''
        # If, at the end of the input file, there are leftover files in the 
        # incomplete batch, index them
        if batch_data != []:

            # Index
            es.bulk(index = index_name, body = batch_data, request_timeout = 500)
            batch_data = []

            # Print progress
            elapsed_time = round(time.time() - start, 3)
            msg = str(doc_count)+' documents indexed in '\
                 +str(round(elapsed_time/60,3))+' minutes'
            logfile.write(msg+'\n')
            print(msg)

            
        '''Print that the process is complete'''
        elapsed_time = round(time.time() - start, 3)
        msg = 'Finished indexing in '+str(round(elapsed_time/60,3))+' minutes'
        logfile.write(msg + '\n')
        print(msg)                    