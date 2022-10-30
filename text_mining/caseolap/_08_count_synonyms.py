import sys, json, time, os
from elasticsearch import Elasticsearch
from elasticsearch_dsl import Search, Q
from multiprocessing import cpu_count, Process

class CountSynonyms(object):
    
    def __init__(self, entity_dict_path, textcube_pmid2category):
        '''
        FUNCTION:
        Initialize class attributes.
        - entity_dict: Dictionary. ID to synonyms/names.
        - synonyms: List. Entity synonyms/names. 
        - pmid_and_cat: List of lists [[PMID1:Cat_0],...,[PMID_N:Cat_n]]
        - relevant_pmids: PMIDs for publications studying the categories
          of interest.
        - pmid_syn_cnt: Dictionary of dictionaries. Synonyms to PMIDs to 
          Synonym Counts. 
        
        PARAMS:
        - entity_dict_path: Input file path. Where the entity_dict is stored.
        - textcube_pmid2category: Input file path. Maps PMIDs to category.
        '''
        # Make entity dict
        self.id2syns = json.load(open(entity_dict_path))
        # Get all synonyms
        synonyms = []
        for syns in list(self.id2syns.values()):
            synonyms += syns
        synonyms = list(set(synonyms))
        self.synonyms = [syn.replace('_',' ').strip('\n') for syn in synonyms]
            
        # PMID-Category ([[PMID1:Cat_0],...,[PMID_N:Cat_n]])
        self.pmid_and_cat = json.load(open(textcube_pmid2category))
        
        # PMIDs of interest
        self.relevant_pmids = set(map(lambda x: x[0], self.pmid_and_cat))
        
        # Synonym Count per PMID
        self.pmid_syn_cnt = dict()
    
    

    '''Find and count synonyms (query)'''
    
    def count_synonym(self, b_id, batch, key, es, logfile,\
                      pmid_syn_cnt_path, index_name, start_year, end_year): 
        '''
        FUNCTION:
        Count batch of synonyms in all indexed publications. (Called in syn_search)
        
        PARAMS:
        - syn: the entity synonym searched for
        - key: 'abstract' or 'full_text' indicating the sections to be searched.
          'abstract' searches the title and abstract. 'full_'text' also searches 
          the full text
        - es: ElasticSearch()
        - temp_logfile: Where the threads' log file parts are temporarily stored
          before being merged later
        - pmid_syn_cnt_path: Path. Nested dictionary mapping pmids to synonyms 
          to counts
        - index_name: Name of the index being searched
        - start_year (int): The earliest year of publications you want to look at
        - end_year (int): The latest year of publications you want to look at
        '''
        
        # Open output path
        temp_logfile = logfile[:-4]+'_'+str(b_id)+'.txt'
        temp_pmid_syn_cnt_path = pmid_syn_cnt_path[:-4]+'_'+str(b_id)+'.txt'
        
        with open(temp_pmid_syn_cnt_path,'w') as fout1,\
               open(temp_logfile,'w') as fout2:
        
            # For each synonym from a list of all synonyms
            total_size = len(batch)
            for syn_index, syn in enumerate(batch):
                                
                # Print progress
                if b_id == 1:
                    print('Batch 1 Progress:', syn_index,'/',total_size, end='\r')
                
                # Define the query
                if key == "abstract":
                    query = Q("match_phrase", title=syn) |\
                            Q("match_phrase", abstract=syn)
                elif key == "full_text":
                    query = Q("match_phrase", abstract=syn) | \
                            Q("match_phrase", title=syn) | \
                            Q("match_phrase", full_text=syn)
                elif key == 'full_text_no_methods':
                     query = Q("match_phrase", abstract=syn) | \
                             Q("match_phrase", title=syn) | \
                             Q("match_phrase", introduction=syn) | \
                             Q("match_phrase", results=syn) | \
                             Q("match_phrase", discussion=syn)
                else:
                    print("ERROR: Key is neither 'abstract' nor 'full_text'")
                    sys.exit


                # Perform the query
                s = Search(using=es, index=index_name).\
                params(request_timeout=300).query(query) 


                # Analyze the query results (synonym counts)
                syn_cnts_all_pubs = 0 # Total hits for a synonym

                # Synonym: mid-phrase hyphens with spaces, keep hyphens at end of word
                orig_syn = syn
                syn = (syn+' ').replace('- ','@$#!').replace('-',' ').replace('@$#!','- ').strip()
                
                
                # For each synonym-containing publication
                for hit in s.scan():
                    
                    # Proceed if the publication topic is of interest
                    pmid = str(hit.pmid)
                    if pmid in self.relevant_pmids:
                                                
                        # Check if publication is in a year of interest
                        if int(hit.year) not in range(start_year, end_year+1):
                            continue      
                    
                        # Document sections (e.g., title, abstract, full_text)
                        sections = []
    
                        # Replace hyphens not at end of word (w/spaces)
                        if type(hit.abstract) == str: 
                            hit.abstract = hit.abstract.replace('- ','@$#!').replace('-',' ').replace('@$#!','- ').strip()
                            sections.append(hit.abstract)
                        if type(hit.title) == str:    
                            hit.title = hit.title.replace('- ','@$#!').replace('-',' ').replace('@$#!','- ').strip()
                            sections.append(hit.title)
                        
                        if key == 'full_text':
                            if type(hit.full_text) == str: 
                                hit.full_text = hit.full_text.replace('- ','@$#!').replace('-',' ').replace('@$#!','- ').strip()
                                sections.append(hit.full_text)
                        
                        if key == 'full_text_no_methods':
                            if type(hit.introduction) == str: 
                                hit.introduction = hit.introduction.replace('- ','@$#!').replace('-',' ').replace('@$#!','- ').strip()
                                sections.append(hit.introduction)  
                            if type(hit.results) == str: 
                                    hit.results = hit.results.replace('- ','@$#!').replace('-',' ').replace('@$#!','- ').strip()
                                    sections.append(hit.results) 
                            if type(hit.discussion) == str: 
                                    hit.discussion = hit.discussion.replace('- ','@$#!').replace('-',' ').replace('@$#!','- ').strip()
                                    sections.append(hit.discussion)      
                        

                        # Count the synonym in each section of the publication
                        syn_cnt_this_pub = 0 
                        for section in sections:
                            syn_cnt_this_pub += section.count(syn)
                        syn_cnts_all_pubs += syn_cnt_this_pub

                        # Save the pmid->syn->counts (thread-safe manner)
                        fout1.write(pmid+'|'+orig_syn+'|'+str(syn_cnt_this_pub)+'\n')

                # Save counts->synonyms (thread-safe manner)
                if syn_cnts_all_pubs > 0:
                    fout2.write('#counts:'+str(syn_cnts_all_pubs)+\
                                  '|'+orig_syn+'\n') 
       
    
    
    
    
    ### Search and count synonyms: to optimize and find count from indexer ###
    def synonym_search(self, key, logfile, syn_pmid_count_in, index_name, 
                       start_year, end_year):          
        '''
        FUNCTION:
        Uses threading to call count_synonym, searching for synonyms in the indexed
        publications.
        
        PARAMS:
        - key: 'abstract' or 'full_text' indicating the sections to be searched.
         'abstract' searches the title and abstract. 'full_'text' also 
          searches the full text
        - logfile: Where the log will be stored. Stores entity count.
        - syn_pmid_count_in: Path. Nested dictionary mapping pmids to synonyms 
          to counts. Aka pmid_syn_cnt.
        - index_name: Name of the index being searched. 
        - start_year (int): The earliest year of publications you want to look at
        - end_year (int): The latest year of publications you want to look at
        '''
        
        # Initialize variables
        print("Synonym count is running.....")
        start_time = time.time()
        es = Elasticsearch(timeout=300)
        procs = cpu_count()
        
        # Clear old output
        for b_id in range(procs):
            temp_path = syn_pmid_count_in[:-4]+'_'+str(b_id)+'.txt'
            open(temp_path,'w')
        
        # List of batches for multiprocessing
        batches = [[] for i in range(procs)] 
        
        # Create batches and send to multiprocessing
        for i, syn in enumerate(self.synonyms):
            
            # Add synonym to a batch
            b_id = i%procs
            batches[b_id].append(syn)
            
        # Create a list of jobs
        print("Running jobs. "+\
              str(round(time.time()-start_time, 1))+' seconds')
        jobs = []
        for b_id, batch in enumerate(batches):
            jobs.append(Process(target = self.count_synonym, 
                                args = [b_id, batch, key, es, logfile, 
                                        syn_pmid_count_in, index_name,
                                       start_year, end_year]))
        
        # Run the jobs
        for j in jobs: j.start()
        for j in jobs: j.join()
        print(len(self.synonyms),' synonyms successfully counted! '+\
                  str(round(time.time()-start_time, 1))+' seconds')
        
    
    '''Finalize'''
    
    def merge_logfiles(self, logfile):
        '''
        FUNCTION:
        - There are temporary files storing counts and synonyms. These
          are from different processes in the multiprocessing step. 
          This function merges those files into one. 
          
        PARAMS:
        - logfile: The final merged logfile
        '''
        # Import temp log files (#counts:_|synonymn)
        procs = cpu_count()
        lines = []
        for b_id in range(procs):
            temp_logfile = logfile.split('.txt')[0]+'_'+str(b_id)+'.txt'
            try: 
                lines += open(temp_logfile,'r').readlines()
                os.remove(logfile.split('.txt')[0]+'_'+str(b_id)+'.txt')
            except: 
                continue
        
        # Export to one log file
        open(logfile, 'w').write('')
        for line in lines:
            open(logfile, 'a').write(line)
    
    
    
    
    def get_synonyms_pmid_counts(self, syn_pmid_count, pmid_syn_count_out):
        '''
        FUNCTION:
        Merges temp files with synonym counts (files produced previously).
        Makes {PMID: {synonym:count,...}
        
        PARAMS:
        - syn_pmid_count: Input file. PMIDs per synonym.
          i.e., {synonym : {pmid:count,...}, ...
        - pmid_syn_count_out:  Output file. Synonyms per PMID.
          i.e., {PMID: {synonym:count,...}, ...
        '''
        psc = dict() 
        procs = cpu_count()
        
        ### Merge synonyms-per-pmid files
        for b_id in range(procs):
            syn_pmid_temp = syn_pmid_count.split('.txt')[0]+'_'+str(b_id)+'.txt'
            
            with open(syn_pmid_temp,'r') as fin:
                for line in fin:
                    line = line.split('|')
                    pmid, syn, cnt = line[0], line[1], int(line[2])
                    if cnt == 0:
                        continue
                    psc[pmid][syn] = psc.setdefault(pmid,{}).setdefault(syn,0)+cnt
            os.remove(syn_pmid_temp)

        
        ### Export synonym-per-pmid files
        self.pmid_syn_cnt = psc
        json.dump(self.pmid_syn_cnt, open(pmid_syn_count_out,'w'))
            
    
    
    
    ### Gets "PMID category category" ###
    def save_synfound_pmid2cat(self, synfound_pmid2cat):
        '''
        FUNCTION:
        Saves mapping from PMIDs with synonyms found --> Category number
        
        PARAMS:
        - synfound_pmid2cat: Output file. PMID->CategoryNum
        '''
        ### Write header
        with open(synfound_pmid2cat, "w") as fout:
            fout.write("doc_id\tlabel_id\n")
            
            ### Write PMID--->CategoryNumber
            for pmid, cat in self.pmid_and_cat:
                if pmid in self.pmid_syn_cnt:
                    fout.write(str(pmid) + "\t" + str(cat) + "\n") 
                    
                    
                    
    
    def finish_synonym_search(self, logfile, syn_pmid_count, \
                              pmid_syn_count_out, synfound_pmid2cat):
        '''
        FUNCTION:
        Finishes the synonym count process: 
        - Merges temp log files into one.
        - Saves synonym-pmid-count mappings
        - Deletes old temp files
        
        PARAMS:
        - syn_pmid_count: Path. Nested dictionary mapping pmids to 
          synonyms to counts.
          Aka pmid_syn_cnt.
        - pmid_syn_count_out: Output path. Nested dict of pmid to synonym to count.
        '''
        
        ### Merge logfiles
        self.merge_logfiles(logfile)
        
        
        ### Merge and export synonyms per PMIDs
        self.get_synonyms_pmid_counts(syn_pmid_count, pmid_syn_count_out)
        
                
        ### Export pmid--->category number
        self.save_synfound_pmid2cat(synfound_pmid2cat)
