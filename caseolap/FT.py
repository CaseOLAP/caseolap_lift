import json
import argparse
import requests as req
import os
from multiprocessing import cpu_count, Process
from elasticsearch import Elasticsearch


def get_relevant_categorized_pmids():
    '''
    FUNCTION:
    - Get the PMIDs from the textcube, i.e., the PubMed articles
      that are considered in this study and are categorized to
      be in your studied categories based on their MeSH labels
    '''
    # PMIDs to label
    list_of_list_of_pmids = json.load(open('data/textcube_category2pmid.json'))
    relevant_pmids = set()
    for pmid_list in list_of_list_of_pmids:
        relevant_pmids = relevant_pmids.union(set(pmid_list))
    relevant_pmids = list(relevant_pmids)
    
    return relevant_pmids


def multiprocess_a_list(thelist, the_function):
    '''
    FUNCTION:
    - This takes a list of strings and splits it into input
      for separate processes. The processes then output 
      their results to temp files which are then merged. 
    
    PARARMS:
    - thelist: the list to be split into input for 
      a multiprocessing function
    - the_function: the function that will use the list
      as input for multiprocessing
    '''
    
    # How many processors can be used
    procs = cpu_count()

    # List of batches for multiprocessing
    batches = [[] for i in range(procs)]   

    # Length of input dictionary 
    tot = len(thelist)

    # Create batches and send to multiprocessing
    for i, item in enumerate(thelist):

        # Add synonym to a batch
        b_id = i%procs
        batches[b_id].append(item)

    # Create a list of jobs 
    print("Running jobs...")
    jobs = []
    for b_id, batch in enumerate(batches):
        jobs.append(Process(target = the_function, \
                            args = [b_id, batch]))

    # Run the jobs
    for j in jobs: j.start()
    for j in jobs: j.join()
    print('Done!')
    
    
def switch_dictset_to_dictlist(the_dict):
    '''
    FUNCTION:
    - Make a new dictionary with values as lists 
      instead of values as sets
      
    PARAMS:
    - the_dict: The initial dict with values of sets
    '''
    
    dictlist = dict()
    
    for k in the_dict.copy():
        dictlist[k] = list(the_dict[k])
        
    return dictlist



def map_pmc_pmids(batch_id, pmids):
    '''
    FUNCTION:
    - Map a batch of PMIDs to PMCs
    - API Source: https://www.ncbi.nlm.nih.gov/pmc/tools/id-converter-api/
    
    PARAMS:
    - batch_id: The batch ID
    - pmids (list): the batch of PMIDs
    '''
    url = 'https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids='
    tot_reqs = int(len(pmids)/200)+1
    temp_outfile = 'data/pmc_pmid_temp_'+str(batch_id)+'.txt'
    
    with open(temp_outfile,'w') as fout:
        for i in range(0, int(len(pmids)/200)+1):

            # Print progress
            if batch_id == 1:
                print(str('Progress of batch 1: '+str(i)+'/'+str(tot_reqs)), end='\r')  

            # Submit batch of up to 200 PMIDs
            if i < int(len(pmids)/200)+1:
                r = req.get(url+','.join(pmids[i*200:(i+1)*200])+'&format=json')
            else:
                r = req.get(url+','.join(pmids[i*200:])+'&format=json')

            # Map all PMID-PMCs
            records = r.json()['records']
            for record in records:
                try:
                    pmc = record['pmcid']
                    pmid = record['pmid']
                    fout.write(pmid+'|'+pmc+'\n')
                except:
                    break
                    
                    
def merge_pmid_pmc_mapping(pmid2pmc, pmc2pmid):
    '''
    FUNCTION:
    Merge the files output by the parallel processes
    
    PARAMS:
    - pmid2pmc: Dictionary to map the PMIDs to the PMCs
    - pmc2pmid: Dictionary to map the PMCs to the PMIDs
    '''
    procs = cpu_count()

    for batch_id in range(procs):
        temp_outfile = 'data/pmc_pmid_temp_'+str(batch_id)+'.txt'
        with open(temp_outfile) as fin:
            for line in fin:
                line = line.split('|')
                pmid = line[0]
                pmc = line[1].replace('PMC','').strip()
                pmid2pmc.setdefault(pmid, set()).add(pmc)
                pmc2pmid.setdefault(pmc, set()).add(pmid)
            
    return pmid2pmc, pmc2pmid


def finish_pmc2pmid_mapping():
    '''
    FUNCTION:
    - Merge the separate files to create a PMCID to PMID mapping
    '''
    pmid2pmc, pmc2pmid = dict(), dict()
    pmid2pmc, pmc2pmid = merge_pmid_pmc_mapping(pmid2pmc, pmc2pmid)

    # Checking that the PMID-PMC is 1-to-1
    for pmid,pmc in pmid2pmc.items():
        assert len(pmc) == 1
    for pmc,pmid in pmc2pmid.items():
        assert len(pmid) == 1


    # Making the PMID-PMC not have list values
    pmc2pmid = switch_dictset_to_dictlist(pmc2pmid)
    pmid2pmc = switch_dictset_to_dictlist(pmid2pmc)
    for pmid,pmc in pmid2pmc.copy().items():
        pmid2pmc[pmid] = pmc[0]
    for pmc,pmid in pmc2pmid.copy().items():
        pmc2pmid[pmc] = pmid[0]  


    '''Export the PMC-PMID mappings'''
    json.dump(pmc2pmid, open('data/pmc2pmid.json','w'))
    json.dump(pmid2pmc, open('data/pmid2pmc.json','w'))

    assert len(pmid2pmc) == len(pmc2pmid)
    print('Mapped', len(pmid2pmc), 'IDs')
    
    
    
def get_pmids_pmcs(save_pmcs_to):
    ''' 
    FUNCTION:
    - Save the PMIDs-PMCIDs mappings, PMCIDs, and PMIDs
    
    PARAMS:
    - save_pmcs_to (str): path where the PMCIDs are saved
    '''
    
    '''Get PMIDs used in this study'''
    relevant_pmids = get_relevant_categorized_pmids()
    json.dump(relevant_pmids, open('data/relevant_pmids.json','w'))

    '''Map PMIDs->PMCIDs'''
    pmids = json.load(open('data/relevant_pmids.json'))
    multiprocess_a_list(pmids, map_pmc_pmids)
    finish_pmc2pmid_mapping()

    '''Save Relevant PMCs'''
    pmc2pmid = json.load(open('data/pmc2pmid.json'))
    pmcs = list(pmc2pmid.keys())
    json.dump(pmcs, open(save_pmcs_to,'w'))
    
    
def get_header_inds(low_text):
    '''
    FUNCTION:
    - Get the indices of the section headers in a PubMed Central full text article.
    
    PARAMS:
    - low_text: The lowercased PubMed Central article separated into lines, 
      i.e., entries in a list
    '''
    intro_header_idx,methods_header_idx,res_header_idx,dis_header_idx = -1,-1,-1,-1

    intro_headers = ['introduction\n', 'intro\n']
    for intro_header in intro_headers:
        if intro_header in low_text:
            intro_header_idx = low_text.index(intro_header)

    methods_headers = ['methods\n', 'materials and methods\n', 'materials & methods']
    for methods_header in methods_headers:
        if methods_header in low_text:
            methods_header_idx = low_text.index(methods_header)

    res_headers = ['results\n','results and discussion\n',
                       'results & discussion\n']
    for res_header in res_headers:
        if res_header in low_text:
            res_header_idx = low_text.index(res_header)

    dis_headers = ['discussion\n', 'discussion and conclusions\n']
    for dis_header in dis_headers:
        if dis_header in low_text:
            dis_header_idx = low_text.index(dis_header)

    header_inds = [intro_header_idx, methods_header_idx, 
                   res_header_idx, dis_header_idx]
    
    return header_inds, intro_header_idx, methods_header_idx, res_header_idx, dis_header_idx




def find_section(section_idx, ordered_inds, fulltext):
    '''
    FUNCTION:
    - Using the indices of the section headers, this returns the section text
      for your section of interest (e.g., Intro)
    
    PARAMS:
    - section_idx: The current section's index
    - ordered_inds: The section indices ordered from first to last
    - fulltext: The full text of the PubMed document
    '''
    section_idx_order = ordered_inds.index(section_idx)

    if section_idx_order < len(ordered_inds)-1:
        next_idx = ordered_inds[section_idx_order+1]
        section = fulltext[section_idx+1:next_idx]
    else:
        section = fulltext[section_idx+1:]
        
    section = ' '.join(section)
    section = section.replace('  ','').replace('\n',' ')  
    
    return section


def remove_references_authors_abbreviation_sections(text):
    '''
    FUNCTION:
    - Removes the References section and the Abbreviation sections
      at the end of the documents/articles
    PARAMS:
    - text (list of str): The document/article full text
    '''
    
    if 'Abbreviations\n' in text:
        abbr_idx = text.index('Abbreviations\n')
        text = text[:abbr_idx]
    
    if 'References\n' in text:
        ref_idx = text.index('References\n')
        text = text[:ref_idx]
        
    if '==== Refs' in text:
        ref_idx = text.index('==== Refs')
        text = text[:ref_idx]
        
    if '==== Body\n' in text:
        body_idx = text.index('==== Body\n')
        text = text[body_idx:]

    return text


def get_pmid2sections_dict():
    '''
    FUNCTION:
    - Access the saved full text and split it into sections.
      Return a dictionary of {PMID : {Intro: ..., Methods: ..., 
      Results: ..., Discussion: ...}}
    '''
    #path = 'pubmed_central_text_download/data'
    path = '../pmc_full_texts'
    pmid2sections, pmc2fulltext = dict(), dict()
    pmc2pmid = json.load(open('data/pmc2pmid.json'))

    # Iterate through PMC full text file directory
    for num_files, filename in enumerate(os.listdir(path)):
        print(str(num_files)+' PMIDs', end='\r')

        # Read in the PMC full text
        with open(os.path.join(path, filename)) as fin:
            try: unprocessed_text = fin.readlines()
            except: continue

            # Get text, PMC, PMID
            text = remove_references_authors_abbreviation_sections(unprocessed_text)
            low_text = [line.lower() for line in text]
            try:
                pmc = filename.split('.txt')[0].split('PMC')[1]
                pmid = pmc2pmid[pmc]
            except:
                print('Couldn\'t map PMCID', pmc, 'to PMID')
                
            # Get indices
            head_inds,intro_idx,method_idx,res_idx,dis_idx = get_header_inds(low_text)
            sort_head_inds = sorted(head_inds)

            # Get sections: Introduction, Methods, Results, Discussion, Full Text
            pmid2sections[pmid] = {'Introduction':'', 'Methods':'', 
                                   'Results':'', 'Discussion':''}    
            if intro_idx != -1:
                pmid2sections[pmid]['Introduction'] = find_section(intro_idx, 
                                                                   sort_head_inds, 
                                                                   low_text)
            if method_idx != -1:
                pmid2sections[pmid]['Methods'] = find_section(method_idx, 
                                                              sort_head_inds, 
                                                              low_text)
            if res_idx != -1:
                pmid2sections[pmid]['Results'] = find_section(res_idx, 
                                                              sort_head_inds, 
                                                              low_text)
            if dis_idx != -1:
                 pmid2sections[pmid]['Discussion'] = find_section(dis_idx, 
                                                                  sort_head_inds,
                                                                  low_text)
            pmid2sections[pmid]['full_text'] = ' '.join(unprocessed_text)
    
    json.dump(pmid2sections, open('data/pmid2fulltext_sections.json','w'))
    print('Saved', len(pmid2sections), 'full text sections out of', num_files, 
          'attempted')
    
    return pmid2sections


def get_number_of_sections(pmid2sections):
    '''
    FUNCTION:
    - Gets counts on how many intros, methods, results, discussions,
      and full text were found.
      
    PARAMS:
    - pmid2sections: PubMed mappings to sections
    '''
    intros, methods, results, discussion, full_text = 0,0,0,0,0

    for sections in pmid2sections.values():
        if sections['Introduction'] != '':
            intros += 1
        if sections['Methods'] != '':
            methods += 1
        if sections['Results'] != '':
            results += 1    
        if sections['Discussion'] != '':
            discussion += 1
        if sections['full_text'] != '':
            full_text += 1

    print(str(len(pmid2sections))+' PMIDs searched\n'+\
          str(full_text)+' Full text found (i.e., more than title and abstract)\n'+\
          str(intros)+' Introductions\n'+\
          str(methods)+' Methods\n'+\
          str(results)+' Results\n'+\
          str(discussion)+' Discussion\n')
    
    
def update_indexing_with_article_sections(index_name, index_type):
    '''
    FUNCTION:
    - Iterate through each publication with full text
    
    PARAMS:
    - index_name: Name of Elasticsearch index
    - index_type: Name of Elasticsearch index type
    '''
    es = Elasticsearch()
    pmid2sections = json.load(open('data/pmid2fulltext_sections.json'))
    total_pmids = len(pmid2sections)
    for num_pmids, (pmid, sections) in enumerate(pmid2sections.items()):
        print(str(num_pmids)+'/'+str(total_pmids)+\
              ' publications\' sections indexed', end='\r')
        
        sections = {'introduction': sections['Introduction'].replace('\n', '')\
                    .replace('\t','').replace('  ', ' '),
                    'methods': sections['Methods'].replace('\n','')\
                    .replace('\t','').replace('  ', ' '),
                    'results': sections['Results'].replace('\n', '')\
                    .replace('\t','').replace('  ', ' '),
                    'discussion': sections['Discussion'].replace('\n', '')\
                    .replace('\t','').replace('  ', ' '),
                    'full_text': sections['full_text'].replace('\n', '')\
                    .replace('\t','').replace('  ', ' ')}

        # Update each publication's index
        es.update(index = index_name, 
                  id = pmid,
                  doc_type = index_type,
                  doc = sections)

        


        
''' MAIN '''
# Parser
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index_name', default='pubmed_lift', type=str)
parser.add_argument('-it', '--index_type', default='pubmed_meta_lift', type=str)
parser.add_argument('-d', '--download_ft', default=False, type=bool)
parser.add_argument('-p', '--print_section_counts', default=True, type=bool)
args = parser.parse_args()

index_name = args.index_name  # ElasticSearch index name
index_type = args.index_type  # ElasticSearch index type name
download_full_text = args.download_ft # Download full text from PMC
print_section_counts = args.print_section_counts
pmc_path = 'pubmed_central_text_download/pmcs.json' # Save full text here



# Get PMIDS, PMCIDs, PMID-PMCID mappings
print('Getting PMCIDs, PMIDs, and PMID-PMCID mappings')
get_pmids_pmcs(save_pmcs_to = pmc_path)



''' Download full texts '''
# This downloads the full texts from PubMed
download_full_text = False
if download_full_text:
    print('Downloading full text articles from PubMed Central')
    os.system('python pubmed_central_text_download/download_pmc_fulltext.py')

    
    
''' Index full text split into sections '''
# Get article sections (e.g., introduction, methods)
print('Sectioning articles into introduction, methods, results, and discussion')
pmid2sections = get_pmid2sections_dict()

# Get number of article sections found
if print_section_counts:
    get_number_of_sections(pmid2sections)

# Update the ElasticSearch index
print('Updating the indexed articles with their sections')
update_indexing_with_article_sections(index_name = index_name, 
                                      index_type = index_type)