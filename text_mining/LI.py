import requests as req
import json
import xml.etree.ElementTree as ET
import itertools
import os
import argparse
import sys
sys.path.append('..')
from utils.biomedkg_utils import switch_dictset_to_dictlist
from elasticsearch import Elasticsearch
from multiprocessing import cpu_count, Process


def es_iterate_all_documents(es, index, pagesize=250, scroll_timeout="1m", **kwargs):
    """
    Helper to iterate ALL values from a single index
    Yields all the documents.
    Source: https://techoverflow.net/2019/05/07/elasticsearch-how-to-iterate-scroll-through-all-documents-in-index/
    """
    is_first = True
    while True:
        # Scroll next
        if is_first: # Initialize scroll
            result = es.search(index=index, scroll="1m", **kwargs, 
                               size = pagesize)
            is_first = False
        else:
            result = es.scroll(
                scroll_id = scroll_id,
                scroll = scroll_timeout)
        scroll_id = result["_scroll_id"]
        hits = result["hits"]["hits"]

        # Stop after no more docs
        if not hits:
            break
        
        # Yield each entry
        yield from (hit['_source'] for hit in hits)


def map_disease_mesh_name_to_id(mesh_desc_file = '../data/MeSH/desc2022.xml', 
                                output_folder='../parsed_mappings/MeSH'):
    '''
    FUNCTION:
    - Map the disease MeSH names/terms to the
      MeSH IDs.
    '''
    print('Mapping MeSH disease names to IDs')
    tree = ET.parse(mesh_desc_file)
    root = tree.getroot()   
    
    name2id = dict()

    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')

            # If Tree is a disease
            for tree_number in tree_numbers:
                if tree_number.text.startswith(('C','F03')):
                    tree_number = tree_number.text

                    # ID to Name
                    try:
                        ID = ele.find('DescriptorUI').text
                        name = ele.find('DescriptorName').find('String').text
                        name2id.setdefault(name,set()).add(ID)
                    except:
                        pass
        except:
            continue        

    name2id = switch_dictset_to_dictlist(name2id)
    json.dump(name2id, open(os.path.join(output_folder,'name2id.json'),'w'))


def map_categories2terms(output_folder,
                         meshterms_per_cat_file='../parsed_mappings/MeSH/meshterms_per_cat.json',
                         category_names_file='../config/textcube_config.json'):
    '''
    FUNCTION:
    - Map disease categories -> MeSH terms
    '''
    print('Mapping disease categories to their terms and subcategory terms')
    terms_lol = json.load(open(meshterms_per_cat_file, 'r'))
    category_names = json.load(open(category_names_file, 'r'))
    category2terms = dict()
    for category_name, term_list in zip(category_names, terms_lol):
        category2terms[category_name] = term_list

    json.dump(category2terms,
              open(os.path.join(output_folder, 'category2terms.json'), 'w'))

    return category2terms


def get_mesh_synonyms_api(category2terms, name2id,
                          output_folder='../parsed_mappings/MeSH/'):
    '''
    FUNCTION:
    - Get MeSH synonyms via API
    '''
    category2synonyms = dict()
    total_categories = len(category2terms)
    for category_num, (category, terms) in enumerate(category2terms.items()):
        print(str(category_num) + '/' + str(total_categories), end='\r')

        # Terms in category
        for term in terms:

            mesh_id = name2id[term]
            assert len(mesh_id) == 1
            r_string = 'https://id.nlm.nih.gov/mesh/lookup/details?descriptor=' + mesh_id[0]
            r = req.get(r_string).json()

            # Add term synonyms to the category
            for entry in r['terms']:
                synonym = entry['label']
                category2synonyms.setdefault(category, set()).add(synonym)

    category2synonyms = switch_dictset_to_dictlist(category2synonyms)
    json.dump(category2synonyms,
              open(os.path.join(output_folder, 'category2synonyms.json'), 'w'))

    return category2synonyms


def map_category_to_permuted_mesh_synonyms(category2synonyms,
                                           output_folder='../parsed_mappings/MeSH/'):
    '''
    FUNCTION:
    - permute the MeSH Synonyms because sometimes they're
      written weirdly with commas (e.g., Order, Messed, Up, Is;
      Messed, Order, Up)

    PARAMS:
    - category2synonyms: category -> MeSH synonyms mapping
    '''
    print('Mapping disease categories to permuted MeSH synonyms')

    # Permuted synonym->category dictionary
    category2permuted_synonyms, permuted_synonyms2category = dict(), dict()
    for category, synonyms in category2synonyms.items():
        synonyms = permute_mesh_synonyms(synonyms)
        for synonym in synonyms:
            category2permuted_synonyms.setdefault(category, set()).add(synonym)
            permuted_synonyms2category.setdefault(synonym, set()).add(category)

    category2permuted_synonyms = switch_dictset_to_dictlist(category2permuted_synonyms)
    permuted_synonyms2category = switch_dictset_to_dictlist(permuted_synonyms2category)

    json.dump(category2permuted_synonyms,
              open(os.path.join(output_folder, 'category2permuted_synonyms.json'), 'w'))
    json.dump(permuted_synonyms2category,
              open(os.path.join(output_folder, 'permuted_synonyms2category.json'), 'w'))

    return category2permuted_synonyms, permuted_synonyms2category


def permute_mesh_synonyms(mesh_synonyms):
    '''
    FUNCTION:
    - Permute MeSH synonyms

    PARAMS:
    - mesh_synonyms: list of mesh synonyms
    '''

    temp_mesh_synonyms = list()
    for mesh_synonym in mesh_synonyms:
        if ',' in mesh_synonym:
            mesh_synonym = mesh_synonym.split(', ')
            mesh_synonym_permuted = [' '.join(list(permuted_synonym)) for permuted_synonym in
                                     list(itertools.permutations(mesh_synonym))]
            temp_mesh_synonyms += mesh_synonym_permuted
        else:
            temp_mesh_synonyms.append(mesh_synonym)

    return temp_mesh_synonyms


def get_relevant_categorized_pmids(list_of_list_of_pmids_file='../data/textcube_category2pmid.json'):
    '''
    FUNCTION:
    - Get the PMIDs from the textcube, i.e., the PubMed articles
      that are considered in this study and are categorized to
      be in your studied categories based on their MeSH labels
    - Purpose: For testing against ground truth labels, or for
      if you think you can impute more labels that the NIH missed
      when they were annotating with MeSH terms.
    '''
    # PMIDs to label
    list_of_list_of_pmids = json.load(open(list_of_list_of_pmids_file))
    relevant_pmids = set()
    for pmid_list in list_of_list_of_pmids:
        relevant_pmids = relevant_pmids.union(set(pmid_list))
    relevant_pmids = list(relevant_pmids)

    return relevant_pmids


def get_all_uncategorized_pmids(index_name):
    '''
    FUNCTION:
    - Get the PMIDs of all the PubMed articles that
      don't have MeSH labels.
    - Purpose: For applying to unlabeled PubMed articles

    PARAMS:
    - index_name: The name of the ElasticSearch index
      where all the PubMed articles are indexed/stored
    '''

    es = Elasticsearch()
    relevant_uncategorized_pmids = set()
    for num_pmids, entry in enumerate(es_iterate_all_documents(es, index_name)):

        # Publication's MeSH (if any)
        unlabeled = (len(entry['MeSH']) == 0)

        # Save PMIDs of unlabeled documents
        if unlabeled:
            relevant_uncategorized_pmids.add(entry['pmid'])

            # Print progress
            if num_pmids % 10000 == 0:
                print(str(num_pmids) + ' PMIDs processed', end='\r')

    return list(relevant_uncategorized_pmids)


def get_all_categorized_pmids(index_name, outpath):
    '''
    FUNCTION:
    - Get the PMIDs of all the PubMed articles that
      don't have MeSH labels.
    - Purpose: For applying to unlabeled PubMed articles

    PARAMS:
    - index_name: The name of the ElasticSearch index
      where all the PubMed articles are indexed/stored
    - outpath: Where the categorized PMIDs will go.
    '''

    es = Elasticsearch()
    relevant_categorized_pmids = set()
    for num_pmids, entry in enumerate(es_iterate_all_documents(es, index_name)):

        # Publication's MeSH (if any)
        labeled = (len(entry['MeSH']) > 0)

        # Save PMIDs of unlabeled documents
        if labeled:
            relevant_categorized_pmids.add(entry['pmid'])

            # Print progress
            if num_pmids % 10000 == 0:
                print(str(num_pmids) + ' PMIDs processed', end='\r')

    relevant_categorized_pmids = list(relevant_categorized_pmids)
    json.dump(relevant_categorized_pmids, open(outpath + 'relevant_categorized_pmids.json', 'w'))

    return relevant_categorized_pmids


def get_document_text(entry):
    '''
    FUNCTION: 
    - Get the full text of the PubMed publication

    PARAMS: 
    - entry (dict): parsed indexed publication 
    '''

    # Title
    title = entry['_source']['title']
    if type(title) != str:
        title = ''
    title.replace('\n', ' ').replace('\t', ' ').replace('   ', ' ')

    # Abstract
    abstract = entry['_source']['abstract']
    if type(abstract) != str:
        abstract = ''
    abstract.replace('\n', ' ').replace('\t', ' ').replace('   ', ' ')

    return abstract, title


def ds_label_matching(batch_id, relevant_pmid_batch,
                      index_name, index_type, 
                      label_unlabeled_only,
                      label_labeled_only,
                      label_all,
                      filter_list, 
                      stop_at_this_many_pmids,
                      run, test,
                      output_folder):
    '''
    FUNCTION:
    - Using MeSH term synonyms, label a document with a MeSH term if the 
      lowercased term matches exactly within the text.
    
    PARAMS:
    - relevant_pmid_batch: PMIDs to look for in this batch
    - temp_outfile: temporary output file of PMID | MeSH Term
    '''
    es = Elasticsearch()
    if run:
        rot = 'run'
    elif test:
        rot = 'test'
    else:
        raise Except('Did not specify if running or testing')

    temp_outfile = os.path.join(output_folder, rot+'_temp_labeling'+str(batch_id)+'.txt')
    procs = cpu_count()
    category2permuted_synonyms = json.load(open(os.path.join(output_folder,'category2permuted_synonyms.json'),'r'))

    with open(temp_outfile,'w') as fout, \
    open(temp_outfile[:-4]+'synonym'+'.txt','w') as fout1:

        # Get PubMed article's text
        for num_pmids, pmid in enumerate(relevant_pmid_batch):
            entry = es.get(id = pmid, index = index_name, doc_type = index_type)

            # Print progress and break early
            if batch_id == 1 and num_pmids % 1000 == 0:
                print(str(num_pmids)+' PMIDs processed in one of the batches',\
                      end='\r')
            if num_pmids > stop_at_this_many_pmids/procs:
                break


            # Determine which publications to find labels for
            labeled_meshes = entry['_source']['MeSH']
            labeled = len(labeled_meshes) > 0
            unlabeled = len(labeled_meshes) == 0
            if label_labeled_only and labeled:
                pass
            elif label_unlabeled_only and unlabeled:
                pass
            elif label_all:
                pass
            else:
                continue
            
                        
            # Publication's text (title, abstract, full text if provided)
            abstract, title = get_document_text(entry)
            document_text = title + ' ' + abstract
            title = ' '+title+' '
            title = title.replace(',',' ')
            title = title.replace(':',' ')
            
            # Only consider documents that have certain broader key words
            # E.g., for cardiovascular disease publications, only look 
            # at publications that say "heart" or "cardiac" 
            #dont_label = True
            #for filter_word in filter_words:
            #    if filter_word in document_text.lower():
            #        dont_label = False
            #        break
            #if dont_label:
            #    continue
              
            #print(dont_label, 'dont_label')
            
            #fout.write(title+' | '+'ignore'+'\n')
            
            # Check if MeSH synonym is in the text
            for category, synonyms in category2permuted_synonyms.copy().items():
                found_syns = set()
                one_syn_in_title = False
                                                    
                # Each synonym in a set category
                for synonym in synonyms:

                    # Synonym in title
                    if ' '+synonym.lower()+' ' in title.lower():
                        one_syn_in_title = True
                        found_syns.add(synonym)
                        #continue
                        
                    # Synonym in abstract
                    
                    if synonym.lower() in abstract.lower():
                        add = True

                        # If similar synonym hasn't been counted already
                        for found_syn in found_syns:
                            fsyn = found_syn.lower().replace('\'','')
                            syn = synonym.lower()
                            if fsyn in syn or syn in fsyn:
                                add = False
                                break
                        if add:
                            found_syns.add(synonym)
                        
                    # Categorize text with 1+ synonym per category in the text
                    # (This could be modified to include confidence levels for
                    #  how many synonyms were found, remove break then)
                    if one_syn_in_title or len(found_syns) > 1:
                        #print(pmid, category, found_syns)
                        fout.write(pmid+'|'+category+'\n')
                        fout1.write(pmid+'|'+synonym+'\n')
                        break
                        
        try:
            print(str(num_pmids)+' PMIDs processed in one of the batches')
        except:
            # Not enough PMIDs to require doing all the batches
            pass
        

                
def multiprocess_ds_label_matching(pmids, index_name, index_type, 
                                   label_unlabeled_only,
                                   label_labeled_only,
                                   label_all,
                                   filter_list,
                                   the_function,
                                   stop_at_this_many_pmids,
                                   run = False,
                                   test = False,
                                   output_folder='../parsed_mappings/MeSH/'):
    '''
    FUNCTION:
    - This takes a list of strings and splits it into input
      for separate processes. The processes then output 
      their results to temp files which are then merged. 
    
    PARARMS:
    - pmids: the list to be split into input for 
      a multiprocessing function
    - the_function: the function that will use the list
      as input for multiprocessing
    '''    
    # How many processors can be used
    procs = cpu_count() if len(pmids) > cpu_count() else len(pmids)

    # List of batches for multiprocessing
    batches = [[] for i in range(procs)]   

    # Length of input list 
    tot = len(pmids)

    # Create batches and send to multiprocessing
    for i, item in enumerate(pmids):

        # Add synonym to a batch
        b_id = i%procs
        batches[b_id].append(item)

    # Create a list of jobs 
    print("Running jobs...")
    jobs = []
    for b_id, batch in enumerate(batches):
        jobs.append(Process(target = the_function,
                            args = [b_id, batch, 
                                    index_name, index_type, 
                                    label_unlabeled_only,
                                    label_labeled_only,
                                    label_all,
                                    filter_list,
                                    stop_at_this_many_pmids, 
                                    run, test, output_folder]))

    # Run the jobs
    for j in jobs: j.start()
    for j in jobs: j.join()
    print('Done!')    
            
                
def merge_pmid2new_mesh_labels(run=False, test=False, output_folder = '../parsed_mappings/MeSH/'):
    '''
    FUNCTION:
    - Merges the separate files containing PMID|category_name
      for the imputed category_name labels.
    '''
    
    pmid2imputed_meshsynonym, pmid2imputed_category = dict(), dict()
    procs = cpu_count()
    if run:
        rot = 'run'
    elif test:
        rot = 'test'

    for batch_id in range(procs):
        ''' PMID - MeSH Category '''
        temp_outfile = os.path.join(output_folder, rot+'_temp_labeling'+str(batch_id)+'.txt')
        with open(temp_outfile) as fin:
            for line in fin:
                line = line.split('|')

                # PMID, MeSH
                assert len(line) == 2
                pmid = line[0]
                mesh = line[1].strip()

                # PMID->MeSH
                pmid2imputed_category.setdefault(pmid, set()).add(mesh)
                
        ''' PMID - MeSH Synonym '''
        with open(temp_outfile[:-4]+'synonym'+'.txt') as fin1:
            for line in fin1:
                line = line.split('|')

                # PMID, MeSH
                assert len(line) == 2
                pmid = line[0]
                mesh = line[1].strip()

                # PMID->MeSH
                pmid2imputed_meshsynonym.setdefault(pmid, set()).add(mesh)                
                
    pmid2imputed_category = switch_dictset_to_dictlist(pmid2imputed_category)
    pmid2imputed_meshsynonym = switch_dictset_to_dictlist(pmid2imputed_meshsynonym)
    print(len(pmid2imputed_category), 'PMIDs with imputed labels')

    # Export PMID-Category mappings to dictionaries
    json.dump(pmid2imputed_category, 
              open(os.path.join(output_folder, 'pmid2imputed_category.json'), 'w'))
    json.dump(pmid2imputed_meshsynonym,
              open(os.path.join(output_folder, 'pmid2imputed_mesh_synonym.json'), 'w'))
    
    return pmid2imputed_category, pmid2imputed_meshsynonym


def get_groundtruth_pmid2categories(relevant_pmid_batch, index_name, 
                                    index_type, permuted_synonyms2category):
    '''
    FUNCTION:
    - This gets the ground truth labels, the MeSH-labeled PubMed documents
    
    PARAMS:
    - relevant_pmid_batch: This is the list of the PubMed IDs you want to 
      get ground truth for
    - index_name: name of the ElasticSearch index
    - index_type: name of the type of ElasticSearch index
    - permuted_synonyms2category: MeSH synonyms -> category
    '''
    pmids2real_categories = dict()
    es = Elasticsearch()
    total_pmids = len(relevant_pmid_batch)
    
    for num_pmids, pmid in enumerate(relevant_pmid_batch):
        if num_pmids % 10000 == 0:
            print('Getting real mappings' + str(num_pmids)+'/'+str(total_pmids),\
                  end='\r')
        entry = es.get(id = pmid, index = index_name, doc_type = index_type)
        pmids2real_categories[pmid] = set()
        
        # Publication's MeSH (if any)
        meshes = entry['_source']['MeSH'] 
        for mesh in meshes:
            try:
                categories = permuted_synonyms2category[mesh]
                for category in categories:
                    pmids2real_categories[pmid].add(category)
            except:
                pass
    pmids2real_categories = switch_dictset_to_dictlist(pmids2real_categories)
    return pmids2real_categories


def evaluate_label_imputation(pmid2imputed_category, pmid2real_categories):
    tp, fp, tn, fn = 0,0,0,0

    for pmid in pmid2real_categories:
        real_categories = pmid2real_categories[pmid]
        try: imputed_categories = pmid2imputed_category[pmid]
        except: imputed_categories = []
        real_and_imputed = real_categories + imputed_categories
        
        for category in real_and_imputed:
            #print(real_and_imputed)
            #print(real_categories)
            #print(imputed_categories)
            
            # Real = Yes
            if category in real_categories:
                
                # Real = Yes, Impute = Yes
                if category in imputed_categories:
                    tp += 1
                    
                # Real = Yes, Impute = No
                else:
                    fn += 1
                    
            elif category not in real_categories:
                
                # Real = No, Impute = Yes
                if category in imputed_categories:
                    fp += 1            

    print('Precision', round(tp/(tp+fp), 4))
    print('Recall', round(tp/(tp+fn), 4))
    print('TP', tp, 'FP', fp, 'FN',fn)
    
    
def index_imputed_mesh_categories(index_name, index_type):    
    '''
    FUNCTION:
    - index the imputed MeSH categories as MeSH terms
      in the ElasticSearch index
    
    PARAMS:
    - index_name: Name of the ElasticSearch index
    - index_type: Type name of the ElastichSearch index
    '''
    
    es = Elasticsearch()
    for pmid, imputed_category in pmid2imputed_category.items():
        
        # Get current MeSH terms
        entry = es.get(id = pmid, 
                       index = index_name, 
                       doc_type = index_type)
        mesh_terms = entry['_source']['MeSH']
        mesh_terms += imputed_category
        mesh_terms = list(set(mesh_terms))
        
        # Update each publication's index
        es.update(index = index_name, 
                  id = pmid,
                  doc_type = index_type,
                  doc = {'MeSH': mesh_terms})



def remove_imputed_category_mesh_terms(index_name, index_type, category_names):
    '''
    FUNCTION:
    - If you want to remove the imputed category names
      from the ElasticSearch index, undoing what you
      did with index_imputed_mesh_categories(), you
      can run this. Ensure that your imputed category 
      names aren't homonyms of actual MeSH terms.
    '''

    es = Elasticsearch()
    for num_pmids, entry in enumerate(es_iterate_all_documents(es, index_name)):
        if num_pmids % 10000 == 0:
            print(str(num_pmids)+' PMIDs checked for undoing label imputation', end='\r')
        
        # Get publication's MeSH terms
        mesh_terms = entry['MeSH']
        mesh_terms = list(set(mesh_terms))
        if len(mesh_terms) == 0:
            continue
        
        # Update MeSH term list
        updated_mesh = False
        for imputed_category_name in category_names:
            try: 
                mesh_terms.remove(imputed_category_name)
                updated_mesh = True
            except:
                continue

        # Update each publication's index
        if updated_mesh:
            pmid = entry['pmid']
            es.update(index = index_name,
                      id = pmid,
                      doc_type = index_type,
                      doc = {'MeSH': mesh_terms})
            
            
def remove_imputed_category_mesh_terms_previous_li(index_name,index_type, pmid2imputed_category):
    '''
    FUNCTION:
    - If you want to remove the imputed category names
      from the ElasticSearch index, undoing what you
      did in the last index_imputed_mesh_categories(), you
      can run this. Note: This does not undo label imputation
      from previous times you ran that command. This just undos
      the last time you did.
    '''

    es = Elasticsearch()
    for pmid, imputed_categories in pmid2imputed_category.items():

        # Get current MeSH terms
        entry = es.get(id=pmid,
                       index=index_name,
                       doc_type=index_type)
        mesh_terms = entry['_source']['MeSH']
        mesh_terms = list(set(mesh_terms))
        
        # Remove imputed category MeSH terms
        for imputed_category in imputed_categories:
            try:
                mesh_terms.remove(imputed_category)
            except:
                continue

        # Update each publication's index
        es.update(index=index_name,
                  id=pmid,
                  doc_type=index_type,
                  doc={'MeSH': mesh_terms})


def update_textcube_files(data_folder='../data/',
                          config_folder='../config/',
                          parsed_mapping_folder='../parsed_mappings'):
    '''
    FUNCTION:
    - Add category names to the considered MeSH Terms
    - Add pmid-category to pmid2category mapping files
    '''
    meshterms_per_cat = json.load(open(os.path.join(data_folder, 'meshterms_per_cat.json'), 'r'))
    meshterms_per_cat = [set(meshlist) for meshlist in meshterms_per_cat]

    category_names = json.load(open(os.path.join(config_folder, 'textcube_config.json'), 'r'))

    for i in range(0, len(category_names)):
        meshterms_per_cat[i].add(category_names[i])
    meshterms_per_cat = [list(meshlist) for meshlist in meshterms_per_cat]
    json.dump(meshterms_per_cat, open(os.path.join(data_folder, 'meshterms_per_cat.json'), 'w'))

    textcube_pmid2category = json.load(open(os.path.join(data_folder, 'textcube_pmid2category.json'), 'r'))
    textcube_category2pmid = json.load(open(os.path.join(data_folder, 'textcube_category2pmid.json'), 'r'))

    ''' Update textcube_category2pmid '''
    category_names = json.load(open(os.path.join(config_folder, 'textcube_config.json'), 'r'))
    category_name2num = {name: num for num, name in enumerate(category_names)}

    for pmid, imputed_categories in pmid2imputed_category.items():
        for imputed_category in imputed_categories:
            cat_num = category_name2num[imputed_category]
            textcube_category2pmid[cat_num].append(pmid)

    for i in range(0, len(textcube_category2pmid)):
        textcube_category2pmid[i] = list(set(textcube_category2pmid[i]))

    ''' Update textcube_pmid2category '''
    new_textcube_pmid2category = list()
    for cat_num, pmid_list in enumerate(textcube_category2pmid):
        for pmid in pmid_list:
            new_textcube_pmid2category.append([pmid, cat_num])

    json.dump(new_textcube_pmid2category, open(os.path.join(data_folder, 'textcube_pmid2category.json'), 'w'))
    json.dump(textcube_category2pmid, open(os.path.join(data_folder, 'textcube_category2pmid.json'), 'w'))


def get_relevant_all_categorized_pmids(index_name):
    '''
    FUNCTION:
    - Get the PMIDs of all the PubMed articles that
      don't have MeSH labels.
    - Purpose: For applying to unlabeled PubMed articles
    
    PARAMS:
    - index_name: The name of the ElasticSearch index
      where all the PubMed articles are indexed/stored
    '''
    
    es = Elasticsearch()
    relevant_categorized_pmids = set()
    for num_pmids, entry in enumerate(es_iterate_all_documents(es, index_name)):

        # Publication's MeSH (if any)
        meshes = entry['MeSH']
        labeled = len(meshes) > 0
        if not labeled:
            continue
        else:
            pmid = entry['pmid']
            relevant_categorized_pmids.add(pmid)
            
        if num_pmids % 1000 == 0:
            print(str(num_pmids)+' PMIDs processed',end='\r')
    relevant_categorized_pmids = list(relevant_categorized_pmids) 
    
    return relevant_categorized_pmids










##############
#### MAIN ####
##############
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--index_name', default='pubmed_lift', type=str)
parser.add_argument('-it', '--index_type', default='pubmed_meta_lift', type=str)
parser.add_argument('-li', '--label_imputation', default=True, type=bool)
parser.add_argument('-u', '--undo_label_imputation', default=False, type=bool)
parser.add_argument('-ul', '--undo_last_label_imputation', default=False, type=bool)
parser.add_argument('-r', '--run', default=True, type=bool)
parser.add_argument('-t', '--test', default=False, type=bool)
parser.add_argument('-p', '--consider_this_many_pmids', default=9999999999, type=int)
args = parser.parse_args()

index_name = args.index_name
index_type = args.index_type
undo_category_label_imputation = args.undo_label_imputation
undo_last_category_label_imputation_only = args.undo_last_label_imputation
label_imputation = args.label_imputation
run = args.run
test = args.test

if undo_category_label_imputation:
    label_imputation = False

if test == True:
    run = False
STOP_AT_THIS_MANY_PMIDS = int(args.consider_this_many_pmids)
filter_words = ['placeholder', 'these arent used anymore']
#['heart', 'cardiac', 'cardiovascular', 'cardiopulmonary']
              
              
              
root_directory = '/caseolap_lift_shared_folder'
root_directory = '../'
data_folder=os.path.join(root_directory,"data")
config_folder=os.path.join(root_directory,"config")
mapping_folder=os.path.join(root_directory,"parsed_mappings")
output_folder = os.path.join(mapping_folder,'MeSH')
    
    
    
    
    
    
    
    
    
    
    
    
''' Label imputation '''
if label_imputation:

    '''Map MeSH ID - Name'''
    #map_disease_mesh_name_to_id()


    meshterms_per_cat_file = os.path.join(data_folder,'meshterms_per_cat.json')
    category_names_file = os.path.join(config_folder,'textcube_config.json')
              
    ''' Get MeSH Synonyms '''
    # Category - MeSH Terms
    category2main_terms = map_categories2terms(output_folder = output_folder,
                            meshterms_per_cat_file = meshterms_per_cat_file,
                            category_names_file = category_names_file)

    # Category - MeSH Terms' Synonyms (including terms)
    name2id = json.load(open(os.path.join(mapping_folder,'MeSH/meshterm-IS-meshid.json'),'r'))
    category2main_terms = json.load(open(os.path.join(mapping_folder, 'MeSH/category2terms.json')))
    try:  category2synonym_terms = json.load(open(os.path.join(mapping_folder,'MeSH/category2synonyms.json')))
    except: category2synonym_terms = get_mesh_synonyms_api(category2main_terms, name2id,output_folder=output_folder)
              
    # Category - permuted MeSH Terms' Synonyms
    temp1, temp2 =  map_category_to_permuted_mesh_synonyms(category2synonym_terms,output_folder=output_folder)
    category2permuted_synonyms = temp1
    permuted_synonyms2category = temp2

    
    
    
    
    
    
    
    
    
    if run:
        print('Running label imputation')
        
        '''Get relevant PMIDs'''
        #try:
        #    relevant_pmids = json.load(open(data_folder+'all_uncategorized_pmids.json'))
        #except:
        relevant_pmids = get_all_uncategorized_pmids(index_name)
        json.dump(relevant_pmids, open(os.path.join(data_folder,'all_uncategorized_pmids.json'),'w'))
    
        print(len(relevant_pmids), 'relevant pmids')

            
        ''' Impute missing MeSH labels '''
        # Impute PMIDs' Categories (i.e., Impute missing MeSH labels)
        multiprocess_ds_label_matching(pmids = relevant_pmids, 
                                       index_name = index_name, 
                                       index_type = index_type, 
                                       label_unlabeled_only = True,
                                       label_labeled_only = False,
                                       label_all = False,
                                       filter_list = filter_words,
                                       stop_at_this_many_pmids = STOP_AT_THIS_MANY_PMIDS,
                                       the_function = ds_label_matching,
                                       output_folder = output_folder,
                                       run = True)
        pmid2imputed_category, pmid2imputed_meshsynonym = merge_pmid2new_mesh_labels(run=True,output_folder=output_folder)


        # Index the imputed MeSH categories into their PMID entries
        index_imputed_mesh_categories(index_name=index_name, index_type=index_type)

        # Update the MeSH Terms Per Category file for the textcube
        update_textcube_files(data_folder=data_folder, config_folder=config_folder, parsed_mapping_folder=mapping_folder)
        
        
        
        
        
        
        
        
        
        
    if test:
        print('Testing label imputation')
        
        '''Get relevant PMIDs'''
        # Ground truth: PMIDs in your study already labeled by NIH with MeSH terms
        #relevant_pmids = json.load(open('../caseolap/data/pmids.json'))
        try:
            relevant_pmids = json.load(open(data_folder+'/relevant_categorized_pmids.json'))
        except:
            relevant_pmids_path = data_folder+'/relevant_categorized_pmids.json'
            relevant_pmids = get_all_categorized_pmids('pubmed_lift', relevant_pmids_path)
            json.dump(relevant_pmids, open(relevant_pmids_path, 'w'))

        ''' Impute missing MeSH labels '''
        # Impute PMIDs' Categories (i.e., Impute missing MeSH labels)
        multiprocess_ds_label_matching(pmids = relevant_pmids, 
                                       index_name = index_name, 
                                       index_type = index_type, 
                                       label_unlabeled_only = False,
                                       label_labeled_only = True,
                                       label_all = False,
                                       filter_list = filter_words,
                                       stop_at_this_many_pmids = STOP_AT_THIS_MANY_PMIDS,
                                       the_function = ds_label_matching,
                                       test = True)
        pmid2imputed_category, pmid2imputed_meshsynonym = merge_pmid2new_mesh_labels(test=True)

        # PMID - Ground Truth Categories
        try:
            pmid2real_categories = json.load(open(output_folder+'/pmid2real_categories.json'))
        except:
            pmid2real_categories = ['']
        if len(pmid2real_categories) != STOP_AT_THIS_MANY_PMIDS:
            pmid2real_categories = get_groundtruth_pmid2categories(relevant_pmid_batch, index_name, index_type, permuted_synonyms2category)
            json.dump(pmid2real_categories, open(output_folder+'/pmid2real_categories.json','w'))

        # Evaluate results on ground truth
        evaluate_label_imputation(pmid2imputed_category, pmid2real_categories)
        
              
        
        
    
''' Undo label imputation '''
# Undo all LI attemps

config_file = os.path.join(root_directory,'config/textcube_config.json')
pmid2imputed_category_file = os.path.join(root_directory,'data/pmid2imputed_category.json')
if undo_category_label_imputation:
    try:
        category_names = json.load(open(config_file,'r'))
        print('Removing imputed labels:', category_names)
        remove_imputed_category_mesh_terms(index_name, index_type, category_names)
    except:
        raise Exception('Error undoing labels')


# Undo last LI attempt (only works for one time. Repeating this will not undo more)
if undo_last_category_label_imputation_only:
    try:
        pmid2imputed_category = json.load(open(pmid2imputed_category_file))
        remove_imputed_category_mesh_terms_previous_li(index_name, index_type, pmid2imputed_category)
    except:
        raise Exception('Check that you have performed label imputation and thus have the'+\
                        'pmid2imputed_category/parsed_mappings/MeSH/pmid2imputed_category.json file')
