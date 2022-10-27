import sys
# setting path
sys.path.append('..')

# importing
import requests as req
import xml.etree.ElementTree as ET
import itertools
from utils.biomedkg_utils import *
from elasticsearch import Elasticsearch


def es_iterate_all_documents(es, index, pagesize=250, scroll_timeout="1m", **kwargs):
    """
    Helper to iterate ALL values from a single index
    Yields all the documents.
    Source: https://techoverflow.net/2019/05/07/elasticsearch-how-to-iterate-scroll-through-all-documents-in-index/
    """
    is_first = True
    while True:
        # Scroll next
        if is_first:  # Initialize scroll
            result = es.search(index=index, scroll="1m", **kwargs,
                               size=pagesize)
            is_first = False
        else:
            result = es.scroll(
                scroll_id=scroll_id,
                scroll=scroll_timeout)
        scroll_id = result["_scroll_id"]
        hits = result["hits"]["hits"]

        # Stop after no more docs
        if not hits:
            break

        # Yield each entry
        yield from (hit['_source'] for hit in hits)


def map_disease_mesh_name_to_id():
    '''
    FUNCTION:
    - Map the disease MeSH names/terms to the
      MeSH IDs.
    '''
    # ! wget - N - P
    # input / https: // nlmpubs.nlm.nih.gov / projects / mesh / MESH_FILES / xmlmesh / desc2022.xml
    tree = ET.parse('../data/desc2022.xml')
    root = tree.getroot()

    name2id = dict()

    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')

            # If Tree is a disease
            for tree_number in tree_numbers:
                if tree_number.text.startswith(('C', 'F03')):
                    tree_number = tree_number.text

                    # ID to Name
                    try:
                        ID = ele.find('DescriptorUI').text
                        name = ele.find('DescriptorName').find('String').text
                        name2id.setdefault(name, set()).add(ID)
                    except:
                        pass
        except:
            continue

    name2id = switch_dictset_to_dictlist(name2id)
    json.dump(name2id, open('../data/name2id.json', 'w'))


def map_categories2terms():
    '''
    FUNCTION:
    - Map disease categories -> MeSH terms
    '''
    terms_lol = json.load(open('../data/meshterms_per_cat.json')) #TODO where is this from?
    category_names = json.load(open('../data/textcube_config.json'))#TODO this is temporary
    category2terms = dict()
    for category_name, term_list in zip(category_names, terms_lol):
        category2terms[category_name] = term_list

    return category2terms


def get_mesh_synonyms_api(category2terms, name2id):
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
            r = req.get('https://id.nlm.nih.gov/mesh/lookup/details?descriptor=' + mesh_id[0]).json()

            # Add term synonyms to the category
            for entry in r['terms']:
                synonym = entry['label']
                category2synonyms.setdefault(category, set()).add(synonym)

    return category2synonyms


def map_category_to_permuted_mesh_synonyms(category2synonyms):
    '''
    FUNCTION:
    - permute the MeSH Synonyms because sometimes they're
      written weirdly with commas (e.g., Order, Messed, Up, Is;
      Messed, Order, Up)

    PARAMS:
    - category2synonyms: category -> MeSH synonyms mapping
    '''

    # Permuted synonym->category dictionary
    category2permuted_synonyms, permuted_synonyms2category = dict(), dict()
    for category, synonyms in category2synonyms.items():
        synonyms = permute_mesh_synonyms(synonyms)
        for synonym in synonyms:
            category2permuted_synonyms.setdefault(category, set()).add(synonym)
            permuted_synonyms2category.setdefault(synonym, set()).add(category)

    category2permuted_synonyms = switch_dictset_to_dictlist(category2permuted_synonyms)
    permuted_synonyms2category = switch_dictset_to_dictlist(permuted_synonyms2category)

    json.dump(category2permuted_synonyms, open('../data/category2permuted_synonyms.json', 'w'))
    json.dump(permuted_synonyms2category, open('../data/permuted_synonyms2category.json', 'w'))

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


def get_relevant_categorized_pmids():
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
    list_of_list_of_pmids = json.load(open('../data/textcube_category2pmid.json'))
    relevant_pmids = set()
    for pmid_list in list_of_list_of_pmids:
        relevant_pmids = relevant_pmids.union(set(pmid_list))
    relevant_pmids = list(relevant_pmids)

    return relevant_pmids


def get_relevant_uncategorized_pmids(index_name):
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
        meshes = entry['MeSH']
        labeled = len(meshes) > 0
        if labeled:
            continue
        else:
            pmid = entry['pmid']
            relevant_uncategorized_pmids.add(pmid)

        if num_pmids % 1000 == 0:
            print(str(num_pmids) + ' PMIDs processed', end='\r')
    relevant_uncategorized_pmids = list(relevant_uncategorized_pmids)

    return relevant_uncategorized_pmids


def get_document_text(entry, full_text=False):
    '''
    Function: Get the full text of the PubMed publication
    Params: entry = parsed indexed publication
    - full_text (bool): indicate whether full text should be used in label imputation / matching
    '''
    document_text = ''
    text_sections = [entry['_source']['title'], entry['_source']['abstract']]
    if full_text == True:
        text_sections.append(entry['_source']['full_text'])

    for text_section in text_sections:
        if type(text_section) == str:
            document_text += text_section.replace('\n', ' ').replace('\t', ' ').replace('   ', ' ') + ' '

    return document_text


def ds_label_matching(batch_id, relevant_pmid_batch,
                      index_name, index_type,
                      label_unlabeled_only,
                      label_labeled_only,
                      label_all,
                      filter_list,
                      stop_at_this_many_pmids):
    '''
    FUNCTION:
    - Using MeSH term synonyms, label a document with a MeSH term if the
      lowercased term matches exactly within the text.

    PARAMS:
    - relevant_pmid_batch: PMIDs to look for in this batch
    - temp_outfile: temporary output file of PMID | MeSH Term
    '''
    es = Elasticsearch()
    temp_outfile = 'data/temp_labeling' + str(batch_id) + '.txt'
    procs = cpu_count()
    category2permuted_synonyms = json.load(open('data/category2permuted_synonyms.json'))

    with open(temp_outfile, 'w') as fout, open(temp_outfile[:-4] + 'synonym' + '.txt', 'w') as fout1:

        # Get PubMed article's text
        for num_pmids, pmid in enumerate(relevant_pmid_batch):
            entry = es.get(id=pmid, index=index_name, doc_type=index_type)

            # Print progress and break early
            if batch_id == 1 and num_pmids % 1000 == 0:
                print(str(num_pmids) + ' PMIDs processed in one of the batches', end='\r')
            if num_pmids > stop_at_this_many_pmids:
                break

            # Determine which publications to find labels for
            true_meshes = entry['_source']['MeSH']
            labeled = len(true_meshes) > 0
            unlabeled = len(true_meshes) == 0
            if label_labeled_only:
                if labeled:
                    pass
                else:
                    continue
            elif label_unlabeled_only:
                if unlabeled:
                    pass
                else:
                    continue
            elif label_all:
                pass

            # Publication's text (title, abstract, full text if provided)
            document_text = get_document_text(entry, full_text=False)

            # Only consider documents that have certain broader key words
            # E.g., for cardiovascular disease publications, only look
            # at publications that say "heart" or "cardiac"
            dont_label = True
            for filter_word in filter_words:
                if filter_word in document_text.lower():
                    dont_label = False
                    break
            if dont_label:
                continue

            # Check if MeSH synonym is in the text
            for category, synonyms in category2permuted_synonyms.items():
                found_syns = set()

                for synonym in synonyms:
                    if synonym.lower() in document_text.lower():

                        # Check if a similar synonym was found
                        if len(found_syns) == 0:
                            found_syns.add(synonym.lower())
                        else:
                            for syn in found_syns.copy():
                                if syn in synonym.lower() or synonym.lower() in syn:
                                    continue
                                else:
                                    found_syns.add(synonym.lower())

                        # If a category has 2+ synonyms in the text, label as the category
                        if len(found_syns) >= 2:
                            fout.write(pmid + '|' + category + '\n')
                            fout1.write(pmid + '|' + synonym + '\n')


def multiprocess_ds_label_matching(pmids, index_name, index_type,
                                   label_unlabeled_only,
                                   label_labeled_only,
                                   label_all,
                                   filter_list,
                                   the_function,
                                   stop_at_this_many_pmids):
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
    thelist = pmids

    # How many processors can be used
    procs = cpu_count()

    # List of batches for multiprocessing
    batches = [[] for i in range(procs)]

    # Length of input dictionary
    tot = len(thelist)

    # Create batches and send to multiprocessing
    for i, item in enumerate(thelist):
        # Add synonym to a batch
        b_id = i % procs
        batches[b_id].append(item)

    # Create a list of jobs
    print("Running jobs...")
    jobs = []
    for b_id, batch in enumerate(batches):
        jobs.append(Process(target=the_function, \
                            args=[b_id, batch,
                                  index_name, index_type,
                                  label_unlabeled_only,
                                  label_labeled_only,
                                  label_all,
                                  filter_list,
                                  stop_at_this_many_pmids]))

    # Run the jobs
    for j in jobs: j.start()
    for j in jobs: j.join()
    print('Done!')


def merge_pmid2new_mesh_labels():
    '''
    FUNCTION:
    - Merges the separate files containing PMID|category_name
      for the imputed category_name labels.
    '''

    pmid2imputed_mesh_synonym, pmid2imputed_category = dict(), dict()
    procs = cpu_count()

    for batch_id in range(procs):
        temp_outfile = 'data/temp_labeling' + str(batch_id) + '.txt'
        with open(temp_outfile) as fin:
            for line in fin:
                line = line.split('|')

                # PMID, MeSH
                assert len(line) == 2
                pmid = line[0]
                mesh = line[1].strip()

                # PMID->MeSH
                pmid2imputed_category.setdefault(pmid, set()).add(mesh)

        with open(temp_outfile[:-4] + 'synonym' + '.txt') as fin1:
            for line in fin1:
                line = line.split('|')

                # PMID, MeSH
                assert len(line) == 2
                pmid = line[0]
                mesh = line[1].strip()

                # PMID->MeSH
                pmid2imputed_mesh_synonym.setdefault(pmid, set()).add(mesh)

    pmid2imputed_category = switch_dictset_to_dictlist(pmid2imputed_category)
    pmid2imputed_mesh_synonym = switch_dictset_to_dictlist(pmid2imputed_mesh_synonym)

    return pmid2imputed_category, pmid2imputed_mesh_synonym


def get_groundtruth_pmid2categories(relevant_pmid_batch, index_name, index_type, permuted_synonyms2category):
    '''
    FUNCTION:
    - This gets the ground truth labels, the MeSH-labeled PubMed documents

    PARAMS:
    - relevant_pmid_batch: This is the list of the PubMed IDs you want to get ground truth for
    - index_name: name of the ElasticSearch index
    - index_type: name of the type of ElasticSearch index
    - permuted_synonyms2category: MeSH synonyms -> category
    '''
    pmids2real_categories = dict()
    es = Elasticsearch()
    total_pmids = len(relevant_pmid_batch)

    for num_pmids, pmid in enumerate(relevant_pmid_batch):
        print('Getting real mappings' + str(num_pmids) + '/' + str(total_pmids), end='\r')
        entry = es.get(id=pmid, index=index_name, doc_type=index_type)
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
    tp, fp, tn, fn = 0, 0, 0, 0

    for pmid in pmid2real_categories:
        real_categories = pmid2real_categories[pmid]
        try:
            imputed_categories = pmid2imputed_category[pmid]
        except:
            imputed_categories = []
        real_and_imputed = real_categories + imputed_categories

        for category in real_and_imputed:
            # print(real_and_imputed)
            # print(real_categories)
            # print(imputed_categories)

            # Real = Yes
            if category in real_categories:

                # Real = Yes, Impute = Yes
                if category in imputed_categories:
                    tp += 1

                # Real = Yes, Impute = No
                else:
                    fn += 1

            # Real = No
            elif category not in real_categories:

                # Real = No, Impute = Yes
                if category in imputed_categories:
                    fp += 1

    print('Precision', round(tp / (tp + fp), 4))
    print('Recall', round(tp / (tp + fn), 4))
    print('TP', tp, 'FP', fp, 'FN', fn)


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
        entry = es.get(id=pmid,
                       index=index_name,
                       doc_type=index_type)
        mesh_terms = entry['_source']['MeSH']
        mesh_terms += imputed_category
        mesh_terms = list(set(mesh_terms))

        # Update each publication's index
        es.update(index=index_name,
                  id=pmid,
                  doc_type=index_type,
                  doc={'MeSH': mesh_terms})


def remove_imputed_category_mesh_terms(index_name,index_type,pmid2imputed_category):
    '''
    FUNCTION:
    - If you want to remove the imputed category names
      from the ElasticSearch index, undoing what you
      did with index_imputed_mesh_categories(), you
      can run this
    '''

    es = Elasticsearch()
    for pmid, imputed_categories in pmid2imputed_category.items():

        # Get current MeSH terms
        entry = es.get(id=pmid,
                       index=index_name,
                       doc_type=index_type)
        mesh_terms = entry['_source']['MeSH']
        mesh_terms = list(set(mesh_terms))
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


def update_textcube_files():
    '''
    FUNCTION:
    - Add category names to the considered MeSH Terms
    - Add pmid-category to pmid2category mapping files
    '''
    meshterms_per_cat = json.load(open('../data/meshterms_per_cat.json'))
    meshterms_per_cat = [set(meshlist) for meshlist in meshterms_per_cat]

    category_names = json.load(open('config/textcube_config.json'))

    for i in range(0, len(category_names)):
        meshterms_per_cat[i].add(category_names[i])
    meshterms_per_cat = [list(meshlist) for meshlist in meshterms_per_cat]
    json.dump(meshterms_per_cat, open('../data/meshterms_per_cat.json', 'w'))

    textcube_pmid2category = json.load(open('../data/textcube_pmid2category.json'))
    textcube_category2pmid = json.load(open('../data/textcube_category2pmid.json'))

    ''' Update textcube_category2pmid '''
    category_names = json.load(open('config/textcube_config.json'))
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

    json.dump(new_textcube_pmid2category, open('data/textcube_pmid2category.json', 'w'))
    json.dump(textcube_category2pmid, open('data/textcube_category2pmid.json', 'w'))

'''Map MeSH ID - Name'''
map_disease_mesh_name_to_id()

''' Get MeSH Synonyms '''
# Category - MeSH Terms
category2terms = map_categories2terms()

# Category - MeSH Terms' Synonyms (including terms)
name2id = json.load(open('../data/name2id.json'))
category2synonyms = get_mesh_synonyms_api(category2terms, name2id)

# Category - permuted MeSH Terms' Synonyms
category2permuted_synonyms, permuted_synonyms2category =  map_category_to_permuted_mesh_synonyms(category2synonyms)

relevant_categorized_pmids = get_relevant_uncategorized_pmids(index_name = 'pubmed')

'''Load data'''
category2permuted_synonyms = json.load(open('../data/category2permuted_synonyms.json'))
#relevant_pmids = get_relevant_categorized_pmids()
relevant_uncategorized_pmids = get_relevant_uncategorized_pmids(index_name = 'pubmed')
json.dump(list(relevant_uncategorized_pmids), open('../data/relevant_uncategorized_pmids.json','w'))
relevant_uncategorized_pmids = json.load(open('../data/relevant_uncategorized_pmids.json'))


''' Impute missing MeSH labels '''
STOP_AT_THIS_MANY_PMIDS = 9999999999
filter_words = ['heart', 'cardiac', 'cardiovascular', 'cardiopulmonary']

# Impute PMIDs' Categories (i.e., Impute missing MeSH labels)
multiprocess_ds_label_matching(pmids = relevant_categorized_pmids,
                               index_name = 'pubmed',
                               index_type = 'pubmed_meta',
                               label_unlabeled_only = True,
                               label_labeled_only = False,
                               label_all = False,
                               filter_list = filter_words,
                               stop_at_this_many_pmids = STOP_AT_THIS_MANY_PMIDS,
                               the_function = ds_label_matching)
pmid2imputed_category, pmid2imputed_mesh_synonym = merge_pmid2new_mesh_labels()


# Export PMID-Category mappings to dictionaries
json.dump(pmid2imputed_category, open('../data/pmid2imputed_category.json','w'))
json.dump(pmid2imputed_mesh_synonym, open('../data/pmid2imputed_mesh_synonym.json','w'))

'''Update existing external files with imputed labels: index, textcube_config'''
# Index the imputed MeSH categories into their PMID entries
index_imputed_mesh_categories(index_name='pubmed', index_type='pubmed_meta')
print(len(pmid2imputed_category), 'PMIDs with imputed labels')

# Optional: Remove the imputed mesh terms from the index
#remove_imputed_category_mesh_terms(pmid2imputed_category)

# Update the files for the textcube based on the new imputed labels
update_textcube_files()