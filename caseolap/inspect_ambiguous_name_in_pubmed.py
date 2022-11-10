'''
When searching for proteins in the PubMed articles, you will match protein names with the article text.
The protein names may be non-specific, referring to things other than the protein. For example, CAD may refer to
a protein or it may refer to a disease, Coronary Artery Disease. To spot check how frequently a word or phrase refers
to the protein you have in mind or to something else, you can run this script by passing in a protein name. For example,
run "python inspect_ambiguous_name_in_pubmed.py --n CAD" on the command line. 
'''

import argparse
import json
from elasticsearch import Elasticsearch
from elasticsearch_dsl import Search, Q


relevant_pmids = json.load(open('data/textcube_relevant_pmids.json'))

# Parser
parser = argparse.ArgumentParser()
parser.add_argument('-n', '--name', type=str)
parser.add_argument('-i', '--index_name', default='pubmed_lift', type=str)
parser.add_argument('-s', '--section', default='abstract', type=str)
args = parser.parse_args()
key = args.section
syn = args.name
index_name = args.index_name

es = Elasticsearch()

# For each synonym from a list of all synonyms
# Define the query
if key == "abstract":
    query = Q("match_phrase", title=syn) |\
            Q("match_phrase", abstract=syn)
elif key == "full_text":
    query = Q("match_phrase", abstract=syn) | \
            Q("match_phrase", title=syn) | \
            Q("match_phrase", full_text=syn)
else:
    print("ERROR: Key is neither 'abstract' nor 'full_text'")
    sys.exit

# Perform the query
s = Search(using=es, index=index_name).params(request_timeout=300).query(query) 

# For each synonym-containing publication
print(syn)
num_print = 0
for num_hits, hit in enumerate(s.scan()):                             

    # Proceed if the publication is of interest
    pmid = str(hit.pmid)
    if pmid in relevant_pmids:
        print('='*20,'\n',hit.title, hit.abstract, '\n')
        num_print += 1

        if num_print > 10:
            break