import json, csv, pandas as pd, numpy as np, requests, subprocess, os,  urllib.parse, urllib.request, html
from bs4 import BeautifulSoup

import re, time, json, zlib
from xml.etree import ElementTree
from urllib.parse import urlparse, parse_qs, urlencode
import requests
from requests.adapters import HTTPAdapter, Retry


POLLING_INTERVAL = 3

API_URL = "https://rest.uniprot.org"


retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
session_UniProtAPI = requests.Session()
session_UniProtAPI.mount("https://", HTTPAdapter(max_retries=retries))


def submit_id_mapping_UniProtAPI(from_db, to_db, ids):
    '''
    FUNCITON:
    - This submits the post request to UniProt's API.
    
    PARAMS:
    - from_db (string): The database to map IDs from
    - to_db (string): The database to map IDs to
    - ids (list of strings): The IDs to map from 
    '''
    request = requests.post(
        f"{API_URL}/idmapping/run",
        data={"from": from_db, "to": to_db, "ids": ",".join(ids)},
    )
    request.raise_for_status()
    return request.json()["jobId"]


def get_next_link_UniProtAPI(headers):
    re_next_link = re.compile(r'<(.+)>; rel="next"')
    if "Link" in headers:
        match = re_next_link.match(headers["Link"])
        if match:
            return match.group(1)


def check_id_mapping_results_ready_UniProtAPI(job_id):
    '''
    FUNCTION:
    - This checks if the submitted job is ready.
      If the job is not ready, it will tell you and try again.
      If the job is ready, it will tell you it is ready.
      
    PARAMS:
    - job_id: This is the ID of the job submitted to UniProt
    '''
    while True:
        request = session_UniProtAPI.get(f"{API_URL}/idmapping/status/{job_id}")
        request.raise_for_status()
        j = request.json()
        if "jobStatus" in j:
            if j["jobStatus"] == "RUNNING":
                print(f"Job still running. Retrying in {POLLING_INTERVAL}s")
                time.sleep(POLLING_INTERVAL)
            else:
                raise Exception(request["jobStatus"])
        else:
            return bool(j["results"] or j["failedIds"])


def get_batch_UniProtAPI(batch_response, file_format, compressed):
    batch_url = get_next_link_UniProtAPI(batch_response.headers)
    while batch_url:
        batch_response = session_UniProtAPI.get(batch_url)
        batch_response.raise_for_status()
        yield decode_results_UniProtAPI(batch_response, file_format, compressed)
        batch_url = get_next_link_UniProtAPI(batch_response.headers)


def combine_batches_UniProtAPI(all_results, batch_results, file_format):
    if file_format == "json":
        for key in ("results", "failedIds"):
            if key in batch_results and batch_results[key]:
                all_results[key] += batch_results[key]
    elif file_format == "tsv":
        return all_results + batch_results[1:]
    else:
        return all_results + batch_results
    return all_results


def get_id_mapping_results_link_UniProtAPI(job_id):
    '''
    FUNCTION:
    - This gets the link where the job results can
      be accessed and then later downloaded by another
      function here.
    
    PARAMS:
    - job_id: This is the ID of the job submitted to UniProt
    '''
    url = f"{API_URL}/idmapping/details/{job_id}"
    request = session_UniProtAPI.get(url)
    request.raise_for_status()
    return request.json()["redirectURL"]


def decode_results_UniProtAPI(response, file_format, compressed):
    if compressed:
        decompressed = zlib.decompress(response.content, 16 + zlib.MAX_WBITS)
        if file_format == "json":
            j = json.loads(decompressed.decode("utf-8"))
            return j
        elif file_format == "tsv":
            return [line for line in decompressed.decode("utf-8").split("\n") if line]
        elif file_format == "xlsx":
            return [decompressed]
        elif file_format == "xml":
            return [decompressed.decode("utf-8")]
        else:
            return decompressed.decode("utf-8")
    elif file_format == "json":
        return response.json()
    elif file_format == "tsv":
        return [line for line in response.text.split("\n") if line]
    elif file_format == "xlsx":
        return [response.content]
    elif file_format == "xml":
        return [response.text]
    return response.text


def get_xml_namespace_UniProtAPI(element):
    m = re.match(r"\{(.*)\}", element.tag)
    return m.groups()[0] if m else ""


def merge_xml_results_UniProtAPI(xml_results):
    merged_root = ElementTree.fromstring(xml_results[0])
    for result in xml_results[1:]:
        root = ElementTree.fromstring(result)
        for child in root.findall("{http://uniprot.org/uniprot}entry"):
            merged_root.insert(-1, child)
    ElementTree.register_namespace("", get_xml_namespace_UniProtAPI(merged_root[0]))
    return ElementTree.tostring(merged_root, encoding="utf-8", xml_declaration=True)


def print_progress_batches_UniProtAPI(batch_index, size, total):
    n_fetched = min((batch_index + 1) * size, total)
    print(f"Fetched: {n_fetched} / {total}", end='\r')


def get_id_mapping_results_search_UniProtAPI(url):
    '''
    FUNCTION:
    - Download the API results from a url
    
    PARAMS:
    - url: the link where the job results can
      be accessed and downloaded here.
    '''
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    if "size" in query:
        size = int(query["size"][0])
    else:
        size = 500
        query["size"] = size
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    parsed = parsed._replace(query=urlencode(query, doseq=True))
    url = parsed.geturl()
    request = session_UniProtAPI.get(url)
    request.raise_for_status()
    results = decode_results_UniProtAPI(request, file_format, compressed)
    total = int(request.headers["x-total-results"])
    print_progress_batches_UniProtAPI(0, size, total)
    for i, batch in enumerate(get_batch_UniProtAPI(request, file_format, compressed), 1):
        results = combine_batches_UniProtAPI(results, batch, file_format)
        print_progress_batches_UniProtAPI(i, size, total)
    if file_format == "xml":
        return merge_xml_results_UniProtAPI(results)
    return results


def get_id_mapping_results_stream_UniProtAPI(url):
    if "/stream/" not in url:
        url = url.replace("/results/", "/stream/")
    request = session_UniProtAPI.get(url)
    request.raise_for_status()
    parsed = urlparse(url)
    query = parse_qs(parsed.query)
    file_format = query["format"][0] if "format" in query else "json"
    compressed = (
        query["compressed"][0].lower() == "true" if "compressed" in query else False
    )
    return decode_results_UniProtAPI(request, file_format, compressed)


def get_possible_fields_UniProtAPI():
    '''
    FUNCTION:
    - This prints the possible fields including the
      databases that can be mapped from and to.
    '''
    # Raw results
    all_fields = requests.get('https://rest.uniprot.org/configure/idmapping/fields').json()

    # Database Category to Database dictionary
    database_dict = dict()

    for group in all_fields['groups']:

        # Database category
        database_category = group['groupName']

        # Databases
        databases = list()
        for item in group['items']:
            databases.append(item['name'])

        # Database category to Database
        database_dict[database_category] = databases
    

    return all_fields, database_dict



def get_to_uniprotid_from_genename_mapping_dict_UniProtAPI(results, your_taxa = [''], filter_by_reviewed = True):
    '''
    FUNCTION:
    - Convert the raw API call's mapping results to a dictionary
      where the keys are the "from" IDs and the values are the
      "to" IDs
    
    PARAMS:
    - results (dict): The UniProt API's raw results from an API call
    - your_taxa (list of ints): a list of taxa you want to include
    - filter_by_reviewed (bool): indicates if you only want reviewed entries
    '''
    name2id_up = dict()
    filter_by_taxa = your_taxa != ['']
    
    for i in range(0, len(results['results'])):
        # API mapping results for this entity
        r = results['results'][i]

        # Gene ID
        from_gene_name = r['from']
        try: 
            gene_name = r['to']['genes'][0]['geneName']['value']
        except: 
            continue
        if gene_name != from_gene_name:
            continue
            
        # Protein ID
        to_protein_id = r['to']['primaryAccession']

        
        # Organism taxon number
        taxon_num = r['to']['organism']['taxonId']

        # Is the entry reviewed
        reviewed = ' reviewed ' in r['to']['entryType']

        # If you're not filtering by 'reviewed' or if the entry is reviewed
        if not filter_by_reviewed or reviewed == True:

            # If you're not filtering by taxa or if the taxa is your taxa
            if not filter_by_taxa or taxon_num in your_taxa:
                name2id_up.setdefault(from_gene_name, set()).add(to_protein_id)
                
    return name2id_up





# Source for STRING API code: Alexander Pelletier
def get_filtered_string_interactors_STRINGAPI(string_ids, score_thresh):
    '''
    Get interactors for STRING IDs using the API, 
    only return those with combined_score >= score_thresh
    '''

    ### PPI API 
    string_api_url = "https://version-11-5.string-db.org/api"
    output_format = "tsv-no-header"
    method = "interaction_partners"
    request_url = "/".join([string_api_url, output_format, method])
    params = {
      "identifiers" : "%0d".join(string_ids), # your proteins
      "species" : 9606,                       # 9606 = humans 
      "required_score" : score_thresh*1000, # only report those with a score above 0.85
      "caller_identity" : "www.awesome_app.org" # your app name
    }
    response = requests.post(request_url, data=params)

    ### Save PPI results
    results = []
    for line in response.text.strip().split("\n"):
        l = line.strip().split("\t")

        query_ensp = l[0]
        query_name = l[2]
        partner_ensp = l[1]
        partner_name = l[3]
        combined_score = float(l[5])
        
        if combined_score >= score_thresh:
            results += [[query_ensp, query_name, partner_ensp, partner_name, combined_score]]
    return results



def batch_filtered_queries_STRINGAPI(string_ids, batch_size = 100, score_thresh = 0.9):
    '''
    API crashes if you put the whole list. batch the queries
    '''
    string_ids = list(string_ids)
    batched_results = []
    for i in range(0,len(string_ids),batch_size):
        i_end = min(i+batch_size, len(string_ids))
        ids = string_ids[i:i_end]
        # print(len(ids))
        results = get_filtered_string_interactors_STRINGAPI(ids,score_thresh=score_thresh)
        batched_results += results
        if i % 1000 == 0:
            print(i,'/',len(string_ids), end='\r')

    ### If proteins are found 
    if len(batched_results) > 0:
        result_df = pd.DataFrame(batched_results)
        result_df.columns = ['query_ensp', 'query_name', 
                            'partner_ensp', 'partner_name', 'combined_score']
        return result_df
    else:
        print("No results!!")
        return None
    
def k_hop_interactors_STRINGAPI(string_ids, k, score_thresh, debug=False, return_counts=False):
    '''
    Get interactors for STRING IDs using the API within k-hops, 
    only including interacting partners above score_threshold 
    '''
    all_proteins = set(string_ids)
    next_query_proteins = all_proteins
    return_df = None
    counts = []
    for i in range(k):
        # get interactors for the new proteins
        string_interactors_df = batch_filtered_queries_STRINGAPI(next_query_proteins, score_thresh=score_thresh)
        if string_interactors_df is not None:
            # append df
            if return_df is not None:
                return_df.append(string_interactors_df, ignore_index=True)
            else:
                return_df = string_interactors_df

            all_resulting_proteins = set(string_interactors_df['query_ensp']).union(set(string_interactors_df['partner_ensp']))

            # see how many we added
            new_all_proteins = all_proteins.union(all_resulting_proteins)
            next_query_proteins = new_all_proteins.difference(all_proteins)
            all_proteins = new_all_proteins
            if debug:
                print("%d total\t%d new proteins added for %d hop"%(len(all_proteins),len(next_query_proteins), i+1))
            counts += [(i+1,len(all_proteins),len(next_query_proteins))]
    if debug:
        print("%d total proteins"%(len(all_proteins)))
    if return_counts:
        return return_df, counts    
    return return_df

# k_hop_interactors_STRINGAPI(string_ids)  

