import xml.etree.ElementTree as ET
import shutil
import gzip
import os
import sys

# setting path
sys.path.append('..')

from utils.biomedkg_utils import *
from utils.biomed_apis import *
from utils.other_functions import *


def process_mesh_mapping(mesh_tree_file="../data/desc2023.xml", output_folder="../parsed_mappings/MeSH/"):
    '''
    This function parses the mesh tree file and parses three major files:
    1) ../data/edges_meshtree_to_meshtree.csv the hierarchical structure of mesh tree
    2) ../data/meshtree2meshname.json the meshtree code to mesh name mappings
    3) ../data/edges_meshtree-IS-meshid_disease.csv the mesh tree to mesh id edges
    :param mesh_tree_file:
    :return:
    '''
    tree = ET.parse(mesh_tree_file)
    root = tree.getroot()

    name2id, id2name, id2tree, tree2id = dict(), dict(), dict(), dict()
    all_tree_numbers = list()

    for ele in root:
        try:
            # MeSH Tree Number
            tree_numbers = ele.find('TreeNumberList').findall('TreeNumber')

            # If anatomy
            for tree_number in tree_numbers:
                if tree_number.text.startswith(('C', 'F03')):

                    '''
                    Tree
                    '''
                    tree_number = tree_number.text
                    all_tree_numbers.append(tree_number)

                    '''
                    ID to Tree
                    '''
                    try:
                        # MeSH ID
                        ID = ele.find('DescriptorUI').text

                        # MeSH ID -[is]- MeSH Tree
                        id2tree.setdefault(ID, set()).add(tree_number)
                        tree2id.setdefault(tree_number, set()).add(ID)
                    except:
                        procced = True

                    '''
                    ID to Name
                    '''
                    try:
                        # MeSH ID
                        ID = ele.find('DescriptorUI').text

                        # MeSH Term Name
                        name = ele.find('DescriptorName').find('String').text

                        # MeSH Term -[is]- MeSH ID
                        name2id.setdefault(name, set()).add(ID)
                        id2name.setdefault(ID, set()).add(name)
                    except:
                        proceed = True
        except:
            continue

    all_tree_numbers = sorted(all_tree_numbers)
    tree2id = dict(sorted(tree2id.items()))

    for k, v in name2id.copy().items():
        name2id[k] = list(name2id[k])

    '''MeSH Tree is MeSH Term/Name'''
    tree2id = switch_dictset_to_dictlist(tree2id)
    tree2name = dict()

    for tree, ID in tree2id.items():
        assert len(ID) == 1

        name_list = list(id2name[ID[0]])
        assert len(name_list) == 1
        name = name_list[0]

        tree2name.setdefault(tree, list()).append(name)

    for tree, name in tree2name.items():
        assert len(name) == 1

    for tree, name in tree2name.copy().items():
        tree2name[tree] = name[0]

    json.dump(tree2name, open(os.path.join(output_folder,'meshtree2meshname.json'), 'w'))

    tree_nums = list({k: v for k, v in tree2name.items() if k.startswith('C14')})
    total_trees = len(tree_nums)

    id2synonyms = dict()

    for index, tree_num in enumerate(tree_nums):

        mesh_id = tree2id[tree_num][0]
        r = req.get('https://id.nlm.nih.gov/mesh/lookup/details?descriptor=' + mesh_id).json()

        for entry in r['terms']:
            id2synonyms.setdefault(mesh_id, set()).add(entry['label'])

        print("MeSH IDs parsed: ",index, '/', total_trees, end='\r')

    # MeSH Tree Number -[is]- MeSH ID
    output_edgefile_onerel_noweight(outpath=os.path.join(output_folder,'edges_meshtree-IS-meshid_disease.csv'),
                                    columns=['Disease (MeSH Tree)', 'Disease (MeSH)', 'Relationship'],
                                    dictionary=tree2id,
                                    rel='-is-',
                                    prefix_col1='MeSH_Tree_Disease:',
                                    prefix_col2='MeSH_Disease:',
                                    edges_folder=False)

    # df = pd.read_csv('../data/edges_meshtree-IS-meshid_disease.csv')
    #
    # MeSH Term -[is]- MeSH ID
    with open(os.path.join(output_folder,'meshterm-IS-meshid.json'), 'w') as fout:
        json.dump(name2id, fout)

        # mesh_name2id = json.load(open('../data/meshterm-IS-meshid.json')) #TODO bug for some reason

    tree2tree = dict()

    # Tree Number
    for tree_num in all_tree_numbers:
        if '.' in tree_num:

            # Parent of Tree Number
            parent = ''
            for num in tree_num.split('.')[:len(tree_num.split('.')) - 1]:
                parent += num + '.'
            parent = parent.strip('.')

            # Tree Number -[subclass of]-> Tree Number
            tree2tree[tree_num] = [parent]

    # MeSH Tree Number -[subclass of]-> MeSH Tree Number
    output_edgefile_onerel_noweight(outpath=os.path.join(output_folder,'edges_meshtree_to_meshtree.csv'),
                                    columns=['Disease (MeSH Tree)', 'Disease (MeSH Tree)', 'Relationship'],
                                    dictionary=tree2tree,
                                    rel='-subclass_of->',
                                    prefix_col1='MeSH_Tree_Disease:',
                                    prefix_col2='MeSH_Tree_Disease:',
                                    edges_folder=False)

    # df = pd.read_csv('../data/edges_meshtree_to_meshtree.csv')

# TODO remove all hard-coded file-paths


def download_file(url, directory):

    local_filename = url.split('/')[-1]
    PATH = os.path.join(directory,local_filename).replace('?','_') # for grndb files
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(PATH, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return local_filename


def gunzip(file_path, output_path):
    with gzip.open(file_path, "rb") as f_in, open(output_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


def filter_REACTOME_HUMAN(read_file, write_file):
    with open(read_file, "r") as f1:
        with open(write_file, "w") as f2:
            for line in f1:
                line_list = line.split()
                if "R-HSA" in line_list[0] and "R-HSA" in line_list[1]:
                    f2.write(line)


def check_required_files(required_files, directory, debug=False):
    '''
    This function checks to see if all prerequisite input files have been downloaded
    :param data_folder:
    :return:
    '''

    # make data folders if doesn't exist
    if not os.path.exists(directory):
        os.makedirs(directory)
    for resource in required_files.keys():
        resource_folder = os.path.join(directory, resource)
        if not os.path.exists(resource_folder):
            os.makedirs(resource_folder)

    resource_to_data_file_bool = {}
    resources_ready_list = []
    resources_not_ready_list = []
    ready = True # ready if all resources pass check
    for resource, data_files in required_files.items():

        # search if all data files exist for a resource
        exists_list = []
        for data_file in data_files:
            resource_folder = os.path.join(directory, resource)
            data_file_path = os.path.join(resource_folder, data_file)
            exists = os.path.exists(data_file_path)
            exists_list += [exists]

        # resource is ready if all data files are downloaded
        resource_ready = all(exists_list)
        if resource_ready:
            resources_ready_list += [resource]
        else:
            resources_not_ready_list += [resource]
        ready = ready and resource_ready

        if debug:
            print("*** %s ready: %s ***" %(resource, str(resource_ready)))
            for data_file, exists in zip(data_files,exists_list):
                print(data_file,exists)

        resource_to_data_file_bool[resource] = {data_file:exists for data_file, exists in zip(data_files,exists_list)}

    return resource_to_data_file_bool, ready


def map_protein2go_ids(goa_file = './data/GO/goa_human.gaf',output_folder="./parsed_mappings/GO"):
    relations, go_terms, total = dict(), set(), 0
    protein2go = dict()
    go2protein = dict()

    for i, line in enumerate(open(goa_file)):
        if i > 40:
            line = line.split('\t')
            protein = line[1]
            relation = line[3]
            go_term = line[4]

            relations[relation] = relations.get(relation, 0) + 1
            go_terms.add(go_term)
            total += 1

            protein2go.setdefault(protein, list()).append(go_term)
            go2protein.setdefault(go_term, list()).append(protein)

    json.dump(protein2go, open(os.path.join(output_folder, 'protein2go.json'), 'w'))
    json.dump(go2protein, open(os.path.join(output_folder, 'go2protein.json'), 'w'))

    return protein2go, go2protein


def load_protein2pathway_data(uniprot2reactome_file_path='./data/Reactome/UniProt2Reactome.txt', output_folder="./parsed_mappings/Reactome"):
    protein2pathway, pathway2protein = dict(), dict()

    for line in open(uniprot2reactome_file_path):
        line = line.strip().split('\t')

        # Protein -involved in-> Pathway
        try:
            species, protein, pathway = line[5].strip(), line[0], line[1]

            if species.lower() == 'homo sapiens':
                protein2pathway.setdefault(protein, list()).append(pathway)
                pathway2protein.setdefault(pathway, list()).append(protein)
        except:
            pass

    # Export
    protein2pathway = switch_dictlist_to_dictset(protein2pathway)
    protein2pathway = switch_dictset_to_dictlist(protein2pathway)

    pathway2protein = switch_dictlist_to_dictset(pathway2protein)
    pathway2protein = switch_dictset_to_dictlist(pathway2protein)

    json.dump(protein2pathway, open(os.path.join(output_folder, 'protein2pathway.json'), 'w'))
    json.dump(pathway2protein, open(os.path.join(output_folder, 'pathway2protein.json'), 'w'))

    return protein2pathway, pathway2protein


def download_data(resource_to_data_file_bool, data_folder, file_to_link_file):
    '''
    This function downloads a data file if it is not available.
    resource_to_data_file_bool is a dict of resource (e.g. MeSH, GO, Reactome)
    to a dict between its data_file and a boolean if it exists or not.
    :param resource_to_data_file_bool:
    :return:
    '''
    file_to_link = json.load(open(file_to_link_file,'r'))

    for resource, data_file_bool_dict in resource_to_data_file_bool.items():
        for data_file, exists in data_file_bool_dict.items():
            if not exists:
                resource_folder = os.path.join(data_folder,resource)
                data_file_path = os.path.join(resource_folder,data_file)
    
                # create data and resource folder if not exists
                if not os.path.exists(data_folder):
                    os.makedirs(data_folder)
                if not os.path.exists(resource_folder):
                    os.makedirs(resource_folder)

                print("Downloading %s"%(data_file_path))
                # GRNdb handled separately, a list of files
                if data_file == 'GRNdb':
                    if not os.path.exists(data_file_path):
                        os.makedirs(data_file_path)
                    for download_link in file_to_link[data_file]:
                        local_file_name = download_file(download_link,data_file_path)
                # download file
                elif data_file in file_to_link:
                    download_link = file_to_link[data_file]
                    local_file_name = download_file(download_link, resource_folder)

                    # unzip files
                    if local_file_name.endswith('.gz'):
                        zipped_file_path = os.path.join(resource_folder,local_file_name)
                        gunzip(zipped_file_path, data_file_path)
                        os.remove(zipped_file_path)

                # update flag
                resource_to_data_file_bool[resource][data_file] = os.path.exists(data_file_path)
    return resource_to_data_file_bool


def get_short_full_names(names, section_key, section_values):
    '''
    FUNCTION:
    - Get all names from an entry

    PARAMS:
    - names (list)
    '''
    # Full name
    if section_key == 'fullName':
        name = section_values['value']
        names.append(name)

    # Short names
    if section_key == 'shortNames':
        for entry in section_values:
            name = entry['value']
            names.append(name)

            # EC Numbers
    if section_key == 'ecNumbers':
        for entry in section_values:
            name = entry['value']
            names.append(name)

    return names


def get_protein_names(entry, gene_names):
    '''
    FUNCTION:
    - Taking entries in the UniProt API's results, parse the entries
      and extract the protein names.

    PARAMS:
    - entry: The UniProt API's results
    - gene_names (bool): Indicates whether you want the gene_names included
    '''

    protein_id = entry['from']

    ''' Protein Names '''
    names = list()

    ''' Recommended name '''
    try:
        reccomended_name_section = entry['to']['proteinDescription']['recommendedName']
        for section_key, section_values in reccomended_name_section.items():
            names = get_short_full_names(names, section_key, section_values)
    except:
        pass

    ''' Alternative Names '''
    try:
        alternative_name_section = entry['to']['proteinDescription']['alternativeNames']
        for sub_entry in alternative_name_section:
            for section_key, section_values in sub_entry.items():
                names = get_short_full_names(names, section_key, section_values)
    except:
        pass

    ''' Submitted Names '''
    try:
        sub_name_section = entry['to']['proteinDescription']['submissionNames']
        for sub_entry in sub_name_section:
            for section_key, section_values in sub_entry.items():
                names = get_short_full_names(names, section_key, section_values)
    except:
        pass

    '''Gene Names'''
    if gene_names == True:
        # Data for gene names
        try:
            rgenes = entry['to'][0]['genes'][0]
        except:
            pass

        # Get gene names
        for nametype, main_entry in rgenes.items():
            if type(main_entry) == dict:
                name = main_entry['value']
            elif type(main_entry) == list:
                for entry in main_entry:
                    name = entry['value']
            names.append(name)

    return protein_id, names


def curl_uniprot_api(ids, from_id, to_id, jobfile):
    '''Submit job (API)'''
    try:
        os.system('rm %s' % jobfile)
    except:
        pass


    BATCH_SIZE = 1000
    if len(ids) < BATCH_SIZE:
        BATCH_SIZE = len(ids) - 1

    tot_reqs = int(len(ids) / BATCH_SIZE) + 1
    for i in range(0, tot_reqs):

        # Batch
        if i != tot_reqs - 1:
            ids_job = ','.join(list(ids)[i * BATCH_SIZE:(i + 1) * BATCH_SIZE])
        else:
            ids_job = ','.join(list(ids)[i * BATCH_SIZE:])

        # API
        command = 'curl --form "from=%s" --form "to=%s" --form "ids=%s" "https://rest.uniprot.org/idmapping/run" >> %s' % (
            from_id,to_id,ids_job, jobfile)
        #print(command)
        os.system(command)
        if i != tot_reqs - 1:
            # temp for splitting the lines, quick windows/linux patch
            os.system("echo @ >> %s" % jobfile)
            # os.system("echo. >> %s" % jobfile) # windows only, here
    # remove the temp character and replace file
    os.system("cat %s | sed 's/@//' > temp.json" %jobfile)
    os.system("mv temp.json %s"%jobfile)

    ''' Get results of job '''
    results = dict()
    count=0
    with open(jobfile) as fin:
        for line in fin:
            job_id = json.loads(line)['jobId']
            results[job_id] = dict()
            print("Job %d of %d. Job ID:%s"%(count, tot_reqs, job_id))
            if check_id_mapping_results_ready_uniprot_api(job_id):
                link = get_id_mapping_results_link_uniprot_api(job_id)
                results[job_id] = get_id_mapping_results_search_uniprot_api(link)
            count+=1
    return results


def get_protein_id_to_synonyms(protein_ids,mapping_folder='../parsed_mappings/'):
    # Protein IDs to submit
    prot_list = protein_ids

    # Submit IDs to API to get names
    uniprot_results = curl_uniprot_api(prot_list, from_id='UniProtKB_AC-ID', to_id='UniProtKB',
                     jobfile=os.path.join(mapping_folder,'protein_id2entry.json'))

    # Extract IDs->Names
    protein_id2names = dict()
    for job_id, entries in uniprot_results.items():
        for entry in entries['results']:
            protein_id, names = get_protein_names(entry, gene_names=False)
            for name in names:
                protein_id2names.setdefault(protein_id, set()).add(name)
    protein_id2names = switch_dictset_to_dictlist(protein_id2names)
    return protein_id2names


def get_uniprot_from_gene_names(gene_names,mapping_folder="../parsed_mappings/"):
    results = curl_uniprot_api(gene_names, "Gene_Name", "UniProtKB-Swiss-Prot",
                               os.path.join(mapping_folder,'curl_uniprot_job_ids_geneid2uniprotkb.json'))

    '''Dictionary from API data'''
    gene_name_2_protein_id = dict()
    for key in results.keys():
        gene_name_2_protein_id[key] = get_to_uniprotid_from_genename_mapping_dict_uniprot_api(
            results[key],
            [9606],
            filter_by_reviewed=True)

    removed_genes = 0
    for key, batch in gene_name_2_protein_id.items():
        for gene_name, protein_ids in batch.copy().items():
            if len(protein_ids) > 1:
                gene_name_2_protein_id[key].pop(gene_name)
                removed_genes += 1
        print(removed_genes, 'genes removed because they were ambiguous\n' + \
              'Names mapped to multiple proteins (possibly ok, but possibly not)')

    new_gene_name_2_protein_id = dict()

    for each_list in list(gene_name_2_protein_id.values()):
        for gene_name, protein_ids in each_list.items():
            for protein_id in protein_ids:
                new_gene_name_2_protein_id.setdefault(gene_name, set()).add(protein_id)
    gene_name_2_protein_id = new_gene_name_2_protein_id.copy()

    gene_name_2_protein_id = switch_dictset_to_dictlist(gene_name_2_protein_id)

    json.dump(gene_name_2_protein_id, open(os.path.join(mapping_folder,'gene_name_2_protein_id.json'), 'w'))

    return gene_name_2_protein_id


def extract_uniprot_to_entrez_mapping(protein_list, debug=False,
                                      output_folder='./parsed_mappings/Transcription_Factor_Dependence'):
    # Map proteins to GeneID
    job_id = submit_id_mapping_uniprot_api(
        from_db='UniProtKB_AC-ID',
        to_db='GeneID',
        ids=protein_list)

    # This checks on the job until it is finished
    if check_id_mapping_results_ready_uniprot_api(job_id):
        link = get_id_mapping_results_link_uniprot_api(job_id)
        results = get_id_mapping_results_search_uniprot_api(link)

    gene_id_2_protein_id = dict()
    protein_id_2_gene_id = dict()

    for protein_to_gene in results['results']:
        protein_id = protein_to_gene['from'].strip()
        gene_id = protein_to_gene['to'].strip()
        if gene_id != '' and protein_id != '':
            gene_id_2_protein_id.setdefault(gene_id, set()).add(protein_id)
            protein_id_2_gene_id.setdefault(protein_id, set()).add(gene_id)

    if debug:
        print('\nAfter UniProt API')
        print('Entrez-is-UniProt', len(gene_id_2_protein_id),
              'UniProt-is-Entrez', len(protein_id_2_gene_id))

    # output
    gene_id_2_protein_id = switch_dictset_to_dictlist(gene_id_2_protein_id)
    protein_id_2_gene_id = switch_dictset_to_dictlist(protein_id_2_gene_id)

    json.dump(gene_id_2_protein_id, open(os.path.join(output_folder, 'all_entrez2uniprot.json'), 'w'))
    json.dump(protein_id_2_gene_id, open(os.path.join(output_folder, 'all_uniprot2entrez.json'), 'w'))


def extract_proteins_from_fasta(fasta_file):
    headers = [l.strip("\n") for l in open(fasta_file,"r").readlines() if ">" in l]
    proteins = set([h.split("|")[1] for h in headers])
    return proteins


def extract_uniprot_to_entrez_mapping(protein_list, debug=False,
                                      output_folder='./parsed_mappings/Transcription_Factor_Dependence'):
    # Map proteins to GeneID
    job_id = submit_id_mapping_uniprot_api(
        from_db='UniProtKB_AC-ID',
        to_db='GeneID',
        ids=protein_list)

    # This checks on the job until it is finished
    if check_id_mapping_results_ready_uniprot_api(job_id):
        link = get_id_mapping_results_link_uniprot_api(job_id)
        results = get_id_mapping_results_search_uniprot_api(link)

    gene_id_2_protein_id = dict()
    protein_id_2_gene_id = dict()

    for protein_to_gene in results['results']:
        protein_id = protein_to_gene['from'].strip()
        gene_id = protein_to_gene['to'].strip()
        if gene_id != '' and protein_id != '':
            gene_id_2_protein_id.setdefault(gene_id, set()).add(protein_id)
            protein_id_2_gene_id.setdefault(protein_id, set()).add(gene_id)

    if debug:
        print('\nAfter UniProt API')
        print('Entrez-is-UniProt', len(gene_id_2_protein_id),
              'UniProt-is-Entrez', len(protein_id_2_gene_id))

    # output
    gene_id_2_protein_id = switch_dictset_to_dictlist(gene_id_2_protein_id)
    protein_id_2_gene_id = switch_dictset_to_dictlist(protein_id_2_gene_id)

    return gene_id_2_protein_id, protein_id_2_gene_id


def prepare_tfd_mappings(input_folder = '../data/Transcription_Factor_Dependence',
                         output_folder='../parsed_mappings/Transcription_Factor_Dependence',
                         debug=False):

    # required files: GRNdb folder, UniProt Human Proteome fasta file
    grndb_folder = os.path.join(input_folder,'GRNdb')
    fasta_file = os.path.join(input_folder,'UP000005640_9606.fasta')
    # check they exist or exit
    
    # read GRNdb folder
    tf_gene_name_2_target_gene_name = dict()  # TF gene name to target gene name
    gene_names = set()  # Gene Names (Tfs and targets)
    tf_gene_names = set()  # Transcription Factor gene names
    target_gene_names = set()  # Target gene names
    for file in os.listdir(grndb_folder):
        if 'txt' not in file:
            continue

        for line in open(os.path.join(grndb_folder, file)):
            line = line.strip().split('\t')
            if "There is something wrong" in line[0]:
                break

            confidence = line[5]
            if confidence == 'High':
                # Gene names
                tf_gene_name = line[0]
                targ_gene_name = line[1]

                # Save gene names
                gene_names.add(tf_gene_name)
                gene_names.add(targ_gene_name)
                tf_gene_names.add(tf_gene_name)
                target_gene_names.add(targ_gene_name)

                # TF Gene -targets-> Target Gene
                tf_gene_name_2_target_gene_name.setdefault(tf_gene_name, set()).add(targ_gene_name)

    # Change the values from set into a list
    tf_gene_name_2_target_gene_name = switch_dictset_to_dictlist(tf_gene_name_2_target_gene_name)

    # map genes to protein ids
    gene_name_2_protein_id = get_uniprot_from_gene_names(gene_names,mapping_folder=output_folder)

    # extract proteins from uniprot fasta
    all_proteins = extract_proteins_from_fasta(fasta_file)

    # get entrez to uniprot mappings
    gene_id_2_protein_id,protein_id_2_gene_id = extract_uniprot_to_entrez_mapping(all_proteins, debug=True)

    # get uniprot to name mappings
    protein_id2names = get_protein_id_to_synonyms(all_proteins,mapping_folder=output_folder)

    # output
    json.dump(tf_gene_name_2_target_gene_name, open(os.path.join(output_folder,'tf_gene_name_2_target_gene_name.json'),'w'))
    json.dump(gene_name_2_protein_id, open(os.path.join(output_folder, 'gene_name_2_protein_id.json'), 'w'))
    json.dump(protein_id2names, open(os.path.join(output_folder, 'id2synonyms_not_case_varied.json'), 'w'))
    json.dump(gene_id_2_protein_id, open(os.path.join(output_folder, 'all_entrez2uniprot.json'), 'w'))
    json.dump(protein_id_2_gene_id, open(os.path.join(output_folder, 'all_uniprot2entrez.json'), 'w'))


#TODO move to prepare_kg_data.py or utils
def extract_go_hierarchy(go_dict):
    go_id_to_links = {}

    for go_term, go_fields in go_dict.items():

        parents = []
        if 'is_a' in go_fields:
            parents += [p.split(" ! ")[0] for p in go_fields['is_a']]
        if 'relationship' in go_fields:
            parents += [p.split(" ! ")[0].strip('part_of ') for p in go_fields['relationship'] if 'part_of' in p]

        # add to dict
        for g in [go_term] + parents:
            if g not in go_id_to_links:
                go_id_to_links[g] = dict()
                go_id_to_links[g]['parents'] = []
                go_id_to_links[g]['children'] = []
        for p in parents:
            go_id_to_links[p]['children'] += [go_term]
            go_id_to_links[go_term]['parents'] += [p]
    return go_id_to_links


def parse_downloaded_data(resource_to_processed_file_bool, mapping_folder, data_folder, debug=False):

    # find out which resources are missing mappings
    missing_mappings = {}
    for resource, mapping_file_bool_dict in resource_to_processed_file_bool.items():
        missing_files = not(all([b for b in mapping_file_bool_dict.values()]))
        missing_mappings[resource] = missing_files
        if debug and missing_files:
            print("%s missing files!"%resource)

    for resource, missing_file in missing_mappings.items():

        if missing_file:
            if debug:
                print("Missing file for %s"%resource)
            # parse mappings if missing
            if resource == 'GO':
                if debug:
                    print("Parsing GO mappings")
                input_folder = os.path.join(data_folder,'GO')
                output_folder = os.path.join(mapping_folder,'GO')
                prepare_go_mappings(input_folder=input_folder, output_folder=output_folder)

            if resource == 'MeSH':
                if debug:
                    print("Parsing MeSH mappings")
                input_folder = os.path.join(data_folder,'MeSH')
                output_folder = os.path.join(mapping_folder,'MeSH')
                mesh_tree_file = os.path.join(input_folder,'desc2023.xml')
                process_mesh_mapping(mesh_tree_file=mesh_tree_file, output_folder=output_folder)

            if resource == 'Reactome':
                if debug:
                    print("Parsing Reactome mappings")
                input_folder = os.path.join(data_folder,'Reactome')
                output_folder = os.path.join(mapping_folder,'Reactome')
                reactome_to_protein_file = os.path.join(input_folder, 'UniProt2Reactome.txt')
                load_protein2pathway_data(reactome_to_protein_file, output_folder=output_folder)

            if resource == 'Transcription_Factor_Dependence':
                if debug:
                    print("Parsing Transcription_Factor_Dependence mappings")
                input_folder = os.path.join(data_folder, 'Transcription_Factor_Dependence')
                output_folder = os.path.join(mapping_folder, 'Transcription_Factor_Dependence')
                prepare_tfd_mappings(input_folder=input_folder,output_folder=output_folder)

def prepare_go_mappings(input_folder = '../data/GO/', output_folder = '../parsed_mappings/GO/'):
    goa_file = os.path.join(input_folder, 'goa_human.gaf')
    go_basic_file = os.path.join(input_folder,'go-basic.obo')

    protein2go, go2protein = map_protein2go_ids(goa_file=goa_file, output_folder=output_folder)

    # GO Term information
    go_dict = get_go_term_information(gobasic_file=go_basic_file)

    # GO Term ID -> GO Term Name
    go_id2term = map_go_id2term(go_dict)

    # Separate GO terms into 3 lists for each type
    go_bio_proc, go_cell_comp, go_mol_func = separate_go_terms_into_three_lists_for_each_tree(go_dict)

    # Cellular Compartment Names
    cell_comp_name = list()
    cell_comp_id2name = dict()

    for ID in go_cell_comp:
        name = go_id2term[ID]
        cell_comp_name.append(name)
        cell_comp_id2name[ID] = name

    # GO Term Cellular Compartments -> Protein
    go_term_cell_comp2protein, protein2go_term_cell_comp = get_protein2go_term(protein2go, cell_comp_id2name)

    # extract go hierarchy (is_a and part_of)
    go_id_to_links = extract_go_hierarchy(go_dict)

    json.dump(go_id_to_links, open(os.path.join(output_folder,"go_id_to_links.json"),'w'))
    print(len(go_id_to_links))
    # json.dump(go_term_cell_comp2protein, open(os.path.join(output_folder,"go_term_cell_comp2protein.json"),'w')) # Don't think we need this


def prepare_knowledge_base_data(data_folder, mapping_folder, file_to_link_file,
                                include_reactome=True, include_tfd = True,
                                redownload=False, remap=False, debug=False):
    '''
    This function checks the data_folder and mapping_folder if the required files are downloaded and parsed, respectively.
    redownload flag means it will delete/overwrite existing files
    :param data_folder:
    :param mapping_folder:
    :param redownload:
    :return:
    '''
    # make folders if doesn't exist
    if not os.path.exists(data_folder):
        os.makedirs(data_folder)
    if not os.path.exists(mapping_folder):
        os.makedirs(mapping_folder)

    required_files = {'MeSH': ['desc2023.xml', 'mtrees2023.bin'], #TODO MeSH dates should be generalized
                      'GO': ['go-basic.obo', 'goa_human.gaf'],
                      'Reactome': ['UniProt2Reactome.txt', 'ReactomePathwaysRelation.txt','UniProt2Reactome_All_Levels.txt', 'ReactomePathways.txt'],
                      'Transcription_Factor_Dependence': ['GRNdb','UP000005640_9606.fasta']
                      }
    processed_files = {'MeSH': ['meshtree2meshname.json', 'edges_meshtree-IS-meshid_disease.csv',
                                'meshterm-IS-meshid.json','edges_meshtree_to_meshtree.csv'], #TODO meshterms_per_cat.json?
                       'GO': ['go2protein.json','protein2go.json','go_id_to_links.json'],
                       'Reactome': ['pathway2protein.json','protein2pathway.json'],
                       'Transcription_Factor_Dependence': ['all_entrez2uniprot.json','all_uniprot2entrez.json','id2synonyms_not_case_varied.json','gene_name_2_protein_id.json','tf_gene_name_2_target_gene_name.json']
                       }
    # remove reactome and tfd if not included in analysis
    if not include_reactome:
        required_files.pop('Reactome')
        processed_files.pop('Reactome')
    if not include_tfd:
        required_files.pop("Transcription_Factor_Dependence")
        processed_files.pop("Transcription_Factor_Dependence")

    ### Downloading data ###
    # determine which files need to be downloaded
    if redownload:
        # set everything to false, redownload everything
        download_complete = False
        resource_to_data_file_bool = {}
        for resource, data_file_list in required_files.items():
            resource_to_data_file_bool[resource] = {data_file: False for data_file in data_file_list}
    else:
        # check to see which files are downloaded
        resource_to_data_file_bool,download_complete = check_required_files(required_files, data_folder, debug=debug)

    if not download_complete:
        # download files
        download_data(resource_to_data_file_bool, data_folder, file_to_link_file)

        # check to see which downloaded files still are not available
        resource_to_data_file_bool,download_complete_after_download = check_required_files(required_files, data_folder, debug=debug)

        if not download_complete_after_download:
            print("Error! Missing required biomedical knowledgebase data")
            _, _ = check_required_files(required_files, data_folder, debug=True)
            sys.exit(1)
    if debug:
        print("Data downloading phase completed.")

    ### Parsing data ###
    resource_to_processed_file_bool, processed_files_complete = check_required_files(processed_files, mapping_folder, debug=debug)
    # force remapping if flagged
    if remap:
        temp_resource_map = {}
        for resource, file_map_dict in resource_to_processed_file_bool.items():
            temp_resource_map[resource] = {k:False for k in file_map_dict.keys()}
        resource_to_processed_file_bool = temp_resource_map
        # resource_to_processed_file_bool = {resource{mfile:False} for  in resource_to_processed_file_bool.keys()}
        processed_files_complete = False
    if redownload or not processed_files_complete:
        # process files
        parse_downloaded_data(resource_to_processed_file_bool, mapping_folder, data_folder,debug=debug)

        # check to see which processed files still are not available
        resource_to_processed_file_bool, processed_files_complete = check_required_files(processed_files, mapping_folder,
                                                                                            debug=debug)
        if not processed_files_complete:
            print("Error! Could not process input data files")
            _, _ = check_required_files(processed_files, mapping_folder, debug=True)
            sys.exit(1)
    if debug:
        print("Parsing data phase completed.")

    return True
