import xml.etree.ElementTree as ET
import shutil
import gzip
import os
import sys

# setting path
sys.path.append('..')

from utils.biomedkg_utils import *


def process_mesh_mapping(mesh_tree_file="../data/desc2022.xml", output_folder="../parsed_mappings/MeSH/"):
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

    print(tree2name)
    json.dump(tree2name, open(os.path.join(output_folder,'meshtree2meshname.json'), 'w'))

    tree_nums = list({k: v for k, v in tree2name.items() if k.startswith('C14')})
    total_trees = len(tree_nums)

    id2synonyms = dict()

    for index, tree_num in enumerate(tree_nums):

        mesh_id = tree2id[tree_num][0]
        r = req.get('https://id.nlm.nih.gov/mesh/lookup/details?descriptor=' + mesh_id).json()

        for entry in r['terms']:
            id2synonyms.setdefault(mesh_id, set()).add(entry['label'])

        print(index, '/', total_trees, end='\r')

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
    PATH = os.path.join(directory,local_filename)
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


def download_data(resource_to_data_file_bool, data_folder):
    '''
    This function downloads a data file if it is not available.
    resource_to_data_file_bool is a dict of resource (e.g. MeSH, GO, Reactome)
    to a dict between its data_file and a boolean if it exists or not.
    :param resource_to_data_file_bool:
    :return:
    '''
    file_to_link = {'go-basic.obo': "http://purl.obolibrary.org/obo/go/go-basic.obo",
                    'goa_human.gaf': "http://geneontology.org/gene-associations/goa_human.gaf.gz",
                    'UniProt2Reactome.txt': "https://reactome.org/download/current/UniProt2Reactome.txt",
                    'ReactomePathwaysRelation.txt': 'https://reactome.org/download/current/ReactomePathwaysRelation.txt',
                    'mtrees2021.bin':"https://nlmpubs.nlm.nih.gov/projects/mesh/2021/meshtrees/mtrees2021.bin",
                    }

    #TODO all_entrez2uniprot.json not downloaded TODO TODO

    for resource, data_file_bool_dict in resource_to_data_file_bool.items():
        for data_file, exists in data_file_bool_dict.items():
            if not exists:
                resource_folder = os.path.join(data_folder,resource)
                data_file_path = os.path.join(resource_folder,data_file)

                # download file
                print("Downloading %s"%(data_file_path))
                if data_file in file_to_link:
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
            print("Missing file for %s"%resource)
            # parse mappings if missing
            if resource == 'GO':
                if debug:
                    print("Parsing GO mappings")
                input_folder = os.path.join(data_folder,'GO')
                output_folder = os.path.join(mapping_folder,'GO')
                goa_file = os.path.join(input_folder,'goa_human.gaf')
                map_protein2go_ids(goa_file=goa_file, output_folder=output_folder)

            if resource == 'MeSH':
                if debug:
                    print("Parsing MeSH mappings")
                input_folder = os.path.join(data_folder,'MeSH')
                output_folder = os.path.join(mapping_folder,'MeSH')
                mesh_tree_file = os.path.join(input_folder,'desc2022.xml')
                process_mesh_mapping(mesh_tree_file=mesh_tree_file, output_folder=output_folder)

            if resource == 'Reactome':
                if debug:
                    print("Parsing Reactome mappings")
                input_folder = os.path.join(data_folder,'Reactome')
                output_folder = os.path.join(mapping_folder,'Reactome')
                reactome_to_protein_file = os.path.join(input_folder, 'UniProt2Reactome.txt')
                load_protein2pathway_data(reactome_to_protein_file, output_folder=output_folder)

            if resource == 'Transcription_Factor_Dependence':
                input_folder = os.path.join(mapping_folder, 'Transcription_Factor_Dependence')
                #TODO


def prepare_knowledge_base_data(data_folder, mapping_folder, redownload=False, debug=False):
    '''
    This function checks the data_folder and mapping_folder if the required files are downloaded and parsed, respectively.
    redownload flag means it will delete/overwrite existing files
    :param data_folder:
    :param mapping_folder:
    :param redownload:
    :return:
    '''

    #TODO have the required_files change based on user input
    required_files = {'MeSH': ['desc2022.xml', 'mtrees2021.bin'], #TODO MeSH dates should be generalized
                      'GO': ['go-basic.obo', 'goa_human.gaf'],
                      'Reactome': ['UniProt2Reactome.txt', 'ReactomePathwaysRelation.txt'],
                      'Transcription_Factor_Dependence': ['GRNdb','all_entrez2uniprot.json']
                      }

    processed_files = {'MeSH': ['meshtree2meshname.json', 'edges_meshtree-IS-meshid_disease.csv',
                                'meshterm-IS-meshid.json','edges_meshtree_to_meshtree.csv'], #TODO meshterms_per_cat.json?
                       'GO': ['go2protein.json','protein2go.json'],
                       'Reactome': ['pathway2protein.json','protein2pathway.json'],
                       'Transcription_Factor_Dependence': ['tf_protein_id_2_target_gene_id.json',
                                                           'tf_protein_id_2_target_protein_id.json']
                       }

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
        download_data(resource_to_data_file_bool, data_folder)

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
    return


data_folder = "../data"
mapping_folder = "../parsed_mappings/"
prepare_knowledge_base_data(data_folder, mapping_folder,redownload=False,debug=True)
