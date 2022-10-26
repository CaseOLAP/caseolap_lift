import xml.etree.ElementTree as ET
import json
import requests as req

import sys
# setting path
sys.path.append('..')

from utils.biomedkg_utils import *



'''
TODO: 
    - remove hard-coded file paths and replace w/ variable
    - strange bug with reloading json file, make sure it's parseable
    
Note for Irsyad: The files you want are (1) the same name, (2) edges_meshtree_to_meshtree.csv
1) /DYLAN-Projects/Drug Repurposing Knowledge Graph/output/disease2disease/edges_meshtree-IS-meshid.csv 
2) "/DYLAN-Projects/Drug Repurposing Knowledge Graph/output/disease2disease/edges_meshtree2meshtree_hierarchy.csv"
'''

def process_mesh_mapping(mesh_tree_file = "../data/desc2022.xml"):
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
    json.dump(tree2name, open('../data/meshtree2meshname.json','w'))

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
    output_edgefile_onerel_noweight(outpath = '../data/edges_meshtree-IS-meshid_disease.csv',
                                    columns = ['Disease (MeSH Tree)','Disease (MeSH)','Relationship'],
                                    dictionary = tree2id,
                                    rel = '-is-',
                                    prefix_col1 = 'MeSH_Tree_Disease:',
                                    prefix_col2 = 'MeSH_Disease:',
                                    edges_folder=False)

    df = pd.read_csv('../data/edges_meshtree-IS-meshid_disease.csv')
    #
    # MeSH Term -[is]- MeSH ID
    with open('../data/meshterm-IS-meshid.json','w') as fout:
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
    output_edgefile_onerel_noweight(outpath = '../data/edges_meshtree_to_meshtree.csv',
                                    columns = ['Disease (MeSH Tree)','Disease (MeSH Tree)','Relationship'],
                                    dictionary = tree2tree,
                                    rel = '-subclass_of->',
                                    prefix_col1 = 'MeSH_Tree_Disease:',
                                    prefix_col2 = 'MeSH_Tree_Disease:',
                                    edges_folder=False)

    df = pd.read_csv('../data/edges_meshtree_to_meshtree.csv')

#TODO remove all hard-coded file-paths
