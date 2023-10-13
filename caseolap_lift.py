import argparse
import sys
import os
import datetime
import requests
from pathlib import Path
from preprocessing.prepare_knowledge_base_data import  download_data,prepare_knowledge_base_data
import json
import shutil
from preprocessing.entity_set_expansion import prepare_subcellular_compartment_proteins
from text_mining import *
from analysis.run_analyse_results import analyze_results
from text_mining.run_text_mining import run_text_mining
# from kg_creation.caseolapLIFT_kg import *
from utils.FileConverter import *
from kg_creation.caseolapLIFT_kg import caseolapLIFT_knowledge_graph

def parse_abbreviations(input, num_categories, debug=False):

    abbreviations = []
    for abbr in input.split(" "):
        if abbr not in abbreviations:
            abbreviations += [abbr]
        else:
            abbreviations += [abbr+"_1"]
    if debug:
        print("Abbreviation list : %s"%str(abbreviations))

    # error checking
    if len(set(abbreviations)) != num_categories:
        print("Wrong number of abbreviations. Number of abbreviations must match the number of disease categories")
        print("Abbreviation list : %s" % str(abbreviations))
        print("Number of disease categories: "+str(num_categories))
        raise Exception("Invalid abbreviations")

    return abbreviations


def parse_diseases(diseases, valid_mesh_terms, debug=False):
    # check if the user input are all valid MeSH tree ids
    if ".txt" in diseases:
        with open(diseases) as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("Disease Category:"):
                    diseases = parse_diseases(line.split(":")[1], valid_mesh_terms)
    else:
        res = []
        mesh_lists = diseases.split(" ")
        for mesh_list in mesh_lists:
            diseases = mesh_list.split(",")
            print("Disease category %d: %s"%(len(res),str(diseases)))
            for disease in diseases:
                if disease in valid_mesh_terms:
                    if debug:
                        print("The MeSH term: " + disease)
                else:
                    if debug:
                        print("The MeSH term %s is invalid!"%disease)
                    raise Exception("Invalid MeSH term %s"%disease)

            res.append(diseases)
        return res


def parse_subcellular_component(components, valid_go_terms, debug=False):
    # check if user inputs are valid GO terms
    input_components = components.split(" ")
    res = []
    for component in input_components:

        if component in valid_go_terms:
            res.append(component)
            if debug:
                print("The GO cellular component: "+component)
        else:
            if debug:
                print("The GO term %s is invalid!"%(component))
            raise Exception("Invalid GO term %s"%component)

    return res


def parse_protein_list(proteins):

    # do nothing if it's empty
    if proteins == None:
        return

    if ".txt" in proteins:
        with open(proteins) as f:
            lines = f.readlines()
        for line in lines:
            #change it for entities.txt not parameters.txt
            #get rid of the header (protein list)
            if line.startswith("Protein list:"):
                    diseases = parse_protein_list(line.split(":")[1])

    else:

        input_proteins = proteins.split(" ")
        res = []
        for protein in input_proteins:
            url = "https://rest.uniprot.org/uniprotkb/" + protein + ".txt"
            r = requests.head(url)
            if r.status_code == 200:
                res.append(protein)
                print("The protein: " + protein)
                pass
            else:
                print("The wrong url is: " + url)
                raise Exception("Invalid UniProtKB accession %s"%protein)
        return res


def parse_include_synonyms(synonyms): #TODO move to utils
    return synonyms.lower().startswith('y') or synonyms.lower().startswith('t')


def parse_output_folder(folder):
    output_file = Path(folder)
    output_file.parent.mkdir(exist_ok=True, parents=True)
    return folder
    # output_file.write_text("output")


def load_mesh_terms(output_folder,file_to_link_file,debug=False):

    if debug:
        print("Importing list of valid MeSH terms")
    data_folder = os.path.join(output_folder,"data")
    
    # Check if mtrees file exists, otherwise download
    input_file = os.path.join(data_folder, 'MeSH/mtrees2023.bin') # TODO do not hard code
    file_exists = check_file((input_file))
    if not file_exists:
        if not os.path.exists(data_folder):
            os.makedirs(data_folder)
        if not os.path.exists(os.path.join(data_folder,'MeSH')):
            os.makedirs(os.path.join(data_folder,'MeSH'))
        ret = download_data({'MeSH':{'mtrees2023.bin':False}}, data_folder, file_to_link_file)

    # parse list of valid mesh terms
    valid_mesh_terms = set([l.strip("\n").split(";")[1] for l in open(input_file,"r").readlines()])
    print("%d MeSH terms parsed"%len(valid_mesh_terms))
    return valid_mesh_terms


def load_go_terms(output_folder, file_to_link_file, debug=False):
    if debug:
        print("Importing list of valid GO-terms")
    # Check if go term file exists, otherwise download
    data_folder = os.path.join(output_folder,"data")
    input_file = os.path.join(data_folder, 'GO/go-basic.obo')
    file_exists = check_file(input_file)
    if not file_exists:
        if not os.path.exists(data_folder):
            os.makedirs(data_folder)
        if not os.path.exists(os.path.join(data_folder,'GO')):
            os.makedirs(os.path.join(data_folder,'GO'))
        ret = download_data({'GO':{'go-basic.obo':False}}, data_folder, file_to_link_file)

    # parse list of valid go terms
    GO_TAG = 'id: '
    valid_go_terms = set([l.strip("\n").strip(GO_TAG) for l in open(input_file,"r").readlines() if GO_TAG in l])
    print("%d GO terms parsed"%len(valid_go_terms))
    return valid_go_terms


def save_parameters(output_folder, parameters,debug=False):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    output_file = os.path.join(output_folder,"parameters.json")
    json.dump(parameters, open(output_file, 'w'))
    if debug:
        print("Saved parameters as %s"%output_file)


def copy_folder(input_folder, output_folder):

    # create the output folder if it doesn't exist
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # copy all items in input folder
    for item in os.listdir(input_folder):
        item_path = os.path.join(input_folder, item)
        # copy if it is a file
        if os.path.isfile(item_path):
            shutil.copy(item_path,output_folder)
        # recursive copy if subfolder
        elif os.path.isdir(item_path):
            subfolder_output = os.path.join(output_folder, item)
            copy_folder(item_path, subfolder_output)

def preprocessing(args, debug=False):

    print("Preproccessing branch")

    # check if output folder is valid, otherwise create the folder
    output_folder = parse_output_folder(args.output_folder)
    # prepare data folders
    data_folder = os.path.join(output_folder,'data')
    mapping_folder = os.path.join(output_folder,'parsed_mappings')
    analysis_output_folder = os.path.join(output_folder,'output')
    config_folder = os.path.join(output_folder,'config')
    file_to_link_file = os.path.join(config_folder,'knowledge_base_links.json')
    # make directories if not exist
    for d in [data_folder, mapping_folder, analysis_output_folder]:
        if not os.path.exists(d):
            os.makedirs(d)
    if not os.path.exists(config_folder):
        original_config = '/workspace/caseolap_lift/config'
        copy_folder(original_config, config_folder)


    # import list of valid mesh and go terms to check user input against
    valid_mesh_terms = load_mesh_terms(output_folder,file_to_link_file, debug=debug)
    valid_go_terms = load_go_terms(output_folder,file_to_link_file, debug=debug)

    # Check that all required arguments are provided
    # At minimum, user must provide a parameter file OR ((protein_list OR GO-term) AND (disease_list))
    has_parameter_file = (args.parameters is not None)
    if not has_parameter_file:
        has_protein_list = (args.protein_list is not None)
        has_cellular_component = (args.subcellular_component is not None)
        if has_protein_list or has_cellular_component:
            pass
        else:
            print("Error! Please provide a custom list of proteins OR a GO term of interest")
            raise Exception("No protein list or go-term provided")
        has_disease_list = (args.disease_list is not None)
        if not has_disease_list:
            print("Error! Please provide a disease list of interest")
            raise Exception("No disease list provided")

    if has_parameter_file:
        param_file_name = args.parameters
        parameters = parse_preprocessing_parameters_file(param_file_name, valid_mesh_terms, valid_go_terms)
    else:

        # parse input files
        diseases = parse_diseases(args.disease_list, valid_mesh_terms)
        abbreviations = parse_abbreviations(args.abbreviations, len(diseases), debug=debug)
        cellular_component = parse_subcellular_component(args.subcellular_component, valid_go_terms)
        protein_list = parse_protein_list(args.protein_list)
        include_synonyms = getattr(args,'include-synonyms')
        include_ppi = getattr(args,'include-ppi')
        ppi_k = args.ppi_k
        ppi_thresh = args.ppi_score_thresh
        include_reactome = getattr(args, 'include-pw')
        pathway_count_thresh = args.pathway_count_thresh
        pathway_prop_thresh = args.pathway_prop_thresh
        include_tfd = getattr(args,'include-tfd')
        filter_against_proteome = getattr(args,'include-filter-against-proteome')
        redownload = args.redownload
        remap = args.remap
        if redownload:
            remap = True

        parameters = {'disease_categories':diseases,
                   'abbreviations':abbreviations,
                   'cellular_components':cellular_component,
                   'protein_list':protein_list,
                   'include_synonyms':include_synonyms,
                   'include_ppi':include_ppi,
                   'ppi_k':ppi_k,
                   'ppi_thresh':ppi_thresh,
                   'include_reactome':include_reactome,
                   'pathway_count_thresh':pathway_count_thresh,
                   'pathway_prop_thresh':pathway_prop_thresh,
                   'include_tfd':include_tfd,
                   'filter_against_proteome':filter_against_proteome,
                   'redownload':redownload,
                   'remap':remap
                   }



    # save parameters as json in output folder
    save_parameters(analysis_output_folder, parameters, debug=True)
    
    # Run the proprocessing module
    successful = prepare_knowledge_base_data(data_folder, mapping_folder, file_to_link_file,
                                             include_reactome=parameters['include_reactome'],
                                             include_tfd=parameters['include_tfd'],
                                             redownload=parameters['redownload'],
                                             remap=parameters['remap'], debug=debug)

    if successful:
        print("Knowledge base data successfully downloaded and mapped.")
    else:
        print("Problem with downloading knowledge base data!")
        sys.exit(1)
    # Run entity_set_expansion
    ent_parameters = {'go-term': parameters['cellular_components'],
                        'include_ppi': parameters['include_ppi'],
                        'ppi_k':  parameters['ppi_k'],
                        'ppi_score_thresh':  parameters['ppi_thresh'],
                        'include_pathways':  parameters['include_reactome'],
                        'pw_count_thresh':  parameters['pathway_count_thresh'],
                        'pw_proportion_thresh': parameters['pathway_prop_thresh'],
                        'include_transcription_factor_dependence': parameters['include_tfd'],
                        'filter_against_proteome': parameters['filter_against_proteome']
                      }
    prepare_subcellular_compartment_proteins(ent_parameters,
                                             mapping_folder=mapping_folder,
                                             data_folder=data_folder,
                                             output_folder=analysis_output_folder, debug=False)

    print("Done with preprocessing module.")
    categories = {k:v for k,v in zip(abbreviations,diseases)}
    with open(os.path.join(analysis_output_folder,'categories.txt'),'w') as file:
        file.write(json.dumps(categories,indent=4))

    return True


def check_date(date):
    # Checks to see if date a valid input.

    #string entered as YYYY-MM-DD
    year, month, day, = date.split('-')
    valid = True
    try:
        datetime.datetime(year=int(year), month=int(month), day=int(day))
    except ValueError:
        valid = False

    return valid


def parse_range(date_range):

    if date_range is None:
        return

    # account for a none value getting put in from the default
    # example date range: 2022-10-24,2025-11-23
    date_range = date_range.replace('"', '').replace('\'', '')
    date_from, date_to = date_range.split(",")

    if check_date(date_from) and check_date(date_to):
        return (date_from, date_to)
    else:
        raise Exception("The date range is invalid")


def parse_text(full_text):
    return full_text.lower().startswith('y') or full_text.lower().startswith('t')


def impute_label(labels):
    return labels.lower().startswith('y') or labels.lower().startswith('t')


def text_mining(args):
    print("Textmining branch")

    # check if output folder is valid, otherwise create the folder
    output_folder = parse_output_folder(args.output_folder)

    # parse parameters
    has_parameter_file = (args.parameters is not None)

    if has_parameter_file:
        param_file_name = args.parameters
        parameters = parse_text_mining_parameters_file(param_file_name)
    else:
        #check if the date is valid
        date_range = parse_range(args.date_range)

        include_full_text = getattr(args,'include-full-text')

        include_label_imputation = getattr(args, 'include-label-imputation')
        # impute_labels = args.impute_labels

        check_synonyms = args.check_synonyms

        rerun_scoring = args.rerun_scoring

        parameters = { 'date_range':date_range,
                       'include_full_text':include_full_text,
                       'include_label_imputation':include_label_imputation,
                       'check_synonyms':check_synonyms,
                       'rerun_scoring':rerun_scoring
                   }
    print(parameters)

    # prepare data folders
    data_folder = os.path.join(output_folder,'data')
    mapping_folder = os.path.join(output_folder,'parsed_mappings')
    analysis_output_folder = os.path.join(output_folder,'output')

    # save parameters as json in output folder
    save_parameters(analysis_output_folder, parameters,debug=True)

    # Run the text mining modules
    run_text_mining(output_folder, data_folder, mapping_folder, analysis_output_folder,
                    date_range=parameters['date_range'],
                    include_full_text=parameters['include_full_text'],
                    include_label_imputation=parameters['include_label_imputation'],
                    check_synonyms=parameters['check_synonyms'],
                    rerun_scoring=parameters['rerun_scoring'])

    print("Done with text mining module.")
    return True


def analysis(args):
    print("Analysis branch, still in development")

    # check if output folder is valid, otherwise create the folder
    output_folder = parse_output_folder(args.output_folder)

    # parse parameters
    has_parameter_file = (args.parameters is not None)

    if has_parameter_file:
        param_file_name = args.parameters
        parameters = parse_analyze_results(param_file_name)
    else:
        analyze_all_proteins = args.analyze_all_proteins
        z_score_thresh = args.z_score_thresh

        parameters = { 'analyze_all_proteins':analyze_all_proteins,
                       'z_score_thresh':z_score_thresh
        }
    print(parameters)

    # prepare data folders
    data_folder = os.path.join(output_folder, 'data')
    mapping_folder = os.path.join(output_folder, 'parsed_mappings')
    analysis_output_folder = os.path.join(output_folder, 'output')

    # save parameters as json in output folder
    save_parameters(analysis_output_folder, parameters, debug=True)

    # Run analyze_results module
    analyze_results(output_folder, z_score_thresh = parameters['z_score_thresh'],
                    merge_proteins=True,
                    use_core_proteins=(not parameters['analyze_all_proteins']),
                    debug=True)

    print("Done with analysis module.")
    return True


def prepare_kg(args):
    print("Prepare knowledge graph branch, still in development")

    # check if output folder is valid, otherwise create the folder
    output_folder = parse_output_folder(args.output_folder)

    # parse parameters
    has_parameter_file = (args.parameters is not None)

    if has_parameter_file:
        param_file_name = args.parameters
        parameters = parse_analyze_results(param_file_name)
    else:
    #    caseolap_score_type = args.caseolap_scores
        caseolap_score_type = 'raw'
        if args.use_z_score:
            caseolap_score_type = 'z_score'
        elif args.scale_z_score:
            caseolap_score_type = 'scaled_z_score'
    #prepare_knowledge_graph.add_argument('--caseolap_scores', default='z_score', choices=['raw','z_score','scaled_z_score'],
        
        include_all_proteins = args.include_all_proteins
        include_mesh = getattr(args, 'include-mesh')
        include_ppi = getattr(args, 'include-ppi')
        include_pw = getattr(args, 'include-pw')
        include_tfd = getattr(args, 'include-tfd')

        parameters = { "include_all_proteins":include_all_proteins,
                       "include_mesh": include_mesh,
                       "include_ppi": include_ppi,
                       "include_pw": include_pw,
                       "include_tfd": include_tfd,
                       "caseolap_score_type":caseolap_score_type
                       }
    print(parameters)

    # prepare data folders
    data_folder = os.path.join(output_folder, 'data')
    mapping_folder = os.path.join(output_folder, 'parsed_mappings')
    analysis_output_folder = os.path.join(output_folder, 'result')

    # save parameters as json in output folder
    save_parameters(analysis_output_folder, parameters, debug=True)

    caseolapLIFT = caseolapLIFT_knowledge_graph(parameters['include_ppi'],
                                                parameters['include_pw'],
                                                parameters['include_mesh'],
                                                parameters['include_tfd'],
                                                parameters['caseolap_score_type'],
                                                root_directory=output_folder,
                                                output_directory=analysis_output_folder)
    # graph = caseolapLIFT.knowledge_graph()

    print("Done with prepare knowledge graph module.")
    return True

def check_file(file_name):
    # Checks to see if input file is valid.
    return os.path.isfile(file_name)


def parse_preprocessing_parameters_file(file_name, valid_mesh_terms, valid_go_terms):
    # User provided input as a parameter file.
    is_valid = check_file(file_name)

    if is_valid:
        print("Parsing parameter file")
        # TODO read parameters (as a json or a text field, TBD).
        parameters = {}
        num_categories = 0
        with open(file_name) as f:
            lines = f.readlines()

        # read through parameters file to find out how many disease categories there are
        for line in lines:
            if line.startswith("disease_categories:"):
                num_categories += 1

        include_ppi = False
        ppi_k = False
        ppi_score_thresh = False
        include_reactome = False
        pathway_count_thresh = False
        pathway_prop_thresh = False
        include_tfd = False

        for line in lines:
            line = line.replace("\n", "")

            if line.startswith("disease_categories: "):
                dis = parse_diseases(line.split(": ")[1], valid_mesh_terms)

            # written for 1 line of abbreviations
            if line.startswith("abbreviations: "):
                # if they don't give abbreivations, add default abbrevs (mesh codes from below (list of 8))
                abbrev = parse_abbreviations(line.split(": ")[1], num_categories)
            # written for 1 line of cellular components right now, can change later
            if line.startswith("cellular_components: "):
                # required, unless they give a protein list
                # no default list
                # have to have one or the other, otherwise throw an error
                comp = parse_subcellular_component((line.split(": ")[1]), valid_go_terms)
            # written for 1 line of proteins
            if line.startswith("protein_list: "):
                # if they give cell comp then don't need it
                # no default list
                protein_string = line.split(": ")[1]
                if protein_string != "None":
                    proteins = parse_protein_list()
                else:
                    protein = None
            if line.startswith("include_synonyms: "):
                # maake the default true if they don't give anythign
                synonyms = parse_include_synonyms(line.split(": ")[1])
            if line.startswith("include_ppi: "):
                include_ppi = True
            if line.startswith("ppi_k: "):
                ppi_k = True
            if line.startswith("ppi_thresh: "):
                ppi_score_thresh = True
            if line.startswith("include_reactome: "):
                include_reactome = True
            if line.startswith("pathway_prop_thresh: "):
                pathway_prop_thresh = True
            if line.startswith("pathway_count_thresh: "):
                pathway_count_thresh = True
            if line.startswith("include_tfd: "):
                include_tfd = True
            # include parsing things for textmining as well

        # add the respective lists into parameters
        # have

        with open(file_name, 'a') as f:
            if not include_ppi:
                f.write("include_ppi: False")
            if not ppi_k:
                f.write("ppi_k: 1")
            if not ppi_score_thresh:
                f.write("ppi_thresh: 0.9")
            if not include_reactome:
                f.write("include_reactome: False")
            if not pathway_count_thresh:
                f.write("pathway_count_thresh: 4")
            if not pathway_prop_thresh:
                f.write("pathway_prop_thresh: 0.5")
            if not include_tfd:
                f.write("include_tfd: False")
        f.close()

        # Convert into json
        converter = FileConversion(file_name, "outputfile.json")
        x = converter.txt_to_json()
        y = json.dumps(x)
        return x
    else:
        raise Exception("The parameter file is invalid")


def parse_analyze_results(file_name):
    is_valid = check_file(file_name)

    if is_valid:
        print('Parsing analyze results parameters file')
        with open(file_name) as f:
            lines = f.readlines()

        include_z_score_thresh = False
        include_analyze_all_proteins = False
        include_analyze_core_proteins = False

        for line in lines:
            line = line.replace("\n", "")

            if line.startswith("z_score_thresh: "):
                include_z_score_thresh = True
            if line.startswith("analyze_all_proteins: "):
                include_analyze_all_proteins = True
            if line.startswith("analyze_core_proteins: "):
                include_analyze_core_proteins = True

        with open(file_name, 'a') as f:
            if not include_z_score_thresh:
                f.write("z_score_thresh: 3.0\n")
            if not include_analyze_all_proteins:
                f.write("analyze_all_proteins: False\n")
            if not include_analyze_core_proteins:
                f.write("analyze_core_proteins: False\n")
        f.close()

        converter = FileConversion(file_name, "outputfile.json")
        print(converter.txt_to_json())
        return converter.txt_to_json()
    else:
        raise Exception("Invalid file")

def parse_text_mining_parameters_file(file_name):
    is_valid = check_file(file_name)
    if is_valid:
        print("Parsing text mining parameter file")
        with open(file_name) as f:
            lines = f.readlines()
        date_range = False
        full_text = False
        impute_labels = False
        check_synonyms = False
        rerun_scoring = False
        include_all_proteins = False
        for line in lines:
            #line = line.replace("\n", "")
            if line.startswith("date_start "):
                range_string = line.split(": ")[1]
                if range_string != "None":
                    date_range = True
                    if parse_range(range_string):
                        date_range = True
                else:
                    date_range = False
                    range = None
            if line.startswith("full_text: "):
                full_text = True
            if line.startswith("impute_labels: "):
                impute_labels = True
            if line.startswith("check_synonyms: "):
                check_synonyms = True
            if line.startswith("rerun_scoring: "):
                rerun_scoring = True
        with open(file_name, 'a') as f:
            if not date_range:
                f.write("date_start date_end: None")
            if not full_text:
                f.write("full_text: False")
            if not impute_labels:
                f.write("impute labels: False")
            if not check_synonyms:
                f.write("check_synonyms: False")
            if not rerun_scoring:
                f.write("rerun_scoring: False")
        f.close()
        converter = FileConversion(file_name, "outputfile.json")
        x = converter.txt_to_json()
        y = json.dumps(x)
        return x
    else:
        raise Exception("The parameter file is invalid")

def parse_prepare_knowledge_graph_file(file_name):
    is_valid = check_file(file_name)
    if is_valid:
        print("Parsing prepare knowledge graph parameter file")
        with open(file_name) as f:
            lines = f.readlines()
        include_tfd = False
        include_ppi = False
        include_pw = False
        include_mesh = False
        use_z_score = False
        scale_z_score = False
        include_all_proteins = False
        include_core_proteins = False
        for line in lines:
            if line.startswith("include_tfd: "):
                include_tfd = True
            if line.startswith("include_ppi: "):
                include_ppi = True
            if line.startswith("include_pw: "):
                include_pw = True
            if line.startswith("include_mesh: "):
                include_mesh = True
            # if line.startswith("use_z_score: "):
            #     use_z_score = True
            # if line.startswith("scale_z_score: "):
            #     scale_z_score = True
            if line.startswith("include_all_proteins: "):
                include_all_proteins = True
            if line.startswith("include_core_proteins: "):
                include_core_proteins = True
        with open(file_name, 'a') as f:
            if not include_tfd:
                f.write("include_tfd: True")
            if not include_ppi:
                f.write("include_ppi: True")
            if not include_pw:
                f.write("include_pw: True")
            if not include_mesh:
                f.write("include_mesh: True")
            if not include_all_proteins:
                f.write("include_all_proteins: True")
            if not include_core_proteins:
                f.write("include_core_proteins: False")
        f.close()
        converter = FileConversion(file_name, "outputfile.json")
        x = converter.txt_to_json()
        y = json.dumps(x)
        return x
    else:
        raise Exception("The parameter file is invalid")


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def add_bool_arg(parser, name, default=False, single_flag=None, help=''):
    ''' source: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse'''
    help_str_feature = help+' Default: %s' % (str(default))
    help_str_no_feature = help+' Default: %s'%(str(not default))
    group = parser.add_mutually_exclusive_group(required=False)
    if single_flag is not None:
        group.add_argument(single_flag,'--' + name, dest=name, action='store_true', help=help_str_feature)
    else:
        group.add_argument('--' + name, dest=name, action='store_true', help=help_str_feature)
    group.add_argument('--no-' + name.replace('include-',''), dest=name, action='store_false',help=help_str_no_feature)
    parser.set_defaults(**{name:default})


def args_parser():
    '''
    Make subparsers for different components of the pipeline
    source: https://towardsdatascience.com/a-simple-guide-to-command-line-arguments-with-argparse-6824c30ab1c3
    '''

    # each will get their own subparser
    # Create a description for the CLI
    DESCRIPTION = ("There are three components to the CaseOLAP LIFT pipeline which must be performed in order\n"
                   "preprocessing      Prepare diseases and proteins of interest\n"
                   "text_mining        Perform text mining analysis\n"
                   "kg_analysis        Analyze the knowledge graph help\n")
    default_output_folder = '/caseolap_lift_shared_folder'
    if not(os.path.exists(default_output_folder) and os.path.isdir(default_output_folder)):
        default_output_folder = '.'

    # Create parser object
    parser = MyParser()

    # adding subparsers
    subparser = parser.add_subparsers(dest="command")
    preprocessing = subparser.add_parser("preprocessing")
    text_mining = subparser.add_parser("text_mining")
    analyze_results = subparser.add_parser("analyze_results")
    prepare_knowledge_graph = subparser.add_parser("prepare_knowledge_graph")

    # Add preprocess flags
    preprocessing.add_argument('-a', '--abbreviations', type=str, required=False,
                               help="Abbreviations for diseases")
    preprocessing.add_argument('-d', '--disease_list', type=str, required=False,
                               help="List of diseases of interest")
    preprocessing.add_argument('-c', '--subcellular_component', type=str, required=False,
                               help="Subcellular component of interest")
    preprocessing.add_argument('-l', '--protein_list', type=str, required=False,
                               help="Include a custom list of proteins of interest")
    preprocessing.add_argument('-o', '--output_folder', type=str, required=False,
                               help='Output directory to program results')
    preprocessing.add_argument('-p', '--parameters', type=str, required=False,
                               help='Input .json or .txt file specifying parameters to run')
    add_bool_arg(preprocessing, 'include-synonyms', default=True,
                 help='Include UniProt synonyms of proteins with text mining step')
    add_bool_arg(preprocessing, 'include-ppi', default=False,
                 help='Include proteins with STRING protein-protein interactions with entity set expansion')
    preprocessing.add_argument('-k', '--ppi_k', type=int, required=False,
                               help='k-hop neighbors of STRING PPI to include. Default:1')
    preprocessing.add_argument('-s', '--ppi_score_thresh', type=float, required=False,
                               help='STRING score threshold, only include interactors above this threshold. Default: 0.9')
    add_bool_arg(preprocessing, 'include-pw', default=False,
                 help='Include proteins with shared Reactome pathways with entity set expansion')
    preprocessing.add_argument('-n', '--pathway_count_thresh', type=int, required=False,
                               help='Minimum number of subcellular component proteins required to consider a pathway as significant. Default:4')
    preprocessing.add_argument('-r', '--pathway_prop_thresh', type=float, required=False,
                               help='Minimum proportion of subcellular component proteins required to consider a pathway as significant. Default: 0.5')
    preprocessing.add_argument('--redownload', action='store_true', required=False, default=False,
                               help="Redownload knowledge base information")
    preprocessing.add_argument('--remap', action='store_true', required=False, default=False,
                               help="Re-calculate parsed mappings")
    add_bool_arg(preprocessing, 'include-tfd', default=False,
                 help='Include proteins with transcription factor dependence from GRNdb with entity set expansion')
    add_bool_arg(preprocessing, 'include-filter-against-proteome', default=True,
                 help='Whether to filter UniProt IDs against the UniProt Human reference proteome.')
    preprocessing.set_defaults(output_folder=default_output_folder)
    preprocessing.set_defaults(ppi_k=1)
    preprocessing.set_defaults(ppi_score_thresh=0.9)
    preprocessing.set_defaults(pathway_count_thresh=4)
    preprocessing.set_defaults(pathway_prop_thresh=0.5)

    text_mining.add_argument('-d', '--date_range', type=str, required=False,
                             default=None, help='Specify the date range for PubMed documents which will be downloaded')
    add_bool_arg(text_mining, 'include-full-text', default=False, single_flag='-t',
                 help='Specify to use full-text in text mining analysis or not')
    add_bool_arg(text_mining, 'include-label-imputation', default=False, single_flag='-l',
                 help='Whether to impute missing labels on text')
    text_mining.add_argument('-c', '--check_synonyms', type=bool, required=False, default=False,
                             help='Indicating the software to halt function to screen ambiguous protein names')
    text_mining.add_argument('-r', '--rerun_scoring', type=bool, required=False, default=False,
                             help='Indicating the software to re-run CaseOLAP score calculation after updating synonym list')
    text_mining.add_argument('-o', '--output_folder', type=str, required=False,
                             help='Directory where all data and results will be stored. Please make sure to use the same output folder for all steps, as future steps rely on files output from previous steps')
    text_mining.add_argument('-p', '--parameters', type=str, required=False,
                             help='Will bypass all options and run based on parameters.txt or parameters.json file')
    text_mining.set_defaults(output_folder=default_output_folder)

    analyze_results.add_argument('-z', '--z_score_thresh', type=float, required=False, default=3.0,
                                 help='The threshold to consider a protein as obtaining a significant score with respect to a disease category')
    analyze_protein_group = analyze_results.add_mutually_exclusive_group(required=False)
    analyze_protein_group.add_argument('--analyze_all_proteins', dest='analyze_all_proteins', action='store_true',
                                       help='Analyze CaseOLAP results from all functionally related proteins. Default: False')
    analyze_protein_group.add_argument('--analyze_core_proteins', dest='analyze_all_proteins', action='store_false',
                                       help='Analyze CaseOLAP results from only core proteins related to the cellular component. Default: True')
    analyze_results.set_defaults(**{'analyze_all_proteins': False})
    analyze_results.add_argument('-o', '--output_folder', type=str, required=False,
                                 help='Directory where all data and results will be stored. Please make sure to use the same output folder for all steps, as future steps rely on files output from previous steps')
    analyze_results.add_argument('-p', '--parameters', type=str,
                                 help='Will bypass all options and run based on parameters.txt or parameters.json file')
    analyze_results.set_defaults(output_folder=default_output_folder)



    add_bool_arg(prepare_knowledge_graph, 'include-tfd', default=True,
                 help='Indicating the software to include edges between proteins with transcription factor dependence')
    add_bool_arg(prepare_knowledge_graph, 'include-ppi', default=True,
                 help='Indicating the software to include edges between proteins with protein-protein interaction')
    add_bool_arg(prepare_knowledge_graph, 'include-pw', default=True,
                 help='Indicating the software to include Reactome pathway nodes and edges between proteins and pathways')
    add_bool_arg(prepare_knowledge_graph, 'include-mesh', default=True,
                 help='Indicating the software to include the MeSH disease hierarchy')
    add_bool_arg(prepare_knowledge_graph, 'use_z_score', default=False,
                 help='Indicating the software to include z-score scale CaseOLAP scores')
    add_bool_arg(prepare_knowledge_graph, 'scale_z_score', default=True,
                 help='Indicating the software to scale the z-score transformed CaseOLAP scores')
    # which caseolap scores to include
    #prepare_knowledge_graph.add_argument('--caseolap_scores', default='z_score', choices=['raw','z_score','scaled_z_score'],
    #                                     help='Specifies which CaseOLAP scores to include in the knowledge graph. Raw scores (between 0.0 and 1.0), z-score per disease category (mean 0), or positive scaled z-scores (all scores > 0.0)')
    # mutual exclusive group on include_all_proteins
    kg_protein_group = prepare_knowledge_graph.add_mutually_exclusive_group(required=False)

    kg_protein_group.add_argument('--include_all_proteins', dest='include_all_proteins', action='store_true',
                                       help='Analyze CaseOLAP results from all functionally related proteins. Default: False')
    kg_protein_group.add_argument('--include_core_proteins', dest='include_all_proteins', action='store_false',
                                       help='Analyze CaseOLAP results from only core proteins related to the cellular component. Default: True')
    prepare_knowledge_graph.set_defaults(**{'include_all_proteins': False})
    prepare_knowledge_graph.add_argument('-o', '--output_folder', type=str,
                                         help='Directory where all data and results will be stored. Please make sure to use the same output folder for all steps, as future steps rely on files output from previous steps')
    prepare_knowledge_graph.add_argument('-p', '--parameters', type=str,
                                         help='Will bypass all options and run based on parameters.txt or parameters.json file')
    prepare_knowledge_graph.set_defaults(output_folder=default_output_folder)

    return parser


def main():

    # Set up argument parsing
    parser = args_parser()
    args = parser.parse_args()
    print(args)
    
    # Print help if no arguments given 
    if len(sys.argv) < 2:
        parser.error("Incorrect usage.")
        sys.exit(0)
    
    # Execute sub-branches of program depending on command
    command = args.command 
    if command == 'preprocessing':
        preprocessing(args)
    elif command == 'text_mining':
        text_mining(args)
    elif command == 'analyze_results':
        analysis(args)
    elif command == 'prepare_knowledge_graph':
        prepare_kg(args)
    else:
        parser.error("Mode not found: %s" % sys.argv)
        sys.exit(0)


if __name__ == "__main__":
    main()
