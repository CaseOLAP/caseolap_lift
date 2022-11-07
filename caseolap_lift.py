
import argparse
import sys
import os
import datetime
import requests
from pathlib import Path
from preprocessing.prepare_knowledge_base_data import  download_data,prepare_knowledge_base_data
import json
from preprocessing.entity_set_expansion import prepare_subcellular_compartment_proteins


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
            # TODO change it to be from categories.txt not parameter.txt
            # TODO there will be no disease category heading, jst straight dieases ids
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


def parse_include_synonyms(synonyms): #TODO remove
    return synonyms.lower().startswith('y') or synonyms.lower().startswith('t')


def parse_output_folder(folder):
    output_file = Path(folder)
    output_file.parent.mkdir(exist_ok=True, parents=True)
    return folder
    # output_file.write_text("output")


def load_mesh_terms(output_folder,debug=False):
    if debug:
        print("Importing list of valid MeSH terms")

    # Check if mtrees file exists, otherwise download
    data_folder = os.path.join(output_folder,"data")
    input_file = os.path.join(data_folder, 'MeSH/mtrees2021.bin')
    file_exists = check_file((input_file))
    if not file_exists:
        ret = download_data({'MeSH':{'mtrees2021.bin':False}}, data_folder)

    # parse list of valid mesh terms
    valid_mesh_terms = set([l.strip("\n").split(";")[1] for l in open(input_file,"r").readlines()])
    print("%d MeSH terms parsed"%len(valid_mesh_terms))
    return valid_mesh_terms


def load_go_terms(output_folder, debug=False):
    if debug:
        print("Importing list of valid GO-terms")
    # Check if go term file exists, otherwise download
    data_folder = os.path.join(output_folder,"data")
    input_file = os.path.join(data_folder, 'GO/go-basic.obo')
    file_exists = check_file(input_file)
    if not file_exists:
        ret = download_data({'GO':{'go-basic.obo':False}}, data_folder)

    # parse list of valid go terms
    GO_TAG = 'id: '
    valid_go_terms = set([l.strip("\n").strip(GO_TAG) for l in open(input_file,"r").readlines() if GO_TAG in l])
    print("%d GO terms parsed"%len(valid_go_terms))
    return valid_go_terms


def save_parameters(output_folder, parameters,debug=False):
    output_file = os.path.join(output_folder,"parameters.json")
    json.dump(parameters, open(output_file, 'w'))
    if debug:
        print("Saved parameters as %s"%output_file)


def preprocessing(args, debug=False):

    print("Preproccessing branch")

    # check if output folder is valid, otherwise create the folder
    output_folder = parse_output_folder(args.output_folder)

    # import list of valid mesh and go terms to check user input against
    valid_mesh_terms = load_mesh_terms(output_folder, debug=debug)
    valid_go_terms = load_go_terms(output_folder, debug=debug)

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
        parameters = parse_parameters_file(param_file_name)
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
                   'include_tfd':include_tfd
                   }

    # prepare data folders
    data_folder = os.path.join(output_folder,'data')
    mapping_folder = os.path.join(output_folder,'parsed_mappings')
    analysis_output_folder = os.path.join(output_folder,'output')

    # save parameters as json in output folder
    save_parameters(analysis_output_folder, parameters,debug=True)

    # Run the proprocessing module
    successful = prepare_knowledge_base_data(data_folder, mapping_folder,
                                             include_reactome=parameters['include_reactome'],
                                             include_tfd=parameters['include_tfd'],
                                             redownload=False, debug=debug)

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
                        'filter_against_proteome': True #TODO
                      }
    prepare_subcellular_compartment_proteins(ent_parameters,
                                             mapping_folder=mapping_folder,
                                             output_folder=analysis_output_folder, debug=False)

    print("Done with preprocessing module.")
    return True


def check_date(date):
    # Checks to see if date a valid input.
    # TODO (completed?)

    #string entered as mm-dd-yyyy or mm/dd/yy
    month, day, year = date.split('-')

    valid = True
    try:
        datetime.datetime(int(year), int(month), int(day))
    except ValueError:
        valid = False

    return valid


def parse_range(date_range):
    # account for a none value getting put in from the default
    # example date range: 10-24-2022 11-23-2025
    date_from, date_to = date_range.split(" ")

    if check_date(date_from) and check_date(date_to):
        return (date_from, date_to)
    else:
        raise Exception("The date range is invalid")


def parse_text(full_text):
    return full_text.lower().startswith('y') or full_text.lower().startswith('t')


def impute_label(labels):
    return labels.lower().startswith('y') or labels.lower().startswith('t')


def textmining(args):
    print(args)
    print("Textmining branch, still in development")

    # flags: date range, full-text flag (boolean), impute-labels (boolean), parameters file, 
    # TODO add: input-folder from previous step; OR list of proteins and diseases from prev step.

    #date range: default is all pubs - represent it as none 
    ##full-text: default is false 
    #impute-

    #parameters file: use other function, specify which branch it came from
    #searching for a specific date range, use full text where available, 
    #just parsing the input from the user 

    date_range = parse_range(args.date_range)
    #check if the date is valid

    
    full_text = parse_text(args.full_text)

    impute_labels = impute_label(args.impute_labels)

    #take parameter files for all 3
    #an option to run the entire pipeline
    #need an example 

    #currently in 2 dif txtfiles that have to be edited before runnng 
    #goal is to make the inpt easier and more intuitive

    #have it as a parameters textfle 

    # Check arguments
    ## TODO same as above

    date_from,date_to = args.date_range # TODO or something like that
    if check_date(date_from) and check_date(date_to):
        pass
    else:
        # throw exception
        print("Invalid date")

    return


#STILL LEFT TO DO 
def analysis(args):
    print(args)
    print("Analysis branch, still in development")
    # TODO add: path to input-folder from previous step; OR add knowledge graph edges
    return


def check_file(file_name):
    # Checks to see if input file is valid.
    return os.path.isfile(file_name)
    

def parse_parameters_file(file_name):
    # User provided input as a parameter file.
    is_valid = check_file(file_name)

    if is_valid:
        print("Parsing parameter file")
        #TODO read parameters (as a json or a text field, TBD).
        parameters = {}

        with open(file_name) as f:
            lines = f.readlines()
        for line in lines:
            if line.startswith("Abbreviations:"):
                #if they don't give abbreivations, add default abbrevs (mesh codes from below (list of 8))
                abbreviations = parse_abbreviations(line.split(":")[1])
            if line.startswith("Disease Category:"):
                diseases = parse_diseases(line.split(":")[1])
            if line.startswith("Cellular component:"):
                #required, unless they give a protein list
                #no default list
                #have to have one or the other, otherwise throw an error
                components = parse_subcellular_component(line.split(":")[1])
            if line.startswith("Protein list:"):
                #if they give cell comp then don't need it
                #no default list
                proteins = parse_protein_list(line.split(":")[1])
            if line.startswith("Include synonyms:"):
                #maake the default true if they don't give anythign
                synonyms = parse_include_synonyms(line.split(":")[1])
            #include parsing things for textmining as well 

        #add the respective lists into parameters
        #have
        return parameters
    else:
        #TODO throw an error, tell user the parameter file is invalid. (completed??)
        raise Exception("The parameter file is invalid")
        print("Error")


class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)


def add_bool_arg(parser, name, default=False, help=''):
    ''' source: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse'''
    help_str_feature = help+' Default: %s' % (str(default))
    help_str_no_feature = help+' Default: %s'%(str(not default))
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument('--' + name, dest=name, action='store_true', help=help_str_feature)
    group.add_argument('--no-' + name, dest=name, action='store_false',help=help_str_no_feature)
    parser.set_defaults(**{name:default})


def args_parser():

    '''
    Make subparsers for different components of the pipeline
    source: https://towardsdatascience.com/a-simple-guide-to-command-line-arguments-with-argparse-6824c30ab1c3
    '''
    
    #each will get their own subparser
    # Create a description for the CLI
    DESCRIPTION = ("There are three components to the CaseOLAP LIFT pipeline which must be performed in order\n" 
                    "preprocessing      Prepare diseases and proteins of interest\n" 
                    "text_mining        Perform text mining analysis\n" 
                    "kg_analysis        Analyze the knowledge graph help\n")

    # Create parser object
    parser = MyParser()

    #adding subparsers
    subparser = parser.add_subparsers(dest = "command")
    preprocessing = subparser.add_parser("preprocessing")
    text_mining = subparser.add_parser("text_mining")
    kg_analysis = subparser.add_parser("kg_analysis")

    # Add preprocess flags
    preprocessing.add_argument('-a', '--abbreviations', type=str, required=False,
                               help="Abbreviations for diseases")
    preprocessing.add_argument('-d', '--disease_list', type=str, required=False,
                               help = "List of diseases of interest")
    preprocessing.add_argument('-c', '--subcellular_component', type=str, required=False,
                               help="Subcellular component of interest")
    preprocessing.add_argument('-l', '--protein_list', type = str, required=False,
                               help="Include a custom list of proteins of interest")
    preprocessing.add_argument('-o', '--output_folder', type = str, required=False,
                               help='Output directory to program results')
    preprocessing.add_argument('-p', '--parameters', type=str, required=False,
                               help='Input .json or .txt file specifying parameters to run')
    add_bool_arg(preprocessing, 'include-synonyms',default=True,
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
    add_bool_arg(preprocessing, 'include-tfd', default=False, help = 'Include proteins with transcription factor dependence from GRNdb with entity set expansion')
    preprocessing.set_defaults(output_folder='.')
    preprocessing.set_defaults(ppi_k=1)
    preprocessing.set_defaults(ppi_score_thresh=0.9)
    preprocessing.set_defaults(pathway_count_thresh=4)
    preprocessing.set_defaults(pathway_prop_thresh=0.5)

    # Add text mining flags
    # add default values here
    text_mining.add_argument('-d date_start date_end', '--date_range date_start date_end', type=str, required=False,
                             default=None, help='Specify the date range for PubMed documents which will be downloaded')
    text_mining.add_argument('-f', '--full_text', type=str, required=False, default=False,
                             help='Specify to use full-text in text mining analysis or not')
    text_mining.add_argument('-i', '--impute_labels', type=str, required=False, default=False,
                             help = 'Whether to impute missing labels on text')

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
        textmining(args)
    elif command == 'analysis':
        analysis(args)
    else:
        parser.error("Mode not found: %s" % sys.argv)
        sys.exit(0)


if __name__ == "__main__":
    main()
