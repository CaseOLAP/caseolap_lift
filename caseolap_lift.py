
import argparse
from operator import truediv
import sys
import os
import datetime
from turtle import resetscreen
import requests
from pathlib import Path
from preprocessing.prepare_knowledge_base_data import  download_data,prepare_knowledge_base_data
from preprocessing.entity_set_expansion import prepare_subcellular_compartment_proteins

#want user input matching
#list of mesh codes and their synonyms, go cellular comps, find the closest match to what the user inputs
#if its too complicated, use parameters.txt 

def parse_abbreviations(input): 
    #user can write in any abbreviation that they want but it has to match the number of categgories
    #has to be unique abbreviations - if not unique, add a (1) at the end to make it unique
    correct_abrev = ['CM', 'ARR', 'CHD', 'VD', 'IHD', 'CCD', 'VOO', 'OTH']
    res = []
    abbreviations = input.split(' ')
    print("Abbreviation list : %s"%str(abbreviations))

    for abrev in abbreviations:
        if abrev in correct_abrev:
            res.append(abrev)
            print("The abbreviation is: " + abrev)
            pass
        else:
            print("The wrong abbreviation is: ")
            print(abrev)
            raise Exception("Invalid abbreviation")
    return res
#     #what are we checking the abbreviations against? 


def parse_diseases(diseases, valid_mesh_terms, debug=False):
#     #check if the user input are all valid MeSH tree ids
    if ".txt" in diseases:
        with open(diseases) as f:
            lines = f.readlines()
        for line in lines:
            #change it to be from categories.txt not parameter.txt 
            #there will be no disease category heading, jst straight dieases ids
            if line.startswith("Disease Category:"):
                    diseases = parse_diseases(line.split(":")[1], valid_mesh_terms)
    else:
        #if diseases 
        res = []
        #if the user does not give a text file
        #else: 
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
                # url = 'https://id.nlm.nih.gov/mesh/lookup/details?descriptor=' + disease
                # r = requests.get(url)
                # # print(r)
                # if r.status_code == 200:
                #     print("The mesh term: " + disease)
                #     pass
                # else:
                #     print("The wrong url is: " + url)
                #     raise Exception("Invalid mesh ID")
            res.append(diseases)
        return res

def parse_subcellular_component(components, valid_go_terms, debug=False):
#     #check if user inputs are valid GO terms
#http://api.geneontology.org/api/bioentity/function/GO:0006915
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

        # url = "http://api.geneontology.org/api/bioentity/function/" + component
        # r = requests.head(url)
        # #print()
        # if r.status_code == 200:
        #     res.append(component)
        #     print("The GO cellular component: " + component)
        #     pass
        # else:
        #     print("The wrong url is: " + url)
        #     raise Exception("Invalid GO ID")
    return res

def parse_protein_list(proteins):
#     #check if they are valid proteins with uniprot????
#COME BACK TO THIS 
#if it takes too long, try batching it or skip it for proteins (need it for diseases)
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

def parse_include_synonyms(synonyms):
    return synonyms.lower().startswith('y') or synonyms.lower().startswith('t')


def parse_output_folder(folder):
    output_file = Path(folder)
    output_file.parent.mkdir(exist_ok=True, parents=True)
    return folder
    # output_file.write_text("output")

def load_mesh_terms(output_folder):

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


def load_go_terms(output_folder):

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


def preprocessing(args, debug=False):
    print(args)
    print("Preproccessing branch")

    # check if output folder is valid, otherwise create the folder
    output_folder = parse_output_folder(args.output_folder)

    # import list of valid mesh and go terms to check user input against
    print("Importing list of valid MeSH terms")
    valid_mesh_terms = load_mesh_terms(output_folder)
    print("Importing list of valid GO-terms")
    valid_go_terms = load_go_terms(output_folder)
    
    ## TODO
    # flags: abbreviations, diseases, cellcomp, proteinlist, output folder, parameters, -s synonyms, 
    # Check arguments
    #TODO: check that all required argumetns are there 
    #eiether (cell comp or protein list), required: disease list
    #create a default output folder, same default for synonyms
    #figure out if we have parameter file 
    # check if args.parameters not none
    has_parameter_file = False #TODO
    if has_parameter_file:
        param_file_name = "TODO" #TODO
        parameters = parse_parameters_file(param_file_name)
    else:
        parameters = {}
        # TODO

        #check if you haev diseases, either cell comp or protein list
        #if they don't have it, raise exception and tell them the requirements

        # parse input files
        abbreviations = parse_abbreviations(args.abbreviations) #TODO make these functions, args.abbrev might not be right
        #valve disease = VD, ischemic heart disease = IHD 
        #ned to have a lit of mesh codes whch are given in as diseases (array of 8 dif mesh codes)
        #abrvs will be what we use to make it easier to understand 
        #take in as text file or read in as parameters
        #ex: (-d disease code 1, DC2, ...) input a list of mesh codes
        # TODO: how to parse the command line input (will it grab everything after -d and before next argument)
        #can also have it as an input textfile (-d diseases.txt(list of msh ids n a textfile)) and then parse th txtfile 
        diseases = parse_diseases(args.disease_list, valid_mesh_terms) # TODO same as above, also need to check against MeSH tree, if it is valid input
        cellular_component = parse_subcellular_component(args.subcellular_component, valid_go_terms) #TODO same, need to check against GO terms to see if it is valid input
        protein_list = parse_protein_list(args.protein_list) #TODO
        
        # TODO synonym list for proteins
        #make this true by default 
        include_synonyms = parse_include_synonyms(args.include_synonyms)


    #wrie a parameters.txt file tht will have the same output as the flag input
    print(parameters)
    # Run the proprocessing module
    print("TODO running preprocessing module")
    data_folder = os.path.join(output_folder,'data')
    mapping_folder = os.path.join(output_folder,'parsed_mappings')
    analysis_output_folder = os.path.join(output_folder,'output')
    #TODO handle variable downloads (i.e. if not using Reactome or TFD, don't need to download it)
    successful = prepare_knowledge_base_data(data_folder, mapping_folder, redownload=False, debug=debug)

    if successful:
        print("Knowledge base data successfully downloaded and mapped.")
    else:
        print("Problem with downloading knowledge base data!")
        sys.exit(1)

    # Run entity_set_expansion
    parameters = {'go-term': cellular_component,
                  'include_ppi': True,
                  'ppi_k': 1, 'ppi_score_thresh': 0.99,
                  'include_pathways': True,
                  'pw_count_thresh': 4,
                  'pw_proportion_thresh': 0.50,
                  'include_transcription_factor_dependence': True}
    prepare_subcellular_compartment_proteins(parameters, mapping_folder=mapping_folder, output_folder=analysis_output_folder, debug=False)

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
    #account for a none value getting put in from the default 
    #example date range: 10-24-2022 11-23-2025
    date_from, date_to = date_range.split(" ")

    if check_date(date_from) and check_date(date_to):
        return (date_from, date_to)
    else:
        raise Exception("The date range is invalid")

'''
For 10/25/2022:
preprocessing:
read in textfiles for both disease_list and protein_list
parse through the textfiles and return ids
synonym list- either return boolean or list of synonyms - DONE

textmining:
return some form of date range - DONE
return boolean for full text - DONE
return boolean for impute labels - DONE

parse_parameter_file - DONE
'''

def parse_text(full_text):
    return full_text.lower().startswith('y') or full_text.lower().startswith('t')

def impute_label(labels):
    return labels.lower().startswith('y') or labels.lower().startswith('t')


#pause on this right now
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
    # TODO (completed??)
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

def args_parser():
    ## TODO implement sub-parsers (completed??)
    
    '''
    TODO: implement sub-parsers. See 'login' vs 'register' in the link below. Three subparsers needed: 
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
    #parser = argparse.ArgumentParser(description = DESCRIPTION, formatter_class=argparse.RawTextHelpFormatter)
    

    #adding subparsers
    subparser = parser.add_subparsers(dest = "command")

    preprocessing = subparser.add_parser("preprocessing")
    text_mining = subparser.add_parser("text_mining")
    kg_analysis = subparser.add_parser("kg_analysis")

    #TODO: implement mutually exclusive arguments for all of these

    #check the logic 

    # Add preprocess flags
    preprocessing.add_argument('-a', '--abbreviations', type = str, required = False,help = "Abbreviations for diseases")
    #changed disease_list type to string instead of list
    preprocessing.add_argument('-d', '--disease_list', type = str, required = False, help = "List of diseases of interest")
    preprocessing.add_argument('-c', '--subcellular_component', type = str, required = False, help = "Subcellular component of interest")
    preprocessing.add_argument('-l', '--protein_list', type = str, required = False, help = "Include a custom list of proteins of interest")
    # unsure if this is right since there are a lot of spaces in the argument
    #changed this argument to include a full name
    preprocessing.add_argument('-o', '--output_folder', type = str, required = False, help = 'Output directory to program results') 
    preprocessing.add_argument('-p parameters parameters.txt', type = str, required = False, help = 'Input .json or .txt file specifying parameters to run')
    ##TODO actually, this should not have 'boolean', it should be a flag only. Same for all boolean (-f,-i)
    #removed boolean term from thoe flags. is that all that needs to be done?

    #by default is
    #test
    preprocessing.add_argument('-s', '--include_synonyms', type = str, required = False, help = 'Include synonyms for proteins')


    # Add text mining flags
    #add default values here
    text_mining.add_argument('-d date_start date_end', '--date_range date_start date_end', type = str, required = False, default = None, help = 'Specify the date range for PubMed documents which will be downloaded')
    text_mining.add_argument('-f', '--full_text', type = str, required = False, default = False, help = 'Specify to use full-text in text mining analysis or not')
    text_mining.add_argument('-i', '--impute_labels', type = str, required = False, default = False, help = 'Whether to impute missing labels on text')

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
    elif command =='analysis':
        analysis(args)
    else:
        parser.error("Mode not found: %s"%sys.argv)
        sys.exit(0)

if __name__ == "__main__":
    main()




#example parameters.txt file
