# CaseOLAP LIFT

This repository contains the code associated with the paper, “A knowledge graph approach to elucidate the role of organellar pathways in disease via biomedical reports". CaseOLAP LIFT is a new computation pipeline with novel functions an optimized existing workflow, to investigate disease and cellular component associations by extracting user-selected information from text datasets (e.g. biomedical literature). 

![workflow](assets/workflow.png)
![workflow](assets/diagram.png)



## Additional details from the protocol (link)

### Step 2: Preparing diseases and proteins
Knowledge base data and their corresponding file names are shown below.
GO Ontologies: go-basic.obo
GO to Proteins: goa_human.gaf
UniProt Human Proteome: UP000005640_9606.fasta
Proteins to Pathways (leaf node pathways in Reactome’s pathway hierarchy): UniProt2Reactome.txt
Proteins to Pathways (all pathways): UniProt2Reactome_All_Levels.txt
Pathways: ReactomePathwaysRelation.txt
Pathway names: ReactomePathways.txt
MeSH Term Data: desc2022.xml 
MeSH Tree Data: mtrees2022.bin
GRNdb/GTEx Transcription Factors: Files from GRNdb


‘python caseolap_lift.py preprocessing --help’. 

The following flags are required:

-d disease list: a list of MeSH tree numbers for your diseases of interest. Separate disease categories by spaces. Separate tree numbers of the same category by commas.

-c cellular component: a list of space-separated GO terms representing cellular component(s), molecular function(s), and/or biological process(es).
	The following flags are optional:

-a disease abbreviations: User-named abbreviations for disease categories of interest. Must match the number of disease categories passed in.

-l protein list: list of additional proteins to include in this analysis. Only UniProtKB IDs are accepted; unmapped IDs will be discarded. A list of proteins may be provided in place of GO term.

--include-ppi: Flag indicates proteins with known protein-protein interaction via STRING will be added to functionally related proteins.

-k ppi_k: indicates how many hops of protein-protein interactions to include. (e.g. k = 2 will include neighbors of neighbors) Default: 1.

-s ppi_score_thresh: only protein-protein interactions with a score greater than this threshold will be included.

--include-pw: Flag indicates proteins which share significant biological pathways via Reactome to be added to functionally related proteins.

-n pathway_count_thresh: minimum number of cellular component proteins required to consider a pathway as significant. Default: 4.

-r pathway_prop_thresh: minimum proportion of cellular component proteins required to consider a pathway as significant. Default: 0.5.

--include-tfd: Flag indicates that proteins related via transcription factor from GRNdb will be added to the functionally related proteins.

-o output_folder: directory where all data and results will be stored. Please make sure to use the same output folder for all steps, as future steps rely on files output from previous steps.

-p parameter_file: will bypass all options and run based on parameters.txt or parameters.json file.

