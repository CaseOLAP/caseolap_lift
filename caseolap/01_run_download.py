'''
The purpose of this file is to download the zipped files containing
the PubMed publications (i.e. documents). These will be mined later.
'''
import os, json, sys
from caseolap._01_download import *


'''
Parameters
'''
# Input
data_dir = './'
logFilePath = './log/download_log.txt'
download_config_file_path = './config/download_config.json'
ftp_config_file_path = './config/ftp_config.json'
baseline_dir = os.path.join(data_dir, 'ftp.ncbi.nlm.nih.gov/pubmed/baseline/')
update_files_dir = os.path.join(data_dir,'ftp.ncbi.nlm.nih.gov/pubmed/updatefiles/')



'''
Main Code
'''
# Start the download, verification, and extraction process
if __name__ == '__main__':
    
    # Open log file
    logfile = open(logFilePath, 'w') 
    
    # Load config files
    download_config = json.load(open(download_config_file_path))
    ftp_config = json.load(open(ftp_config_file_path))
    
    # Check main directory
    if not os.path.isdir(data_dir):
        print('Directory not found:', data_dir) 
    
    # Start download 
    download_pubmed(data_dir, download_config, ftp_config, logfile)
    
    # Verify download: 'baseline files' & 'update files'
    check_all_md5_in_dir(baseline_dir, logfile, linux = True, mac = False)
    check_all_md5_in_dir(update_files_dir, logfile, linux = True, mac = False)

    # Extract downloaded files: 'baseline files' & 'update files'
    extract_all_gz_in_dir(baseline_dir, logfile)
    extract_all_gz_in_dir(update_files_dir, logfile)
    
    logfile.close()
    
    
    
    