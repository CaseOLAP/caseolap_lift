import os, sys, re, time, subprocess

'''
Download
'''
def download_pubmed(data_dir, download_config, ftp_config, logfile):
    '''
    FUNCTION:
    - Downloads PubMed's 'baseline' files and 'update' files. Baseline
      files are from all previous years. Update files are from this year.
      
    PARAMS:
    - data_dir (str): The main directory which will contain folders for both
      the baseline and update files.
    - ftp_config (json): The ftp links to the baseline and update file downloads
    - download_config (json): Indicates whether to download baseline and/or 
      update. True/False  
    - logfile (file): The place to write the progress of downloading
    '''
     
        
    ''' Download Baseline'''
    if download_config['baseline']:

       # Time before download
        t1 = time.time()

        # Print progress
        msg = str('Downloading baseline files from '+ ftp_config['baseline'])
        print(msg +' For details, see the download logfile')
        logfile.write(msg +'\n')

        # Download files
        return_code = os.system('wget -N -q -r --directory-prefix='+data_dir+\
                       ' --no-parent '+ftp_config['baseline'])

        # Check for download error
        if return_code != 0:
            msg = str('Error (Baseline file download): Return code:'+str(return_code)+\
                  '\n Link: ' + ftp_config['baseline']+'\n')
            logfile.write(msg)
            exit(return_code)

        # Print progress 
        msg = str('Finished downloading PubMed baseline files. '+\
                  str(round(time.time() - t1))+' seconds')
        logfile.write(msg +'\n') 
        print(msg)



    ''' Download updatefiles'''
    if download_config['update']:

        # Time before download
        t2 = time.time()        

        # Print progress
        msg = str('Downloading update files from '+ ftp_config['update'])
        print(msg + ' For details, see the download logfile')
        logfile.write(msg+'\n')

        # Download files
        return_code = os.system(str('wget -N -q -r --directory-prefix='+data_dir+\
                               ' --no-parent '+ftp_config['update']))

        # Check for download error
        if return_code != 0:
            msg = str('Error: Update file download. Return code:'+str(return_code)+\
                  '\n Link: ' + ftp_config['update']+'\n')
            logfile.write(msg)
            exit(return_code)

        #  Report progress
        msg = str('Finished downloading PubMed update files.'+\
              str(round(time.time() - t1))+' seconds')
        logfile.write(msg+'\n')
        print(msg)
       

    
    
    
'''
MD5-checksum 
'''
def check_all_md5_in_dir(data_dir, logfile, mac, linux):
    '''
    FUNCTION:
    - Verify the downloaded PubMed article files via md5sum, verifies MD5 hashes.
    
    PARAMS:
    - data_dir (str): The main directory which will contain folders for both
      the baseline and update files.
    - logfile (file): The place to write the progress of downloading
    - mac (bool): Indicate if this is a MacOS
    - linux (bool): Indicate if this is a LinuxOS
    '''
    
    # Current year, last two suffix
    year = '22'
    
    # Check if the md5sum exists
    if os.system("which md5sum 1>/dev/null") != 0:
        print("md5sum not found")           
        # Continue executing
        return
    
    # Print progress starting checking md5
    msg = str('==== Start checking md5 in '+data_dir+' ====')
    logfile.write(msg+'\n')
    print(msg)
    
    # If the main 'data' directory exists
    if os.path.isdir(data_dir):
        file_count = 0

        # Iterate through each downloaded PubMed file
        for file in os.listdir(data_dir):
            
            # Check if the downloaded file is a zipped file containing PubMed articles
            if re.search('^pubmed'+year+'n\d\d\d\d.xml.gz$', file):
                
                # If this .xml.gz file and its .md5 file pass the md5 check, proceed
                check_md5(os.path.join(data_dir, file), mac, linux)
                
                # Print progress checking md5
                file_count += 1
                if file_count % 100 == 0:
                    print(str(file_count) +' files checked')
        
        # Print final progress of md5 checking
        msg = str('==== All md5 checks succeeded ('+str(file_count)+' files) ====\n')
        logfile.write(msg)
        print(msg)
        
        
    # If the main 'data' directory doesn't exist, print that. 
    else:
        print("Directory not found: %s (for md5 check)" % data_dir)
        logfile.write("==== Succeessful md5 checks (%d files) ====" % file_count+'\n')
        

def check_md5(file, mac, linux):
    '''
    FUNCTON:
    - Verify one file via md5 
    
    PARAMS:
    - file (file): Zipped file containing PubMed articles
    - mac (bool): Indicate if this is a MacOS
    - linux (bool): Indicate if this is a LinuxOS
    '''
    if os.path.isfile(file) and os.path.isfile(file + ".md5"):
        
        # Get md5 output from the zipped .xml
        if linux == True:
            stdout = subprocess.check_output('md5sum '+file, shell=True).decode('utf-8')
        elif mac == True:
            stdout = subprocess.check_output('md5 '+file, shell=True).decode('utf-8')
        md5_calculated = re.search('[0-9a-f]{32}', stdout).group(0)

        # Get md5 output from the corresponding .md5 file
        md5 = re.search('[0-9a-f]{32}', open(file + '.md5', 'r').readline()).group(0)

        # Check if the md5 are equal
        if md5 != md5_calculated:
            raise Exception('Error: md5 check failed for '+ file)
            exit(1)
            
            
            
            
            
            
'''
Extract data from downloaded files
'''
def extract(file, logfile):
    '''
    FUNCTION:
    - Extract data from the zipped downloaded PubMed data.
    
    PARAMS:
    - file (file): .gz file, zipped PubMed publications
    - logfile (file): The place to write the progress of downloading
    '''
    return_code = os.system('gunzip -fq '+ file) 
    if return_code != 0:
        msg = str('Error: Extraction. Return code for '+file+' is '+str(return_code))
        logfile.write(msg+'\n')
        print(msg)
        exit(return_code)
    return return_code


def extract_all_gz_in_dir(data_dir, logfile):
    '''
    FUNCTION:
    - Extract all .gz zipped files. The PubMed publications in in them. 
    
    PARAMS:
    - data_dir (str): The main directory which will contain folders for both
      the baseline and update files.
    - logfile (file): The place to write the progress of downloading
    '''
    
    # Check if the main 'data' directory exists
    if os.path.isdir(data_dir):
        count = 0
        print('==== Start extracting in '+data_dir+' ====')
        t1 = time.time()
        
        # For each downloaded file under the data directory
        for file in os.listdir(data_dir):
            
            # Extract it if it's a .gz file
            if re.search('.*\.gz$', file):
                extract(os.path.join(data_dir, file), logfile)
                
                # Print progress every 50 files extracted
                count += 1
                if count % 50 == 0:
                    msg = str(str(count)+' files extracted,'+str(time.time()-t1)+' seconds')
                    logfile.write(msg+'\n')
                    print(msg)
                
        # Print final progress of extracted files
        print('==== All files extracted ('+str(count)+' files). '\
              +'Total time: '+str(time.time() - t1)+' ====')

    else:
        print('Directory not found during extraction: '+data_dir) 