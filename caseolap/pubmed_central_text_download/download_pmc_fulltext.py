import requests
from bs4 import BeautifulSoup
import pandas as pd
import tarfile
import time
import os
from multiprocessing import Pool
import logging
import datetime
from pprint import pprint
from utils import make_full_text_processing_info
import json

def parse_nih_oa_bulk_page(page: str):
    """
    Parses the nih.gov page and returns a list of files to download
    """
    soup = BeautifulSoup(page, 'html.parser')
    all_files = []
    for link in soup.find_all('a'):
        all_files.append(link.get('href'))
    if not all_files:
        raise Exception("Did not find any links in page")
    return all_files

def download_file(url: str, file_name: str, tmp_path: str, log_path:str):
    logger = logging.getLogger(file_name)
    logger.setLevel(logging.INFO)
    main_logger = logging.getLogger('main')
    fh = logging.FileHandler(log_path + file_name + '.log')
    logger.addHandler(fh)
    try:
        response = requests.get(url, stream=True)
        number_of_bytes = int(response.headers["Content-length"])
        logger.info(f"[INFO] Number of bytes to download for {file_name}: {number_of_bytes} ")
        chunk_size = 512
        last_percentage = 0
        if response.status_code == 200:
            with open(tmp_path+file_name, 'wb') as f:
                for chunk_index, chunk in enumerate(response.iter_content(chunk_size=chunk_size)):
                    f.write(chunk)
                    percent = (((chunk_index+1) * chunk_size) / number_of_bytes) * 100
                    if percent - last_percentage >= 1:
                        msg = f"[INFO] At {percent:.2f}% for download of {file_name}"
                        last_percentage = percent
                        logger.info(msg)
        logger.info(f"[INFO] Finished downloading {number_of_bytes} bytes")
        main_logger.info(f"[INFO] Finished downloading {file_name}")
    except Exception as error:
        print(f"Failed to download {file_name}")
        print(f"Response status code: {response.status_code}")
        logger.error(f"[ERROR] Failed to download file {file_name}")
        logger.error(f"[ERROR] Error message: {error}")
    
def get_full_text_data(tar_path:str, full_text_names:str, logger):
    pmcids_json = json.load(open('pmcs.json','r'))
    pmcids = set(pmcids_json)
    main_logger = logging.getLogger('main')
    seen_issue = True
    with tarfile.open(tar_path) as tar:
        for member in full_text_names:
            try:
                file_name = member.split('/')[1]
                pmcid = file_name.replace("PMC",'').replace('.txt','')
                if pmcid in pmcids:
                    f = tar.extractfile(member)
                    content = f.read()
                    yield file_name, content
            except Exception as error:
                if seen_issue:
                    seen_issue = False
                    main_logger.error(f"[ERROR] {tar_path} had an extraction failure")
                    main_logger.error(f"[ERROR] Error message for {tar_path}: {error}")
                logger.error(f"[ERROR] Failed to extract {member} from {tar_path}")
                logger.error(f"[ERROR] Error message: {error}")
    

def get_full_text_file_names(csv_path):
    df = pd.read_csv(csv_path, usecols=['Article File'], chunksize=100)
    for data in df:
        for file_name in data['Article File']:
            yield file_name

def extract_full_text_files(tar_path, csv_path, data_path, log_path, tar_name):
    logger = logging.getLogger("extract."+tar_path)
    logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_path + f'extract_{tar_name}' + '.log')
    logger.addHandler(fh)
    main_logger = logging.getLogger('main')
    main_logger.info(f'Started extracting full text from {tar_name}')
    full_text_file_names = get_full_text_file_names(csv_path)
    full_text_data = get_full_text_data(tar_path, full_text_file_names, logger)
    files_processed = 0
    total_files = 0
    for full_text_file_name, full_text in full_text_data:
        try:
            with open(data_path + full_text_file_name, 'wb') as f:
                f.write(full_text)
                files_processed += 1
            if files_processed % 1000 == 0:
                msg = f"Processed {files_processed} files"
                logger.info(msg)
        except Exception as error:
            logger.error(f'[ERROR] Could not extract file {full_text_file_name} from {tar_name}')
            logger.error(f'[ERROR] Error message: {error}')
        total_files += 1
    logger.info(f"Finished extracting full text from {tar_path}")
    logger.info(f"Extracted {files_processed} pmc full text articles")
    main_logger.info(f"[INFO] Finished extracting full text from {tar_path}")
    main_logger.info(f"[INFO] Extracted {files_processed} pmc full text articles")
    
    
def download_and_extract_full_text(datum):
    tar_path, csv_path, data_path, log_path, tmp_path, tar_url, csv_url, tar_name, csv_name = datum
    main_logger = logging.getLogger('main')
    download_file(tar_url, tar_name, tmp_path, log_path)
    download_file(csv_url, csv_name, tmp_path, log_path)
    extract_full_text_files(tar_path, csv_path, data_path, log_path, tar_name)
    os.remove(tmp_path+tar_name)
    main_logger.info(f"[INFO] Deleted file {tar_name}")
    os.remove(tmp_path+csv_name)
    main_logger.info(f"[INFO] Deleted file {csv_name}")
    
    
oa_bulk_url = 'https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/'

if __name__ == "__main__":
    #Make these configurable
    current_date_string = str(datetime.datetime.now()).replace(' ', "_").replace(':','_')
    current_date_string = current_date_string[:current_date_string.find('.')]
    log_path = f'logs/{current_date_string}/'
    tmp_path = 'tmp/'
    data_path = 'data/'
    os.makedirs(tmp_path, exist_ok=True)
    os.makedirs(data_path, exist_ok=True)
    os.makedirs(log_path)
    main_logger = logging.getLogger('main')
    main_logger.setLevel(logging.INFO)
    fh = logging.FileHandler(log_path + "main" + '.log')
    main_logger.addHandler(fh)
    pool = Pool(8)
    nih_oa_bulk_page = requests.get(oa_bulk_url)
    all_files = parse_nih_oa_bulk_page(nih_oa_bulk_page.text)
    tar_files = [f for f in all_files if '.tar' in f]
    csv_files = [f for f in all_files if '.csv' in f]
    tar_urls = [oa_bulk_url + '/' + tar_url for tar_url in tar_files]
    csv_urls = [oa_bulk_url + '/' + csv_url for csv_url in csv_files]
    data = make_full_text_processing_info(tar_urls, csv_urls, tmp_path, data_path, log_path)
    pool.map(download_and_extract_full_text, data)
