import os
from multiprocessing import Pool
import logging
import datetime
from pprint import pprint
from utils import make_full_text_processing_info
from utils import get_files_from_nih_page, make_file_logger, set_up
from extract import extract_full_text_files
from download import download_file
import argparse
import time    
    
def download_and_extract_full_text(urls, paths, file_names):
    tar_path, csv_path, data_path, log_path, tmp_path = paths
    csv_url, tar_url = urls
    tar_name, csv_name = file_names
    main_logger = logging.getLogger('main')
    download_file(tar_url, tar_name, tmp_path, log_path)
    download_file(csv_url, csv_name, tmp_path, log_path)
    extract_full_text_files(tar_path, csv_path, data_path, log_path, tar_name)
    os.remove(tmp_path+tar_name)
    main_logger.info(f"[INFO] Deleted file {tar_name}")
    os.remove(tmp_path+csv_name)
    main_logger.info(f"[INFO] Deleted file {csv_name}")
    
oa_bulk_urls = [
    'https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_noncomm/txt/',
    'https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_other/txt/',
    'https://ftp.ncbi.nlm.nih.gov/pub/pmc/oa_bulk/oa_comm/txt/'
]

if __name__ == "__main__":
    start_time = time.time()
    parser = argparse.ArgumentParser()
    parser.add_argument('--log_path', help='path where logs will be stored')
    parser.add_argument('--pmcids', help='path to json with list of pmcids')
    parser.add_argument('--data_path', help='path to where data will be stored')
    parser.add_argument('--tmp_path', help='path where nih tar and csv will be stored temporarily')
    args = parser.parse_args()
    date_string = str(datetime.datetime.now()).replace(' ', "_")
    date_string = date_string.replace(':','_')
    date_string = date_string[:date_string.find('.')]
    log_path = args.log_path if args.log_path else f'logs/{date_string}/'
    pmcids = args.pmcids if args.pmcids else 'pmcids.json'
    data_path = args.data_path if args.data_path else 'data/'
    tmp_path = args.tmp_path if args.tmp_path else 'tmp/'
    set_up(tmp_path, data_path, log_path)
    main_logger = make_file_logger('main', 'info', log_path + "main" + '.log')
    tar_urls = []
    csv_urls = []
    pool = Pool(16)
    for oa_bulk_url in oa_bulk_urls:
        _tar_urls, _csv_urls = get_files_from_nih_page(oa_bulk_url)
        tar_urls.extend(_tar_urls)
        csv_urls.extend(_csv_urls)
    assert len(tar_urls) == len(csv_urls)
    data = make_full_text_processing_info(tar_urls, csv_urls, tmp_path, data_path, log_path)
    pool.starmap(download_and_extract_full_text, data)
    end_time = time.time()
    print(end_time-start_time)
