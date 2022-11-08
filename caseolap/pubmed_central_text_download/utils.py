import requests
from bs4 import BeautifulSoup
import os
from collections import namedtuple
import logging
import datetime

Datum = namedtuple('Datum',
                   [
                       'tar_path', 'csv_path', 'data_path', 'log_path',
                       'tmp_path', 'tar_url', 'csv_url', 'tar_name', 'csv_name'
                   ])

Urls = namedtuple('Urls', ['tar_url', 'csv_url'])
FileNames = namedtuple('FileNames',['tar_name', 'csv_name'])
Paths = namedtuple('Paths',
                   [
                       'tar_path', 'csv_path', 'data_path',
                       'log_path', 'tmp_path'
                   ])





log_levels = {'info': logging.INFO, 'debug': logging.DEBUG,
              'error': logging.ERROR, 'warning': logging.WARNING}

def parse_nih_oa_bulk_page(page: str):
    """
    Parses the nih.gov page and returns a list of files to download
    """
    soup = BeautifulSoup(page, 'html.parser')
    all_files = [link.get('href') for link in soup.find_all('a')]
    if not all_files:
        raise Exception("Did not find any links in page")
    return all_files

def make_full_text_processing_info(tar_urls, csv_urls, tmp_path, data_path, log_path):
    com_types = ['oa_comm_txt', 'oa_other_txt','oa_noncomm_txt']
    data = []
    for tar_url, csv_url in zip(tar_urls, csv_urls):
        for com_type in com_types:
            if com_type in tar_url:
                tar_name = tar_url[tar_url.find(com_type):]
                tar_path = tmp_path+tar_name
                csv_name = csv_url[csv_url.find(com_type):]
                csv_path = tmp_path+csv_name
        urls = Urls(csv_url, tar_url)
        paths = Paths(tar_path, csv_path, data_path, log_path, tmp_path)
        file_names = FileNames(tar_name, csv_name)
        data.append((urls, paths, file_names))
    return data

def set_up(tmp_path, data_path, log_path):
    os.makedirs(tmp_path, exist_ok=True)
    os.makedirs(data_path, exist_ok=True)
    os.makedirs(log_path)


def get_files_from_nih_page(oa_bulk_url):
    nih_oa_bulk_page = requests.get(oa_bulk_url)
    all_files = parse_nih_oa_bulk_page(nih_oa_bulk_page.text) 
    tar_files = [f for f in all_files if '.tar' in f]
    csv_files = [f for f in all_files if '.csv' in f]
    tar_urls = [oa_bulk_url + '/' + tar_url for tar_url in tar_files]
    csv_urls = [oa_bulk_url + '/' + csv_url for csv_url in csv_files]
    return tar_urls, csv_urls


def make_file_logger(name, log_level, log_path):
    logger = logging.getLogger(name)
    logger.setLevel(log_levels[log_level])
    fh = logging.FileHandler(log_path)
    logger.addHandler(fh)
    return logger
    
