import requests
from bs4 import BeautifulSoup
from collections import namedtuple

Datum = namedtuple('Datum',
                   [
                       'tar_path', 'csv_path', 'data_path', 'log_path',
                       'tmp_path', 'tar_url', 'csv_url', 'tar_name', 'csv_name'
                   ])

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

def make_full_text_processing_info(tar_urls, csv_urls, tmp_path, data_path, log_path):
    data = []
    for tar_url, csv_url in zip(tar_urls, csv_urls):
        tar_name = tar_url[tar_url.find('oa_comm_txt'):]
        tar_path = tmp_path+tar_name
        csv_name = csv_url[csv_url.find('oa_comm_txt'):]
        csv_path = tmp_path+csv_name
        datum = Datum(tar_path, csv_path, data_path, log_path,
                      tmp_path, tar_url, csv_url, tar_name, csv_name)
        data.append(datum)
    return data
