import pandas as pd
from utils import make_file_logger
import logging
import json
import tarfile

def get_full_text_data(tar_path:str, full_text_names:str, logger):
    pmcids_json = json.load(open('pmcs.json','r'))
    pmcids = set(pmcids_json)
    main_logger = logging.getLogger('main')
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
                logger.error(f"[ERROR] Failed to extract {member} from {tar_path}")
                logger.error(f"[ERROR] Error message: {error}")
    

def get_full_text_file_names(csv_path):
    df = pd.read_csv(csv_path, usecols=['Article File'], chunksize=100)
    for data in df:
        for file_name in data['Article File']:
            yield file_name

def extract_full_text_files(tar_path, csv_path, data_path, log_path, tar_name):
    logger = make_file_logger("extract."+tar_path, 'info', log_path + f'extract_{tar_name}' + '.log')
    main_logger = logging.getLogger('main')
    main_logger.info(f'[INFO] Started extracting full text from {tar_name}')
    full_text_file_names = get_full_text_file_names(csv_path)
    full_text_data = get_full_text_data(tar_path, full_text_file_names, logger)
    files_processed = 0
    total_files = 0
    first_time_seeing_error = True
    for full_text_file_name, full_text in full_text_data:
        try:
            with open(data_path + full_text_file_name, 'wb') as f:
                f.write(full_text)
                files_processed += 1
            if files_processed % 1000 == 0:
                msg = f"Processed {files_processed} files"
                logger.info(msg)
        except Exception as error:
            if first_time_seeing_error:
                main_logger.warn(f"[WARNING] There are issues extracting file(s) from {tar_path}")
                first_time_seeing_error = False
            logger.error(f'[ERROR] Could not extract file {full_text_file_name} from {tar_name}')
            logger.error(f'[ERROR] Error message: {error}')
        total_files += 1
    logger.info(f"Finished extracting full text from {tar_path}")
    logger.info(f"Extracted {files_processed} pmc full text articles")
    main_logger.info(f"[INFO] Finished extracting full text from {tar_path}")
    main_logger.info(f"[INFO] Extracted {files_processed} pmc full text articlesout of {total_files}")
