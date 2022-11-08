import logging
import requests
from utils import make_file_logger

def download_file(url: str, file_name: str, tmp_path: str, log_path:str) -> None:
    logger = make_file_logger(file_name, 'info', log_path + file_name + '.log')
    main_logger = logging.getLogger('main')
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
        main_logger.info(f"[INFO] Finished downloading {file_name}")
    except Exception as error:
        logger.error(f"[ERROR] Failed to download file {file_name}")
        logger.error(f"[ERROR] Error message: {error}")
