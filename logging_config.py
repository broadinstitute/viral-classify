import logging
import sys
from datetime import datetime

def setup_logging():
    timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    log_filename = f"debug_blast_log_{timestamp}.log"  

    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(filename)s - %(funcName)s [%(lineno)d] - %(message)s')

    handlers = [
        logging.FileHandler(log_filename, 'a'),
        logging.StreamHandler(sys.stdout)  
    ]
    
    for handler in handlers:
        handler.setFormatter(formatter)
    
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    for handler in handlers:
        logger.addHandler(handler)

    return logger