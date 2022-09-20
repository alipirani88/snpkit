__author__ = 'alipirani'

import logging
from datetime import datetime

def keep_logging(pmessage, lmessage, logger, mode):
    #print(pmessage)
    if mode == 'warning':
        logger.warning(lmessage)
    elif mode == 'info':
        logger.info(lmessage)
    elif mode == 'exception':
        logger.exception(lmessage)
    elif mode == 'debug':
        logger.debug(lmessage)
    else:
        logger.error(lmessage)


def generate_logger(output_folder, analysis_name, log_unique_time):
    # Create a new logger; create a file handler; create a logging format; add the handlers to the logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    #handler = logging.FileHandler('{}/{}_{}_.log.txt'.format(args.output_folder, args.analysis_name, datetime.now().strftime('%Y_%m_%d_%H_%M_%S')))
    handler = logging.FileHandler('{}/{}_{}.log.txt'.format(output_folder, log_unique_time, analysis_name))
    handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger