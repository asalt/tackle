# logger.py
import logging


def get_logger(name):
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(ch)
    try:
        fh = logging.FileHandler("tackle.log")
        fh.setFormatter(formatter)
        logger.addHandler(fh)
    except PermissionError:
        pass

    # fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    return logger
