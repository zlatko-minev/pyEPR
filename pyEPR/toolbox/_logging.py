import logging

from .. import config

def set_up_logger(logger):
    # add custom stream handler
    logger.c_handler = logging.StreamHandler()

    # Jupyter notebooks already has a stream handler on the default log,
    # Do not propagate upstream to the root logger.
    # https://stackoverflow.com/questions/31403679/python-logging-module-duplicated-console-output-ipython-notebook-qtconsole
    logger.propagate = False

    logger.c_format = logging.Formatter(config.log.format, config.log.datefmt)
    logger.c_handler.setFormatter(logger.c_format)
    logger.addHandler(logger.c_handler)
    logger.setLevel(getattr(logging, config.log.level))


