# python reaxXtract/logger.py
import logging
import os
from typing import Optional



class DuplicateFilter(logging.Filter):
    def __init__(self):
        self.msgs = set()

    def filter(self, record):
        rv = record.msg not in self.msgs
        self.msgs.add(record.msg)
        return rv

def configure_log(level: Optional[str] = None, force: bool = False) -> None:
    """
    Configure the package logger.

    - level: optional string like "DEBUG", "INFO", ...; if None, used from env REAXXTRACT_LOG_LEVEL or default "INFO".
    - force: if True remove existing handlers before adding new one (useful in notebooks).

    Numeric log levels:
    - CRITICAL = 50
    - ERROR = 40
    - WARNING = 30
    - INFO = 20
    - DEBUG = 10
    - NOTSET = 0
    """
    if level is None:
        level = os.getenv("REAXXTRACT_LOG_LEVEL", "INFO")

    numeric_level = getattr(logging, level.upper(), logging.INFO)

    if force:
        # remove existing handlers to allow reconfiguration (notebooks)
        for h in list(log.handlers):
            log.removeHandler(h)

    if not log.handlers:
        handler = logging.StreamHandler()
        if level.upper() == "INFO":
            formatter = logging.Formatter("[%(levelname)s] %(message)s")
        else:   
            formatter = logging.Formatter("[%(levelname)s] %(module)s.%(funcName)s:%(lineno)d - %(message)s")
        handler.setFormatter(formatter)
        log.addHandler(handler)

    log.setLevel(numeric_level)

log = logging.getLogger("reaxXtract")
log.addFilter(DuplicateFilter())
# configure default on import (non-force) � can be overridden by calling configure(...) from entrypoint
configure_log()