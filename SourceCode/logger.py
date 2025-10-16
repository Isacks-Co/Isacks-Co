import logging
import time

class CustomFormatter(logging.Formatter):
    COLORS =  {
    logging.DEBUG: "\033[36m",   # Grå/vit
    logging.INFO: "\033[32m",    # Grön
    logging.WARNING: "\033[33m", # Gul
    logging.ERROR: "\033[31m",   # Röd
    logging.CRITICAL: "\033[41m" # Röd bakgrund
    }

    RESET = "\033[0m"

    def format(self, record):
        color = self.COLORS.get(record.levelno, "")
        orig_levelname = record.levelname
        record.levelname = f"{color}{orig_levelname}{self.RESET}" if color else orig_levelname
        message = super().format(record)
        return f"{color}{message}{self.RESET}"

def loggerSetup(debug=True):
    log = logging.getLogger()
    log.setLevel(logging.INFO if not debug else logging.DEBUG)
    log.propagate = False        

    if not log.handlers:           
        h = logging.StreamHandler()
        fmt = "%(asctime)s | %(levelname)s: %(message)s"
        datefmt = "%H:%M:%S.%ms"  
        h.setFormatter(CustomFormatter(fmt=fmt, datefmt=datefmt))        
        log.addHandler(h)

