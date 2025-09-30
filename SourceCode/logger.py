import logging
import time

class CustomFormatter(logging.Formatter):
    COLORS =  {
    logging.DEBUG: "\033[37m",   # Grå/vit
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



def logger_setup():
    log = logging.getLogger("MD")
    log.setLevel(logging.INFO)        # byt till DEBUG vid felsökning
    log.propagate = False        

    if not log.handlers:           
        h = logging.StreamHandler()
        fmt = "%(asctime)s | %(levelname)s: %(message)s"
        datefmt = "%H:%M:%S.%ms"  
        h.setFormatter(CustomFormatter(fmt=fmt, datefmt=datefmt))        
        log.addHandler(h)

    return log



def log_timing(msg, t0 ):
    log = logging.getLogger("MD")
    dt = time.perf_counter() - t0
    log.info("%-35s | %8.3f s", msg, dt)