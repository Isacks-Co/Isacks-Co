import logging

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
        color = self.COLORS.get(record.levelno, self.RESET)
        message = super().format(record)
        return f"{color}{message}{self.RESET}"



def logger_setup():
    log = logging.getLogger("MD")
    log.setLevel(logging.INFO)        # byt till DEBUG vid felsökning
    log.propagate = False             # så vi inte bubblar upp till root och får dubletter

    if not log.handlers:              # skydda mot dubbla handlers om filen körs flera gånger
        h = logging.StreamHandler()
        h.setFormatter(CustomFormatter("%(levelname)s: %(message)s"))
        log.addHandler(h)
        return log
