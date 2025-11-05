# MIT License
#
# Copyright (c) 2025 Isacks-Co contributors
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import logging


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

def loggerSetup(debug):
    log = logging.getLogger()
    log.setLevel(logging.INFO if not debug else logging.DEBUG)
    log.propagate = False        

    if not log.handlers:           
        h = logging.StreamHandler()
        fmt = "%(asctime)s | %(levelname)s: %(message)s"
        datefmt = "%H:%M:%S.%ms"  
        h.setFormatter(CustomFormatter(fmt=fmt, datefmt=datefmt))        
        log.addHandler(h)

