import os
from contextlib import contextmanager

@contextmanager
def suppress_cpp_output():
    with open(os.devnull, 'w') as devnull:
        # Duplicate original fds so we can restore later
        old_stdout = os.dup(1)
        old_stderr = os.dup(2)

        try:
            # Redirect stdout and stderr (C++ uses these)
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            # Restore original fds
            os.dup2(old_stdout, 1)
            os.dup2(old_stderr, 2)
            os.close(old_stdout)
            os.close(old_stderr)
