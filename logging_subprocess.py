# Source: https://gist.github.com/hangtwenty/6390750
import sys
if sys.version_info < (3, 2):
    import subprocess32 as sp
else:
    import subprocess as sp
import select
import logging
from logging import DEBUG, ERROR

def call(popenargs, logger, stderr_log_level=logging.ERROR, stdout_log_level=logging.DEBUG, **kwargs):
    """
    Variant of subprocess.call that accepts a logger instead of stdout/stderr,
    and logs stdout messages via logger.debug and stderr messages via
    logger.error.
    """
    try:
        child = sp.Popen(popenargs, stdout=sp.PIPE, stderr=sp.PIPE, shell=True, **kwargs)
        log_level = {child.stdout: stdout_log_level, child.stderr: stderr_log_level}
        def check_io():
            ready_to_read = select.select([child.stdout, child.stderr], [], [], 1000)[0]
            for io in ready_to_read:
                line = io.readline()
                if line:
                    logger.log(log_level[io], line[:-1])
        # keep checking stdout/stderr until the child exits
        while child.poll() is None:
            check_io()
        check_io()  # check again to catch anything after the process exits
    except Exception as e:
        raise(e)
    return child.wait()
