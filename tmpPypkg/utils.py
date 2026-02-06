from typing import List
import sys
import logging
import inspect

def any_in(x: List[str], pool: set) -> bool:
    for i in x:
        if i in pool:
            return True
    return False

def reimport(module: str) -> None:
    import importlib
    importlib.reload(sys.modules[module])

def flatten_list(xss: List[List]) -> List:
    return [x for xs in xss for x in xs]

def unique_list(xs: List) -> List:
    return list(set(xs))

# about log
class StreamToLogger(object):
    """
    Fake file-like stream object that redirects writes to a logger instance.
    Ref:
    https://stackoverflow.com/questions/19425736/how-to-redirect-stdout-and-stderr-to-logger-in-python
    """

    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self.linebuf = ""

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass


def set_file_logger(
    fnm: str, fmode: str = "a", name: str = "sa2_pp",
    log_level: int = logging.DEBUG
) -> logging.Logger:
    logger = logging.getLogger(name)
    logger.setLevel(log_level)
    fh = logging.FileHandler(filename=fnm, mode=fmode)
    fm = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    fh.setFormatter(fm)
    logger.addHandler(fh)
    return logger


def handle_exception(logger, exc_type, exc_value, exc_traceback):
    import traceback

    if issubclass(exc_type, KeyboardInterrupt):
        sys.__excepthook__(exc_type, exc_value, exc_traceback)
        return
    logger.error(
        "".join(
            [
                "Uncaught exception: ",
                *traceback.format_exception(exc_type,
                                            exc_value, exc_traceback),
            ]
        )
    )

def chunk_list(lst, size: int = 1000):
    return [lst[i:i + size] for i in range(0, len(lst), size)]

def show(f) -> None:
    source = inspect.getsource(f)
    print(source)
