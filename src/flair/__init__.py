import os
from collections import namedtuple
from flair.pycbio import NoStackError

VERSION = "3.0.0b1"

# Fixed minimum and maximums size for an intron, set to very
# conservative values.  This could be configurable at some point.
MIN_INTRON_SIZE = 20
MAX_INTRON_SIZE = 2000000

class FlairError(Exception):
    """General error condition in FLAIR"""
    pass

class FlairInputDataError(FlairError, NoStackError):
    """Error in FLAIR input data"""
    pass

def set_unix_path():
    "add programs in package to PATH."

    # FIXME: 2025-04-02 markd
    # this should be replaced with converting the exec use an explicit
    # path to make code more obvious.
    os.environ["PATH"] = os.path.dirname(os.path.realpath(__file__)) + ':' + os.environ["PATH"]


class Range(namedtuple("Range", ("start", "end"))):
    "container for a zero-based, half-open genome range"

    def __len__(self):
        return self.end - self.start
