import os
from collections import namedtuple
from types import NoneType
from flair.pycbio import NoStackError

VERSION = "3.0.0b1"

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


class PosRange(namedtuple("PosRange",
                          ("start", "end"))):
    """0-base, 1/2 open range of sequence positions"""

    def __new__(cls, start, end):
        assert start <= end  # allows zero length
        return super(PosRange, cls).__new__(cls, start, end)

class SeqRange(namedtuple("SeqRange",
                          ("name", "start", "end", "strand"))):
    """Range withing a sequence, including optional strand. Duck-type inheritance
    from PosRange"""

    def __new__(cls, name, start, end, strand=None):
        assert start <= end  # allows zero length
        assert isinstance(strand, (str, NoneType)) and ((strand is None) or (strand in ('+', '-', '.', None)))
        return super(SeqRange, cls).__new__(cls, name, start, end, strand)
