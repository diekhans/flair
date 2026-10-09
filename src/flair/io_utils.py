"""Filesystem / I/O helpers shared across FLAIR pipelines."""

import os
import shutil


def make_temp_dir(out_prefix):
    """The intermediate directory of a run, kept if it is already there, so a stage
    can reuse what an earlier run left.  A stage whose results depend on what the
    directory holds wants make_clean_temp_dir instead."""
    # FIXME: use TMPDIR unless directory explicitly specified
    temp_dir = out_prefix + ".intermediate"
    try:
        os.makedirs(temp_dir, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Creation of the directory `{temp_dir}' failed") from exc
    return temp_dir + '/'


def make_clean_temp_dir(out_prefix):
    """The intermediate directory of a run, emptied first.  For a stage that reads
    back every file it finds there: a file left by an earlier run, which a crash or
    --keep_intermediate can leave behind, would otherwise be counted as this run's."""
    temp_dir = out_prefix + ".intermediate"
    if os.path.exists(temp_dir):
        shutil.rmtree(temp_dir)
    return make_temp_dir(out_prefix)
