#! /usr/bin/env python3
import time
from datetime import timedelta
import argparse
import logging
from flair import VERSION, set_unix_path
from flair.pycbio.sys import cli, loggingOps


###
# definitions of modules
###

ALIGN_MODULE = 'align'
CORRECT_MODULE = 'correct'
TRANSCRIPTOME_MODULE = 'transcriptome'
COLLAPSE_MODULE = 'collapse'
QUANTIFY_MODULE = 'quantify'
COMBINE_MODULE = 'combine'
VARIANTS_MODULE = 'variants'
FUSION_MODULE = 'fusion'
DIFFEXP_MODULE = 'diffexp'
DIFFSPLICE_MODULE = 'diffsplice'


MODULES = (
    ALIGN_MODULE,
    CORRECT_MODULE,
    TRANSCRIPTOME_MODULE,
    COLLAPSE_MODULE,
    QUANTIFY_MODULE,
    COMBINE_MODULE,
    VARIANTS_MODULE,
    FUSION_MODULE,
    DIFFEXP_MODULE,
    DIFFSPLICE_MODULE,
)

MODULES_DESC = {
    ALIGN_MODULE: "Align reads to the genome",
    CORRECT_MODULE: "Correct misaligned splice sites using genome annotations and/or short-read splice junctions",
    TRANSCRIPTOME_MODULE: "Generates a transcriptome of isoforms directly from a long-read BAM",
    COLLAPSE_MODULE: "Define high-confidence isoforms from corrected reads",
    QUANTIFY_MODULE: "Quantify the expression level of the predicted isoforms",
    COMBINE_MODULE: "Combine FLAIR transcriptomes or annotation transcriptomes",
    VARIANTS_MODULE: "Integrate variants from other sources with the predicted isoforms",
    FUSION_MODULE: "Identify gene fusions and generates a fusion transcriptome",
    DIFFEXP_MODULE: "Differential expression and differential usage analysis",
    DIFFSPLICE_MODULE: "Call alternative splicing events from predicted isoforms",
}

###
# Command line parsing
#
# This is bit tricky to achieve:
#   - subcommands for each module
#   - options to the overall to enable logging, print version, and maybe more
#   - provide good help messages for both top-level and subcommands
#   - delay loading modules until needed to allow for optional dependencies
#
# ArgumentParser.parse_known_args() is not well behaved with help.  This avoids
# the issue by pre-parsing to get the module name, the do the real parsing
# in the module.
#
###

def create_top_args_parent_parser():
    """create a parent parser exists and is used only to add common, top level options"""
    top_args_parent = argparse.ArgumentParser(add_help=False)
    top_args_parent.add_argument('--version', action='version', version='FLAIR ' + VERSION,
                                 help="print FLAIR version")
    loggingOps.addCmdOptions(top_args_parent)
    return top_args_parent

def create_module_subparsers(parser, top_args_parent, add_help=True):
    subparsers = parser.add_subparsers(dest="module", required=True,
                                       help="the FLAIR module to run")
    for module in MODULES_DESC:
        subparsers.add_parser(module, help=MODULES_DESC[module], add_help=add_help,
                              parents=[top_args_parent])
    return subparsers

def preparse_module_name(top_args_parent):
    "module name or None if not valid"
    # no help here to prevent printing help and exiting
    pre_parser = argparse.ArgumentParser(add_help=False, exit_on_error=False,
                                         parents=[top_args_parent])
    create_module_subparsers(pre_parser, top_args_parent, add_help=False)
    try:
        known_args, _ = pre_parser.parse_known_args()
        return known_args.module
    except argparse.ArgumentError:
        return None

def setup_parser():
    """Setup argument parser and pre-parse to get module name.
    call of parse_args() is passed to the module
    """
    top_args_parent = argparse.ArgumentParser(add_help=False)
    module = preparse_module_name(top_args_parent)

    # main parser
    desc = 'Run FLAIR modules for the FLAIR analysis pipeline'
    parser = argparse.ArgumentParser(description=desc,
                                     parents=[top_args_parent])
    print("@", "parser", loggingOps.haveCmdOptions(parser))
    subparsers = create_module_subparsers(parser, top_args_parent)
    if module is None:
        # let parser display the error message
        parser.parse_args()
        raise RuntimeError("BUG: setup_parser")

    return parser, subparsers.choices[module], module

def flair_module_run(parser, subparser, module):  # noqa: C901
    start_time = time.time()

    # Delay importing modules until needed to allow for optional dependencies.
    if module == ALIGN_MODULE:
        from flair import flair_align
        flair_align.align_subcommand(parser, subparser)
    elif module == CORRECT_MODULE:
        from flair.flair_correct import correct
        correct()
    elif module == TRANSCRIPTOME_MODULE:
        from flair.flair_transcriptome import collapsefrombam
        collapsefrombam()
    elif module == COLLAPSE_MODULE:
        from flair.flair_collapse import collapse
        collapse()
    elif module == QUANTIFY_MODULE:
        from flair.flair_quantify import quantify
        quantify()
    elif module == COMBINE_MODULE:
        from flair import flair_combine
        flair_combine.combine()
    elif module == VARIANTS_MODULE:
        from flair.flair_variants import getvariants
        getvariants()
    elif module == FUSION_MODULE:
        from flair.flair_fusion import detectfusions
        detectfusions()
    elif module == DIFFEXP_MODULE:
        from flair import flair_diffExp
        flair_diffExp.diffExp()
    elif module == DIFFSPLICE_MODULE:
        from flair import flair_diffSplice
        flair_diffSplice.diffSplice()

    elapsed = time.time() - start_time
    logging.info(f"Flair {module} took " + str(timedelta(seconds=round(elapsed))))

def main():
    set_unix_path()
    parser, subparser, module = setup_parser()
    with cli.ErrorHandler():
        flair_module_run(parser, subparser, module)


if __name__ == '__main__':
    main()
