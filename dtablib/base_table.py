"""
Utilities for constructing the base table.
"""

import dtablib
import os


def get_input_bed(wildcards, config):
    """
    Get input BED.

    :param wildcards: Rule wildcards.

    :return: Path to SV-Pop BED file.
    """

    table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

    fmt_dict = {}
    fmt_dict.update(table_def)
    fmt_dict.update(wildcards)

    return os.path.join(config['svpop_run_dir'], dtablib.svpop.VAR_PATTERN.format(**fmt_dict))

def get_input_fa(wildcards, config):
    """
    Get input FASTA file.

    :param wildcards: Rule wildcards.

    :return: Path to SV-Pop FASTA file.
    """

    table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

    fmt_dict = {}
    fmt_dict.update(table_def)
    fmt_dict.update(wildcards)

    return os.path.join(config['svpop_run_dir'], dtablib.svpop.FA_PATTERN.format(**fmt_dict))
