"""
Utilities for constructing the base table.
"""

import dtablib
import gzip
import os
import pandas as pd
import pickle


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

    return os.path.join(config['svpop_dir'], dtablib.svpop.VAR_PATTERN.format(**fmt_dict))

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

    return os.path.join(config['svpop_dir'], dtablib.svpop.FA_PATTERN.format(**fmt_dict))

def merge_base_table(input, output, wildcards, config):
    """
    Merge base table and write output files.

    :param input: Rule input.
    :param output: Rule output.
    :param wildcards: Rule wildcards.
    :param config: Pipeline config.
    """
    table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

    ### ID list ###
    id_set = set()

    with gzip.open(output.id_list, 'wt') as out_file:
        for file_name in input.id_list:
            if os.stat(file_name).st_size > 0:
                with gzip.open(file_name, 'rt') as in_file:
                    for line in in_file:
                        out_file.write(line)
                        id_set.add(line.strip())

    ### Sample set PKL ###
    obj_list = list()

    for file_name in input.sample_pkl:
        if os.stat(file_name).st_size > 0:
            with open(file_name, 'rb') as in_file:
                obj_list.append(pickle.load(in_file))

    if len(obj_list) == 0:
        raise RuntimeError('Found 0 data files to merge (all empty): sample_set')

    sample_set = pd.concat(obj_list, axis=0)

    if set(sample_set.index) != id_set:
        raise RuntimeError('ID set mismatch: sample_set')

    sample_set.to_pickle(output.sample_pkl)

    del sample_set

    ### Merge Map ###
    obj_list = list()

    for file_name in input.merge_map:
        if os.stat(file_name).st_size > 0:
            obj_list.append(pd.read_csv(file_name, sep='\t', low_memory=False))

    if len(obj_list) == 0:
        raise RuntimeError('Found 0 data files to merge (all empty): merge_map')

    df_sample_map = pd.concat(obj_list, axis=0)

    if set(df_sample_map['ID']) != id_set:
        raise RuntimeError('ID set mismatch: sample_map')

    df_sample_map = df_sample_map[['ID'] + [val for val in df_sample_map.columns if val != 'ID']]

    df_sample_map.to_csv(output.merge_map, sep='\t', index=False, compression='gzip')

    ### BED Base ###
    obj_list = list()

    for file_name in input.bed_base:
        if os.stat(file_name).st_size > 0:
            obj_list.append(pd.read_csv(file_name, sep='\t', low_memory=False))

    if len(obj_list) == 0:
        raise RuntimeError('Found 0 data files to merge (all empty): bed_base')

    df = pd.concat(obj_list, axis=0).sort_values(['#CHROM', 'POS'])

    if set(df['ID']) != id_set:
        raise RuntimeError('ID set mismatch: sample_map')

    # Check match against id_table
    if 'id_table' in table_def:

        # Read set of IDs to retain
        id_set = set(pd.read_csv(
            table_def['id_table'].format(**wildcards),
            sep='\t', usecols=('ID', ), squeeze=True
        ))

        # If retained ID set contains IDs not in the unfiltered table, error
        id_set_missing = id_set - set(df['ID'])

        if id_set_missing:
            raise RuntimeError('Filter ID set contains {} IDs not in the original table: {}{})'.format(
                len(id_set_missing),
                ','.join(sorted(id_set_missing)[:3]),
                '...' if len(id_set_missing) > 3 else ''
            ))

        del id_set

    df.to_csv(output.bed_base, sep='\t', index=False, compression='gzip')