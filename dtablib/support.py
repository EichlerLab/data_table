"""
Functions for external support (validation methods, other callers, etc).
"""

import collections
import os
import pandas as pd
import re
import dtablib


def get_conf(wildcards, config):

    table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

    # Get section
    conf_support = table_def.get(
        'support', dict()
    ).get(
        f'{wildcards.vartype}_{wildcards.svtype}', dict()
    ).get(
        wildcards.support_section, None
    )

    if conf_support is None:
        raise RuntimeError(f'No support section "{wildcards.support_section}" in config')

    conf_support = conf_support.copy()

    # Check section
    missing_index = {'path', 'type'} - set(conf_support.keys())

    if missing_index:
        raise RuntimeError(
            'Missing config elements for support section {}: {}'.format(
                wildcards.support_section,
                ', '.join(sorted(missing_index))
            )
        )

    if 'name' in conf_support:
        raise RuntimeError('Support section {} has reserved keyword "name"'.format(wildcards.support_section))

    conf_support['name'] = wildcards.support_section

    # Check path wildcards
    if '{sample}' not in conf_support['path'] and conf_support['type'] not in {'preformat', 'svpopinter-striphap'}:
        raise RuntimeError(
            'Support section {} "path" element is missing the "{{sample}}" wildcard'.format(wildcards.support_section)
        )

    if conf_support['type'] == 'svpopinter-striphap':
        missing_list = [
            val  for val in ('sample_l', 'sample_r')
                if val not in conf_support['path']
        ]

        if missing_list:
            raise RuntimeError(
                'Support section {} is missing sample wildcard(s): {}'.format(wildcards.support_section, ', '.join(missing_list))
            )

    if '{vartype}' not in conf_support['path']:
        raise RuntimeError(
            'Support section {} "path" element is missing the "{{vartype}}" wildcard'.format(wildcards.support_section)
        )

    if '{svtype}' not in conf_support['path']:
        raise RuntimeError(
            'Support section {} "path" element is missing the "{{svtype}}" wildcard'.format(wildcards.support_section)
        )

    # Make section wildcard substitutions (section has a "wildcard_list" element)
    if table_def['wildcards'] is not None and len(table_def['wildcards']) > 0:
        conf_support['path'] = conf_support['path'].format(**table_def['wildcards'])

    # Check allow-missing
    if 'allow-missing' in conf_support:
        conf_support['allow-missing'] = dtablib.util.get_bool(conf_support['allow-missing'])
    else:
        conf_support['allow-missing'] = False

    # Get column name
    if 'column-name' not in conf_support:
        conf_support['column-name'] = wildcards.support_section.upper()

    # Return
    return conf_support


def get_input_dict(wildcards, config):
    """
    Get input files for each sample. The dictionary returned is keyed by sample name. Values are lists of input files
    for that sample. If the base rule were not run (no base table generated yet), return `None`.

    :param wildcards: Rule wildcards.
    :param config: Pipeline config.

    :return: Dictionary of samples containing lists of input files per sample. Return `None` if base table was not
        generated, which is needed for the sample list.
    """

    conf_support = get_conf(wildcards, config)

    # Get samples from the ID map table header
    base_file_name = 'sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}.tsv.gz'.format(**wildcards)

    if not os.path.isfile(base_file_name):
        raise RuntimeError('Missing base file name: {}'.format(base_file_name))

    sample_list = list(pd.read_csv(base_file_name, sep='\t', nrows=0).squeeze().columns)[1:]

    # Get input files
    input_file_dict = collections.defaultdict(list)

    if conf_support['type'] == 'preformat':
        # If preformat, check for input file (not by sample or type)

        file_name = conf_support['path'].format(**wildcards)

        if not os.path.isfile(file_name):
            raise RuntimeError('Missing preformat table for section {}: {}'.format(conf_support['name'], file_name))

        input_file_dict['preformat'].append(file_name)

    elif conf_support['type'] == 'svpopinter-striphap':
        # Strip the haplotype off the second callset. Use option if samples on the first callset have "-h1" and "-h2",
        # but samples on the second callset does not. For example, intersecting a phased callset (left) with an
        # unphased callset (right).

        # Check sample names
        hap_re = re.compile(r'^(.*)-h[\d]+$')

        bad_list = [
            sample for sample in sample_list if re.search(hap_re, sample) is None
        ]

        if bad_list:
            n_bad = len(bad_list)
            bad_list = ', '.join(bad_list[:3])
            elipses = '...' if n_bad > 3 else ''

            raise RuntimeError(
                f'Sample name error in {wildcards.support_section}: Expected sample names to end with haplotype ("SAMPLE-h1", "SAMPLE-h2"): Found {n_bad} sample names without haplotypes at the end: {bad_list}{elipses}'
            )

        # Get input files per sample and type
        kw_dict = dict(wildcards)

        for sample in sample_list:
            kw_dict['sample_l'] = sample
            kw_dict['sample_r'] = re.search(hap_re, sample)[1]

            for svtype in dtablib.definitions.SVTYPE_EXPAND[wildcards.svtype]:
                kw_dict['svtype'] = svtype

                file_name = conf_support['path'].format(**kw_dict)

                if os.path.isfile(file_name):
                    input_file_dict[sample].append(file_name)
                else:
                    if not conf_support['allow-missing']:
                        raise RuntimeError('Missing support file for section {}: {}'.format(conf_support['name'], file_name))

    else:

        # Get input files per sample and type
        kw_dict = dict(wildcards)

        for sample in sample_list:
            kw_dict['sample'] = sample

            for svtype in dtablib.definitions.SVTYPE_EXPAND[wildcards.svtype]:
                kw_dict['svtype'] = svtype

                file_name = conf_support['path'].format(**kw_dict)

                if os.path.isfile(file_name):
                    input_file_dict[sample].append(file_name)
                else:
                    if not conf_support['allow-missing']:
                        raise RuntimeError('Missing support file for section {}: {}'.format(conf_support['name'], file_name))

    return input_file_dict


def get_input_file_list(wildcards, config):
    """
    Get supporting input files as a flat list for input into a Snakemake rule.

    :param wildcards: Rule wildcards.
    :param config: Pipeline config.

    :return: List of input support files.
    """

    input_file_dict = get_input_dict(wildcards, config).values()

    if input_file_dict is None:
        return []

    return [file for file_list in input_file_dict for file in file_list]


def get_support_svpopinter_lead(id_set, input_dict, support_col_name):
    """
    Get support using an SV-Pop intersect table. The table has two columns, "ID_A" and "ID_B", listing calls in one
    that intersected with the other. This function assumes the base table calls are in "ID_A" and the supporting
    calls are in "ID_B". Supported variants have an ID in both columns.

    :param id_set: Dictionary of variant ID sets (from the base table) keyed by the sample they came from.
    :param input_dict: Dictionary of support input files keyed by the sample name.
    :param support_col_name: Name of the column name to be added.

    :return: DataTable with an "ID" column (for the base IDs) and a column named `conf_support` with the variant
        ID supporting the base variant.
    """

    df_list = list()

    for sample in id_set.keys():
        sample_id_set = id_set[sample]

        # Skip if input file list is empty
        if not input_dict[sample]:
            continue

        # Get support and subset for matched (supported) variant calls
        # ID_A: Variants merged into the data table
        # ID_B: Variants from the supporting callset
        df_support = pd.concat(
            [pd.read_csv(in_file_name, sep='\t', usecols=('ID_A', 'ID_B')) for in_file_name in input_dict[sample]],
            axis=0
        )

        df_support = df_support.loc[
            (
                ~ pd.isnull(df_support['ID_A'])
            ) & (
                ~ pd.isnull(df_support['ID_B'])
            )
        ]

        df_support = df_support.loc[df_support['ID_A'].apply(lambda val: val in sample_id_set)]

        df_list.append(df_support)

    # Merge support tables
    if df_list:
        df_support = pd.concat(df_list, axis=0)

        df_support = df_support[['ID_A', 'ID_B']]
        df_support.columns = [['ID', support_col_name]]
    else:
        df_support = pd.DataFrame([], columns=['ID', support_col_name])

    # Return
    return df_support


def get_support_subseq_lead(id_set, input_dict, support_col_name):
    """
    Get support from a subseq validation table. Tables have variant IDs and a validation column.

    :param id_set: Dictionary of variant ID sets (from the base table) keyed by the sample they came from.
    :param input_dict: Dictionary of support input files keyed by the sample name.
    :param support_col_name: Column name to be added.

    :return: DataTable with an "ID" column (for the base IDs) and a column named `conf_support` with the variant
        ID supporting the base variant.
    """

    df_list = list()

    for sample in id_set.keys():
        sample_id_set = id_set[sample]

        # Skip if input file list is empty
        if not input_dict[sample]:
            continue

        # Get support and subset for matched (supported) variant calls
        # ID_A: Variants merged into the data table
        # ID_B: Variants from the supporting callset
        df_support = pd.concat(
            [pd.read_csv(in_file_name, sep='\t', usecols=('ID', 'VAL')) for in_file_name in input_dict[sample]],
            axis=0
        )

        df_support = df_support[['ID', 'VAL']]

        df_support = df_support.loc[df_support['ID'].apply(lambda val: val in sample_id_set)]

        df_support.columns = ['ID', support_col_name]

        df_list.append(df_support)

    # Merge support tables
    df_support = pd.concat(df_list, axis=0)

    # Return
    return df_support

def get_support_table_lead(id_set, input_dict, support_col_name, conf_support):
    """
    Get support using a collection of pre-formatted tables (per sample, per svtype.

    :param id_set: Dictionary of variant ID sets (from the base table) keyed by the sample they came from.
    :param input_dict: Dictionary of support input files keyed by the sample name.
    :param support_col_name: Name of the column name to be added.
    :param conf_support: Support config section.

    :return: DataTable with an "ID" column (for the base IDs) and a column named `conf_support` with the variant
        ID supporting the base variant.
    """

    if 'table_col' not in conf_support:
        raise RuntimeError(
            'Table config for support section {} is missing field "table_col" (name of table column to extract)'.format(
                conf_support['name']
            )
        )

    table_col = conf_support['table_col']

    df_list = list()

    for sample in id_set.keys():
        sample_id_set = id_set[sample]

        # Skip if input file list is empty
        if not input_dict[sample]:
            continue

        # Get support and subset for matched (supported) variant calls
        # ID_A: Variants merged into the data table
        # ID_B: Variants from the supporting callset
        df_support = pd.concat(
            [pd.read_csv(in_file_name, sep='\t') for in_file_name in input_dict[sample]],
            axis=0
        )

        missing_col = {'ID', table_col} - set(df_support.columns)

        if missing_col:
            raise RuntimeError('Missing columns in support tables for section {}: {}'.format(
                conf_support['name'],
                ', '.join(sorted(missing_col))
            ))

        df_support = df_support.loc[
            df_support['ID'].apply(lambda val: val in sample_id_set),
            ['ID', table_col]
        ]

        df_support.columns = ['ID', support_col_name]

        df_list.append(df_support)

    # Merge support tables
    if df_list:
        df_support = pd.concat(df_list, axis=0)
    else:
        df_support = pd.DataFrame([], columns=['ID', support_col_name])

    # Return
    return df_support