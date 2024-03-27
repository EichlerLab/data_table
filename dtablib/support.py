"""
Functions for external support (validation methods, other callers, etc).
"""

import collections
import os
import pandas as pd
import re
import dtablib


def get_conf(wildcards, config, support_section=None):
    """
    Get support section config.

    :param wildcards: Rule wildcards.
    :param config: Pipeline config.
    :param support_section: Support section or `None` to read section `wildcards.support_section`.

    :return: Dictionary of config.
    """

    if support_section is None:
        support_section = wildcards.support_section

    table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

    section_wildcards = dict()  # Wildcards to parse into the path

    # Get section
    conf_support = table_def.get(
        'support', dict()
    ).get(
        f'{wildcards.vartype}_{wildcards.svtype}', dict()
    ).get(
        support_section, None
    )

    if conf_support is None:
        raise RuntimeError(f'No support section "{support_section}" in config')

    if isinstance(conf_support, str):
        conf_support = [conf_support]

    if isinstance(conf_support, list):
        conf_support_name = conf_support[0]

        for wildcard_avp in conf_support[1:]:
            tok = wildcard_avp.split('=', 1)

            if len(tok) != 2:
                raise RuntimeError(f'Support section "{support_section}" ({wildcards.vartype}_{wildcards.svtype}) contains a wildcard definition with no "=" (expected key=value): "{wildcard_avp}"')

            tok[0] = tok[0].strip()
            tok[1] = tok[1].strip()

            if tok[0] == '':
                raise RuntimeError(f'Support section "{support_section}" ({wildcards.vartype}_{wildcards.svtype}) contains a wildcard definition with no key (empty before "="): "{wildcard_avp}"')

            section_wildcards[tok[0]] = tok[1]

        # Section is a name that refers to support_conf section
        conf_support = table_def.get(
            'support_sections', dict()
        ).get(
            conf_support_name, dict()
        )

        if conf_support is None:
            raise RuntimeError(f'Support section "{support_section}" ({wildcards.vartype}_{wildcards.svtype}) refers to a missisng pre-configured support section "{conf_support_name}" of the configuration file (expected to find in section "support_sections")')

    conf_support = conf_support.copy()

    # Check section
    missing_index = {'path', 'type'} - set(conf_support.keys())

    if missing_index:
        raise RuntimeError(
            'Missing config elements for support section {}: {}'.format(
                support_section,
                ', '.join(sorted(missing_index))
            )
        )

    if 'name' in conf_support:
        raise RuntimeError('Support section {} has reserved keyword "name"'.format(support_section))

    conf_support['name'] = support_section

    # Check path wildcards
    if '{sample}' not in conf_support['path'] and conf_support['type'] not in {'preformat', 'svpopinter-striphap'}:
        raise RuntimeError(
            'Support section {} "path" element is missing the "{{sample}}" wildcard'.format(support_section)
        )

    if conf_support['type'] == 'svpopinter-striphap':
        missing_list = [
            val  for val in ('sample_l', 'sample_r')
                if val not in conf_support['path']
        ]

        if missing_list:
            raise RuntimeError(
                'Support section {} is missing sample wildcard(s): {}'.format(support_section, ', '.join(missing_list))
            )

    if '{vartype}' not in conf_support['path']:
        raise RuntimeError(
            'Support section {} "path" element is missing the "{{vartype}}" wildcard'.format(support_section)
        )

    if '{svtype}' not in conf_support['path']:
        raise RuntimeError(
            'Support section {} "path" element is missing the "{{svtype}}" wildcard'.format(support_section)
        )

    # Make support section wildcard substitutions
    conf_support['path'] = dtablib.util.format_cards(conf_support['path'], **section_wildcards)

    # Make section wildcard substitutions (section has a "wildcard_list" element)
    if table_def['wildcards'] is not None and len(table_def['wildcards']) > 0:
        conf_support['path'] = dtablib.util.format_cards(conf_support['path'], **table_def['wildcards'])

    # Check allow-missing
    if 'allow-missing' in conf_support:
        conf_support['allow-missing'] = dtablib.util.get_bool(conf_support['allow-missing'])
    else:
        conf_support['allow-missing'] = False

    # Get column name
    if 'column-name' not in conf_support:
        conf_support['column-name'] = support_section.upper()

    # Get sample pattern - Find support only for samples matching this pattern
    if 'sample-pattern' not in conf_support:
        conf_support['sample-pattern'] = r'.*'

    conf_support['sample-pattern'] = re.compile(conf_support['sample-pattern'])

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
        hap_re = re.compile(r'^(.+)-([^-]+)$')

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
            if re.search(conf_support['sample-pattern'], sample) is None:
                continue

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

            if re.search(conf_support['sample-pattern'], sample) is None:
                continue

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
    Get support using a collection of pre-formatted tables (per sample, per svtype).

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


def get_support_table_lead_bool(id_set, input_dict, support_col_name, conf_support):
    """
    Get support using a collection of tables containing an ID column (per sample, per svtype).

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

        if 'ID' not in df_support.columns:
            raise RuntimeError('Missing columns in support tables for section {}: ID'.format(
                conf_support['name']
            ))

        df_support = df_support.loc[
            df_support['ID'].apply(lambda val: val in sample_id_set),
            ['ID']
        ]

        df_support[support_col_name] = True

        df_support = df_support.set_index('ID').reindex(sample_id_set, fill_value=False).reset_index()

        df_list.append(df_support)

    # Merge support tables
    if df_list:
        df_support = pd.concat(df_list, axis=0)
    else:
        df_support = pd.DataFrame([], columns=['ID', support_col_name])

    # Return
    return df_support
