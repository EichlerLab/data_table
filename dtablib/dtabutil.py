"""
Get and process table definition configs.
"""

import os
import gzip

import dtablib.svpop


SECTION_DICT = {
    'refseq': ['sections/{tab_name}/genes/refseq-intersect_{vartype}_{svtype}.tsv.gz'],
    'refseq_prox': lambda params: ['sections/{{tab_name}}/genes/refseq-prox-{params}_{{vartype}}_{{svtype}}.tsv.gz'.format(params=params)],
    'ref_sd': ['sections/{tab_name}/regions/sd_{vartype}_{svtype}.tsv.gz'],
    'ref_trf': ['sections/{tab_name}/regions/trf_{vartype}_{svtype}.tsv.gz'],
    'encode': ['sections/{tab_name}/regions/encode_{vartype}_{svtype}.tsv.gz'],
    'dhs2020': ['sections/{tab_name}/regions/dhs2020_{vartype}_{svtype}.tsv.gz'],
    'ccre2020': ['sections/{tab_name}/regions/ccre2020_{vartype}_{svtype}.tsv.gz'],
    'oreganno': ['sections/{tab_name}/regions/oreganno_{vartype}_{svtype}.tsv.gz'],
    'win': lambda params: ['sections/{{tab_name}}/regions/ref_win_{params}_{{vartype}}_{{svtype}}.tsv.gz'.format(params=params)],
    'band': ['{svpop_run_dir}/results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/bands/bands_{vartype}_{svtype}.tsv.gz'],
    'pop': lambda params: ['sections/{{tab_name}}/pop/pop_{params}_{{vartype}}_{{svtype}}.tsv.gz'.format(params=params)],
    'rmsk': ['sections/{tab_name}/interseq/rmsk_{vartype}_{svtype}.tsv.gz'],
    'trf': ['sections/{tab_name}/interseq/trf_{vartype}_{svtype}.tsv.gz'],
    'gc': ['sections/{tab_name}/interseq/gc_{vartype}_{svtype}.tsv.gz'],
    'tandem': ['sections/{tab_name}/mapping/tandem-dup_{vartype}_{svtype}.tsv.gz'],
    'homop': ['sections/{tab_name}/regions/homop-intersect_{vartype}_{svtype}.tsv.gz'],
    'homop_nearest': ['sections/{tab_name}/regions/homop-nearest_{vartype}_{svtype}.tsv.gz'],
    'dinucl': ['sections/{tab_name}/regions/dinucl-intersect_{vartype}_{svtype}.tsv.gz'],
    'dinucl_nearest': ['sections/{tab_name}/regions/dinucl-nearest_{vartype}_{svtype}.tsv.gz'],
    'snv_trans': ['sections/{tab_name}/snv/transitions_{vartype}_{svtype}.tsv.gz'],
    'geneset': lambda params: ['sections/{{tab_name}}/geneset/{params}_{{vartype}}_{{svtype}}.tsv.gz'.format(params=params)],
}

SECTION_TIER = {
    'geneset': 2
}

def get_table_def(tab_name, config):
    """
    Get a table definition and check for required fields.

    :param tab_name: Table definition name (under "table_def" in config).
    :param config: Config dict.

    :return: Table definition dict.
    """

    # Split table name into name and wildcards ('-' delimited)
    tab_name_org = tab_name
    tab_name_tok = tab_name.split('-')

    config_section_name = tab_name_tok[0]
    wildcard_list = tab_name_tok[1:]

    # Find table definition
    if 'table_def' not in config:
        raise RuntimeError('Config is missing section: table_def')

    if 'svpop_run_dir' not in config:
        raise RuntimeError('Config is missing section: svpop_run_dir')

    if config_section_name not in config['table_def']:
        raise RuntimeError(f'No table_def config for table definition: {config_section_name}')

    table_def = config['table_def'][config_section_name].copy()

    # Check for missing fields
    missing_fields = [field for field in ('sourcetype', 'sourcename', 'sample', 'filter', 'svset', 'sections') if field not in table_def]

    if missing_fields:
        raise RuntimeError('Table definition "{}" is missing field(s): {}'.format(config_section_name, ', '.join(missing_fields)))

    # Set name
    if 'tab_name' in table_def:
        raise RuntimeError(f'Table configuration far "{config_section_name}" contains reserved field "tab_name"')

    table_def['tab_name'] = tab_name
    table_def['config_section_name'] = config_section_name

    # Set wildcards
    if wildcard_list:

        # Check lists
        if 'wildcard_list' not in table_def:
            raise RuntimeError('Table name "{}" has {} wildcards ("-" delimited), but table definition "{}" does not contain a "wildcard_list" element'.format(
                tab_name_org, len(wildcard_list), config_section_name
            ))

        if not issubclass(table_def['wildcard_list'].__class__, list):
            raise RuntimeError('Table definition "{}" has a "wildcard_list" element that is not a list'.format(
                config_section_name
            ))

        if len(wildcard_list) != len(table_def['wildcard_list']):
            raise RuntimeError(
                'Wildcard length mismatch for {}: Found {} wildcards in table name but definition "{}" has {} elements'.format(
                    tab_name_org, len(wildcard_list), config_section_name, len(table_def['wildcards'])
                ))

        # Set dictionary
        table_def['wildcards'] = {
            attribute: value for attribute, value in zip(table_def['wildcard_list'], wildcard_list)
        }

        # Do parameter substitutions
        table_def['sample'] = table_def['sample'].format(**table_def['wildcards'])
        table_def['sourcename'] = table_def['sourcename'].format(**table_def['wildcards'])

    else:

        # Check: Definition should have no wildcards list
        if 'wildcard_list' in table_def:
            raise RuntimeError('Table name "{}" has 0 wildcards ("-" delimited), but table definition "{}" contains a "wildcard_list" element'.format(
                tab_name_org, config_section_name
            ))

        # Set empty list
        table_def['wildcards'] = dict()

    # Return table definition
    return table_def


def get_id_list(file_name):
    """
    Get a list of IDs written by base_table (avoids loading the whole dataframe).

    :params file_name: ID list file name.
    """

    if not os.path.isfile(file_name):
        raise RuntimeError(f'Missing ID list file: {file_name}')

    with gzip.open(file_name, 'rt') as in_file:
        id_list = [val.strip() for val in in_file]

    id_list = [val for val in id_list if val]

    return id_list


def section_input_files(wildcards, config, max_tier=None):
    """
    Get input files needed to satisfy each section.

    :param wildcards: Rule wildcards.
    :param config: Config.
    :param max_tier: Get maximum tier allowed. Sections are grouped into tiers. Most are tier 1 (default). Tier 2
        requires tier 1 to complete first (e.g. gene sets). If `None`, get input files for all sections.

    :return: List of input files.
    """

    file_list = list()  # List of files to merge

    table_def = get_table_def(wildcards.tab_name, config)

    varsv_type = '{vartype}_{svtype}'.format(**wildcards)

    if varsv_type not in table_def['sections']:
        raise RuntimeError(f'Cannot find configuration section "{varsv_type}"')

    # Add keywords to dict
    kwd_dict = dict(wildcards)
    kwd_dict.update(table_def)
    kwd_dict['svpop_run_dir'] = config['svpop_run_dir']

    # Get files for support sections
    for support_section in table_def.get('support', {}).get(varsv_type, {}):
        kwd_dict['support_section'] = support_section
        file_list.append('sections/{tab_name}/support/{vartype}_{svtype}/{support_section}.tsv.gz'.format(**kwd_dict))

    # Get files for annotation sections
    for section in table_def['sections'][varsv_type]:

        # Split section into section name and params ("-" delimited)
        section_tok = section.split('-', 1)

        section_name = section_tok[0]
        section_param = section_tok[1] if len(section_tok) > 1 else None

        if section_name not in SECTION_DICT:
            raise RuntimeError('Unknown section "{}" in table definition: {}'.format(section, wildcards.tab_name))

        # Skip if section is is a higher tier than max
        if max_tier is not None and SECTION_TIER.get(section_name, 1) > max_tier:
            continue

        # Get section values (string with wildcards, list of strings with wildcards, or a function returning either a
        # string or list of strings with wildcards.
        section_val = SECTION_DICT[section_name]

        # Get paths if section is defined by a function
        if callable(section_val):
            try:
                section_val = section_val(section_param)
            except Exception as ex:
                raise RuntimeError('Error calling function for section {} (param = "{}"): {}'.format(section, section_param, ex))

        elif section_param is not None:
            raise RuntimeError('Cannot resolve a section with parameters (values after "-") if section is not defined by a function: {}'.format(section))

        # Translate to a list if not already
        if not issubclass(section_val.__class__, list):
            if not issubclass(section_val.__class__, str):
                raise RuntimeError('Section definition must be string or list of string (or a function returning these): {}: class = {}'.format(section, section_val.__class__))

            section_val = [section_val]

        # Fill wildcards and append to list
        for file_name in section_val:
            file_list.append(file_name.format(**kwd_dict))

    # Get external files
    for extern_section in table_def.get('extern', {}).get(varsv_type, {}):
        kwd_dict['extern_section'] = extern_section
        file_list.append('sections/{tab_name}/extern/{vartype}_{svtype}/{extern_section}.tsv.gz'.format(**kwd_dict))

    # Return file list
    return file_list


def sample_info_file_name(tab_name, config):

    table_def = get_table_def(tab_name, config)

    if 'sample_info' in table_def:
        sample_info_file = table_def['sample_info']
        sample_info_source = 'table def "{}"'.format(tab_name)

    elif 'sample_info' in config:
        sample_info_file = config['sample_info']
        sample_info_source = 'global config'.format(tab_name)

    else:
        raise RuntimeError('No "sample_info" section found in table definition or global config')

    if not os.path.isfile(sample_info_file):
        raise RuntimeError(f'Sample info table in {sample_info_source} is not a file: {sample_info_file}')

    return sample_info_file


def get_callable_bed_dict(tab_name, config, hap):
    """
    Get a dictionary keyed by samples to a list of BED files from "callable_filter".

    :param tab_name: Target table name.
    :param config: Config.
    :param hap: Haplotype. Parses into "{hap}" keyword in "callable_filter".

    :return: Dictionary of BED file names for callable regions keyed by the sample.
    """

    # Get table definition
    table_def = get_table_def(tab_name, config)

    # Get BED pattern
    if 'callable_filter' not in table_def:
        raise RuntimeError(
            'Cannot get callable filter file: no "callable_filter" section in table definition "{}"'.format(tab_name)
        )

    bed_pattern = table_def['callable_filter']

    # Get sampleset config
    sampleset_config = dtablib.svpop.get_sampleset_config(table_def, config)

    # Setup parse_dict (wildcard parsing into bed_pattern)
    parse_dict = config.copy()
    parse_dict.update(table_def)
    parse_dict['hap'] = hap
    parse_dict['sample'] = None

    # Parse
    bed_dict = dict()

    for sample in sampleset_config['samples']:
        parse_dict['sample'] = sample

        bed_dict[sample] = bed_pattern.format(**parse_dict)

    # Return
    return bed_dict
