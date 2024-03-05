"""
Functions and paths for SV-Pop integration.
"""

import json
import os

import dtablib.dtabutil

import svpoplib


#
# Locations
#

# Variants
VAR_PATTERN = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/{vartype}_{svtype}.bed.gz'
FA_PATTERN = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/fa/{vartype}_{svtype}.fa.gz'


#
# Rules
#

### Base SV-Pop ###

def get_svpop_config(config):
    """
    Get SV-Pop configuration file.

    :param config: Data table config.

    :return: SV-Pop config dict.
    """

    with open(os.path.join(config['svpop_dir'], 'config/config.json')) as config_in:
        return json.load(config_in)


def resolve_rel_path(file_name, wildcards, config, param_sub=None):
    """
    Resolve a path with wildcards relative to the SV-Pop run directory.

    :param file_name: File name with wildcards relative to the SV-Pop run directory.
    :param wildcards: Rule wildcards.
    :param config: Data table config.
    :param param_sub: Dictionary for substituting `wildcards.params` or a function taking a parameter and returning
        a substituted one.

    :return: Resloved path.
    """

    table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

    fmt_dict = {}
    fmt_dict.update(table_def)
    fmt_dict.update(wildcards)

    if 'params' in fmt_dict.keys() and param_sub is not None:
        if type(param_sub) == dict:
            fmt_dict['params'] = param_sub.get(fmt_dict['params'], fmt_dict['params'])
        elif callable(param_sub):
            fmt_dict['params'] = param_sub(fmt_dict['params'])
        else:
            raise RuntimeError(f'Bad type for parameter substitution argument "param_sub": Expected dict or callable: {type(param_sub)}')

    return os.path.join(config['svpop_dir'], file_name.format(**fmt_dict))


def get_sampleset_config(table_def, config):
    """
    Get a sampleset configuration dictionary.

    :param table_def:
    :return:
    """

    if table_def['sourcetype'] != 'sampleset':
        raise RuntimeError('Cannot get sampleset config for table "{}": "sourcetype" is not "sampleset"'.format(table_def['tab_name']))

    return svpoplib.sampleset.get_config_entry(
        table_def['sourcename'], table_def['sample'], get_svpop_config(config)
    )


def sampleset_source_dict(table_def, file_pattern, config, vartype, svtype):
    """
    Get a dictionary of source input files keyed by sample names for each sample in a sampleset. These are pre-merged
    files before SV-Pop pulls them into a sampleset. They can be used to get some information for each variant in the
    set, such as the GT field from all contributing samples.

    :param table_def: Table definition.
    :param file_pattern: Pattern for source files. Wildcards in the pattern are resolved by table_def, sampleset config,
        and each sample in the sampleset config.
    :param config: Data table configuration dict.
    :param vartype: Variant type.
    :param svtype: SV type.

    :return: Dictionary of files keyed by sample names.
    """

    # Get table definition if "table_def" is a definition name and not the definition itself
    if issubclass(table_def.__class__, str):
        table_def = dtablib.dtabutil.get_table_def(table_def, config)

    # Get sampleset configuration dict
    sampleset_config = get_sampleset_config(table_def, config)

    # Setup format dictionary to read source files
    fmt_dict = table_def.copy()

    fmt_dict['sourcetype'] = sampleset_config['sourcetype']
    fmt_dict['sourcename'] = sampleset_config['sourcename']
    fmt_dict['vartype'] = vartype
    fmt_dict['svtype'] = svtype

    if 'sample' in fmt_dict:
        del(fmt_dict['sample'])

    # Resolve files
    file_dict = dict()

    for sample in sampleset_config['samples']:
        fmt_dict['sample'] = sample

        file_dict[sample] = os.path.join(config['svpop_dir'], file_pattern.format(**fmt_dict))

    return file_dict
