"""
Functions for supporting external annotations
"""

import dtablib


def get_conf(wildcards, config):

    table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

    conf_extern = table_def.get(
        'extern', dict()
    ).get(
        f'{wildcards.vartype}_{wildcards.svtype}', dict()
    ).get(
        wildcards.extern_section, None
    )

    if conf_extern is None:
        raise RuntimeError(f'No extern section "{wildcards.extern_section}" in config for vartype={wildcards.vartype} and svtype={wildcards.svtype}')

    conf_extern = conf_extern.copy()

    # Check section
    missing_index = {'path'} - set(conf_extern.keys())

    if missing_index:
        raise RuntimeError(
            'Missing config elements for extern section {}: {}'.format(
                wildcards.extern_section,
                ', '.join(sorted(missing_index))
            )
        )

    if 'name' in conf_extern:
        raise RuntimeError('Extern section {} has reserved keyword "name"'.format(wildcards.support_section))

    conf_extern['name'] = wildcards.extern_section

    # Check path
    if '{vartype}' not in conf_extern['path']:
        raise RuntimeError(
            'Extern section {} "path" element is missing the "{{vartype}}" wildcard'.format(wildcards.support_section)
        )

    if '{svtype}' not in conf_extern['path']:
        raise RuntimeError(
            'Extern section {} "path" element is missing the "{{svtype}}" wildcard'.format(wildcards.support_section)
        )

    # Check usecols
    if 'usecols' in set(conf_extern.keys()):
        usecols = dict()

        if conf_extern['usecols'].__class__ == list:
            usecols = dict()

            for val in conf_extern['usecols']:
                usecols[val] = val

            conf_extern['usecols'] = usecols

        elif conf_extern['usecols'].__class__ == str:
            conf_extern['usecols'] = {conf_extern['usecols']: conf_extern['usecols']}

        elif conf_extern['usecols'].__class__ != dict:
            raise RuntimeError('Extern config usecols is not list, str, or dict: {}'.format(wildcards.extern_section))

        if 'ID' not in conf_extern['usecols'].values():
            conf_extern['usecols']['ID'] = 'ID'
    else:
        conf_extern['usecols'] = None

    # Fill path wildcards
    wildcard_dict = table_def['wildcards']
    wildcard_dict['vartype'] = wildcards.vartype
    wildcard_dict['svtype'] = wildcards.svtype

    if 'path' in conf_extern:
        conf_extern['path'] = conf_extern['path'].format(**wildcard_dict)

    # Return cconfig
    return conf_extern
