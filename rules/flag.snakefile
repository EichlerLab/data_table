# Flag rules for running multiple jobs

import dtablib

global touch

def _dtab_flag_input_support_matrix(wildcards):
    table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

    support_section = table_def.get('support', {}).get(f'{wildcards.vartype}_{wildcards.svtype}', {})

    sec_list = [sec_name for sec_name in support_section if support_section[sec_name].get('type', None) in {'svpopinter', 'svpopinter-striphap'}]

    return [
        'sections/{tab_name}/support_matrix/{vartype}_{svtype}/{support_section}.tsv.gz'.format(
            tab_name=wildcards.tab_name, vartype=wildcards.vartype, svtype=wildcards.svtype, support_section=sec_name
        )
            for sec_name in sec_list
    ]

rule dtab_flag_support_matrix:
    input:
        _dtab_flag_input_support_matrix
    output:
        txt=touch('flag/sections/{tab_name}/support_matrix/{vartype}_{svtype}/all.txt')
