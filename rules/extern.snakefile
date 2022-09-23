"""
External pre-formatted variants
"""

# dtab_sup_tag_support
#
# Get a table of variant call support from another source.
rule dtab_extern_table:
    input:
        tsv=lambda wildcards: dtablib.extern.get_conf(wildcards, config)['path']
    output:
        tsv='sections/{tab_name}/extern/{vartype}_{svtype}/{extern_section}.tsv.gz'
    run:

        conf_extern = dtablib.extern.get_conf(wildcards, config)

        # Read
        if conf_extern['usecols'] is not None:
            col_list = [key for key in conf_extern['usecols'].keys()]
        else:
            col_list = None

        df = pd.read_csv(input.tsv, sep='\t', usecols=col_list, low_memory=False)

        if conf_extern['usecols'] is not None:
            df.columns = [conf_extern['usecols'][val] for val in df.columns]

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
