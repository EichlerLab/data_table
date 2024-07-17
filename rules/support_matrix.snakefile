# Get matrix of support through a sampleset merge

import pandas as pd

import dtablib


# Generate a support matrix for a support source. Table of variant IDs (rows) by samples (columns) with variant IDs from
# the support source (or empty if no support).
rule dtab_sampleset_support_matrix:
    input:
        merge_map='sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}.tsv.gz',
        tsv=lambda wildcards: dtablib.support.get_input_file_list(wildcards, config)
    output:
        tsv='sections/{tab_name}/support_matrix/{vartype}_{svtype}/{support_section}.tsv.gz'
    run:

        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)
        conf_support = dtablib.support.get_conf(wildcards, config)

        if table_def['sourcetype'] != 'sampleset':
            raise RuntimeError(f'Source type must be "sampleset" in table definition {table_def["tab_name"]}: Found "{table_def["sourcetype"]}"')

        if conf_support['type'] not in {'svpopinter', 'svpopinter-striphap'}:
            raise RuntimeError(f'Support type must be "svpopinter" or "svpopinter-striphap" for support section {conf_support["name"]}: Found "{conf_support["type"]}"')

        input_dict = dtablib.support.get_input_dict(wildcards, config)

        # Read merge map
        df_map = pd.read_csv(
            input.merge_map, sep='\t', index_col='ID',
            dtype={'#CHROM': str}, low_memory=False
        )

        # Read each column
        df_list = list()

        for sample in df_map.columns:

            # Read support table
            if len(input_dict[sample]) > 0:
                df_support = pd.concat(
                    [pd.read_csv(in_file_name, sep='\t', usecols=['ID_A', 'ID_B']) for in_file_name in input_dict[sample]],
                    axis=0
                )

            else:
                df_support = pd.DataFrame([], columns=['ID_A', 'ID_B'])

            df_support = df_support.loc[
                (
                    ~ pd.isnull(df_support['ID_A'])
                ) & (
                    ~ pd.isnull(df_support['ID_B'])
                )
            ]

            df_support = df_support.set_index('ID_A')['ID_B']
            df_support.name = sample

            # Translate variant ID from sample to merge
            id_tr = df_map.loc[:, sample].reset_index()
            #id_tr = id_tr.loc[~ pd.isnull(id_tr[sample])].set_index(sample).squeeze()
            id_tr = id_tr.loc[~ pd.isnull(id_tr[sample])].set_index(sample)['ID']

            df_support.index = [id_tr[val] if val in id_tr.index else np.nan for val in df_support.index]
            df_support = df_support.loc[~ pd.isnull(df_support.index)]

            # Add to support list
            df_list.append(df_support)

        df = pd.concat(df_list, axis=1).reindex(df_map.index)

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')
