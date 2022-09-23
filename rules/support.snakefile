"""
Tag variants by orthogonal support.
"""


# dtab_sup_tag_support
#
# Get a table of variant call support from another source.
rule dtab_sup_tag_support:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        tsv=lambda wildcards: dtablib.support.get_input_file_list(wildcards, config)
    output:
        tsv='sections/{tab_name}/support/{vartype}_{svtype}/{support_section}.tsv.gz'
    run:

        # Get config
        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)
        conf_support = dtablib.support.get_conf(wildcards, config)

        # Read
        if table_def['sourcetype'] == 'sampleset':
            df = pd.read_csv(input.bed_base, sep='\t', usecols=('ID', 'MERGE_SRC'))
            id_set = df.groupby('MERGE_SRC')['ID'].apply(lambda id_list: set(id_list))

        elif table_def['sourcetype'] == 'caller':
            id_set = pd.Series(
                [
                    set(pd.read_csv(input.bed_base, sep='\t', usecols=('ID', ), squeeze=True))
                ],
                index=[table_def['sample']]
            )
        else:
            raise RuntimeError('sourcetype not yet supported: {sourcetype}'.format(**table_def))

        # Get input files
        input_dict = dtablib.support.get_input_dict(wildcards, config)

        # Get support table
        support_type = conf_support.get('type', None)

        if support_type is None:
            raise RuntimeError(f'Missing "type" from support config for {wildcards.support_section}')

        if support_type == 'svpopinter':
            df_support = dtablib.support.get_support_svpopinter_lead(id_set, input_dict, conf_support['column-name'])

        elif support_type == 'subseq':
            df_support = dtablib.support.get_support_subseq_lead(id_set, input_dict, conf_support['column-name'])

        elif support_type == 'rausch_bkpt':
            df_support = dtablib.support.get_support_rauschbkpt_lead(id_set, input_dict, conf_support['column-name'])

        elif support_type == 'rausch_win':
            df_support = dtablib.support.get_support_rauschwin_lead(id_set, input_dict, conf_support['column-name'])

        elif support_type == 'table':
            df_support = dtablib.support.get_support_table_lead(id_set, input_dict, conf_support['column-name'], conf_support)

        elif support_type == 'preformat':
            df_support = pd.read_csv(conf_support['path'].format(**wildcards), sep='\t', low_memory=False)

        else:
            raise RuntimeError(f'Unkown "type" from support config for {wildcards.support_section}: {support_type}')

        # Reformat columns
        if 'usecols' in conf_support:
            usecols_class = conf_support['usecols'].__class__

            if issubclass(conf_support['usecols'].__class__, dict):
                cols = ['ID'] + [col for col in df_support.columns if col in conf_support['usecols']]

                # Check for missing columns
                missing_cols = set(conf_support['usecols']) - set(cols)

                if missing_cols:
                    missing_col_list = ', '.join(sorted(missing_cols))
                    raise RuntimeError(f'Missing columns required for "usecols" in {wildcards.support_section}: {missing_col_list}')

                # Change columns
                df_support = df_support[cols]
                df_support.columns = [conf_support['usecols'].get(col, col) for col in df_support.columns]

            elif issubclass(conf_support['usecols'].__class__, list):
                cols = ['ID'] + [col for col in conf_support['usecols'] if col != 'ID']

                # Check for missing columns
                missing_cols = set(conf_support['usecols']) - set(cols)

                if missing_cols:
                    missing_col_list = ', '.join(sorted(missing_cols))
                    raise RuntimeError(f'Missing columns required for "usecols" in {wildcards.support_section}: {missing_col_list}')

                # Check for duplicates
                if len(cols) != len(set(cols)):
                    raise RuntimeError(f'Duplicate columns in list for "usecols" {wildcards.support_section}')

                # Alter columns
                df_support = df_support[cols]

            else:
                    usecols_class = str(conf_support['usecols'].__class__)
                    raise RuntimeError(f'Unkown class for "usecols" in {wildcards.support_section}: {usecols_class}')


        # Write
        df_support.to_csv(output.tsv, sep='\t', index=False, compression='gzip')
