"""
Annotations based on population metrics.
"""

# dtab_pop_table
#
# Make population table.
rule dtab_pop_table:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        tsv=expand(
            'temp/sections/{{tab_name}}/pop/pop_summary_{{vartype}}_{{svtype}}/{chrom}.tsv.gz',
            chrom=svpoplib.ref.get_df_fai(config['reference'] + '.fai').index
        )
    output:
        tsv='sections/{tab_name}/pop/pop_summary_{vartype}_{svtype}.tsv.gz'
    run:

        id_set = set(pd.read_csv(input.bed_base, sep='\t', usecols=('ID', ), squeeze=True))

        obj_list = list()

        for file_name in input.tsv:
            if os.stat(file_name).st_size > 0:
                obj_list.append(pd.read_csv(file_name, sep='\t', low_memory=False))

        if len(obj_list) == 0:
            raise RuntimeError('Found 0 data files to merge (all empty)')

        df = pd.concat(obj_list, axis=0)

        if set(df['ID']) != id_set:
            raise RuntimeError('ID set mismatch: sample_map')

        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# dtab_pop_table
#
# Make population table.
rule dtab_pop_table_chrom:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        tsv='sections/{tab_name}/base_table/pre_merge/gt_{vartype}_{svtype}.tsv.gz',
        tsv_pop=lambda wildcards: dtablib.dtabutil.sample_info_file_name(wildcards.tab_name, config)
    output:
        tsv=temp('temp/sections/{tab_name}/pop/pop_summary_{vartype}_{svtype}/{chrom}.tsv.gz')
    run:

        # Get IDs for this chrom
        id_table = pd.read_csv(input.bed_base, sep='\t', usecols=('#CHROM', 'ID'))

        id_set = set(id_table.loc[id_table['#CHROM'] == wildcards.chrom, 'ID'])

        # Write an empty output file if no records on this chromosome
        if not id_set:
            with open(output.tsv, 'w') as out_file:
                pass

            return

        # Read samples
        df = pd.read_csv(input.tsv, sep='\t', index_col='ID')
        df_pop = pd.read_csv(input.tsv_pop, sep='\t')

        # Subset to chromosome (by IDs) and check for expected IDs
        df = df.loc[[val in id_set for val in df.index]]

        missing_set = id_set - set(df.index)

        if missing_set:
            raise RuntimeError('Missing {} IDs in GT TSV: {}{}: {}'.format(
                len(missing_set),
                ', '.join(sorted(missing_set)[:3]),
                '...' if len(missing_set) > 3 else '',
                input.tsv
            ))

        # Remove no-pop samples (do not calculate population stats on these)
        if 'no_pop_samples' in config:
            pop_drop_set = set(config['no_pop_samples'])

            df_pop = df_pop[
                [val for val in df_pop.columns if val not in pop_drop_set]
            ]

        # Start list of rows
        df_list = list()

        # Get AN and AC tables
        df_an = df.applymap(lambda val: dtablib.pop.AN_TR[tuple(val.split('|'))])
        df_ac = df.applymap(lambda val: dtablib.pop.AC_TR[tuple(val.split('|'))])

        # Count AN & AC
        df_an_all = df_an.apply(np.sum, axis=1)
        df_ac_all = df_ac.apply(np.sum, axis=1)

        # Add AN and AF
        df_an_all.name = 'POP_ALL_AN'
        df_list.append(df_an_all)

        df_af_all = df_ac_all / df_an_all
        df_af_all.name = 'POP_ALL_AF'
        df_list.append(df_af_all)

        # Get rows for superpops
        for superpop in sorted(set(df_pop['SUPERPOP'])):

            # Get AN and AC for this superpop
            sample_sp = set(df_pop.loc[df_pop['SUPERPOP'] == superpop, 'SAMPLE'])
            df_an_sp = df_an.loc[:, [col for col in df_an.columns if col in sample_sp]].apply(np.sum, axis=1)
            df_an_sp.name = 'AN_SP'

            df_ac_sp = df_ac.loc[:, [col for col in df_ac.columns if col in sample_sp]].apply(np.sum, axis=1)
            df_ac_sp.name = 'AC_SP'

            # Get AN and AC for others (not-superpop)
            df_an_nsp = df_an.loc[:, [col for col in df_an.columns if col not in sample_sp]].apply(np.sum, axis=1)
            df_an_nsp.name = 'AN_NSP'

            df_ac_nsp = df_ac.loc[:, [col for col in df_ac.columns if col not in sample_sp]].apply(np.sum, axis=1)
            df_ac_nsp.name = 'AC_NSP'

            # Get AF
            df_af_sp = df_ac_sp / df_an_sp
            df_af_sp.name = 'POP_{}_AF'.format(superpop)
            df_list.append(df_af_sp)

            # Calculate FST
            fst_sp = pd.concat(
                [df_an_sp, df_ac_sp, df_an_nsp, df_ac_nsp],
                axis=1
            ).apply(
                lambda row: dtablib.fst.fst_wc(
                    row[['AC_SP', 'AC_NSP']].to_numpy(),
                    row[['AN_SP', 'AN_NSP']].to_numpy()
                ),
                axis=1
            )

            fst_sp.name = 'POP_{}_FST'.format(superpop)
            df_list.append(fst_sp)

        # Merge
        df_pop = pd.concat(df_list, axis=1)
        df_pop['POP_ALL_AN'] = df_pop['POP_ALL_AN'].astype(np.int32)
        
        # Write
        df_pop.to_csv(output.tsv, sep='\t', index=True, compression='gzip')