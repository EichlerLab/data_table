"""
Build base variant table.
"""

#
# Definitions
#

def _dtab_base_table_chrom_mem(wildcards):
    """
    Get memory requirements for the initial (per chromosome) data table resources.

    :param wildcards: Rule wildcards.

    :return: Memory usage.
    """

    if wildcards.vartype == 'sv':
        return '6G'

    if wildcards.vartype == 'indel':
        return '16G'

    if wildcards.vartype == 'snv':
        return '64G'

    # Unknown
    return '8G'

def filter_fasta(fa_file_name, id_set, verify_all=True):
    """
    Read records from a FASTA file, `fa_file_name`, and keep only

    :param fa_file_name: Input FASTA file.
    :param id_set: Set of IDs to extract.
    :param verify_all: If `True`, throw an exception if any IDs in `id_set` were not found`.

    :return: An iterator for `SeqRecord` objects.
    """

    id_set_found = set()

    # Iterate over dataframe and return SeqRecords
    for index, row in df.iterrows():

        with gzip.open(input.fa, 'rt') as in_file:
            for record in Bio.SeqIO.parse(in_file, 'fasta'):

                if record.id in id_set:

                    if record.id in id_set_found:
                        raise RuntimeError('Duplicate ID found in input: {}'.format(record.id))

                    id_set_found.add (record.id)



#
# Rules
#

# rule dtab_base_fa:
#     input:
#         tsv='tsv/variants_{tab_name}_{vartype}_{svtype}.tsv.gz',
#         fa=lambda wildcards: dtablib.base_table.get_input_fa(wildcards, config)
#     output:
#         fa='fasta/variants_{tab_name}_{vartype}_{svtype}.fa.gz'
#     run:
#
#         # Read ID set
#         id_set = set(pd.read_csv(input.tsv, sep='\t', usecols=('ID', ), squeeze=True))
#
#         # Process FASTA
#         with gzip.open(input.fa, 'rt') as in_file:
#             for record in Bio.SeqIO.parse(in_file, 'fasta'):

# dtab_base_xlsx
#
# Table to Excel.
rule dtab_base_xlsx:
    input:
        tsv='tsv/variants_{tab_name}_{vartype}_{svtype}.tsv.gz'
    output:
        xlsx='tsv/variants_{tab_name}_{vartype}_{svtype}.xlsx'
    run:

        df = pd.read_csv(
            input.tsv, sep='\t'
        ).to_excel(
            output.xlsx, index=False
        )

# dtab_base_merge
#
# Merge final BED
rule dtab_base_merge:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        tsv_sec=lambda wildcards: dtablib.dtabutil.section_input_files(wildcards, config)
    output:
        tsv='tsv/variants_{tab_name}_{vartype}_{svtype}.tsv.gz'
    run:

        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

        # Read base
        df = pd.read_csv(input.bed_base, sep='\t', index_col='ID', low_memory=False)

        # Read sections
        input_list = list(input.tsv_sec)

        df_list = list()

        for file_name in input.tsv_sec:
            print(file_name)  # DBGTMP
            df_anno = pd.read_csv(file_name, sep='\t', low_memory=False)

            if 'ID' not in df_anno.columns:
                raise RuntimeError('Missing ID column: {}'.format(file_name))

            if df_anno.shape[1] == 1:
                raise RuntimeError('No annotation fields (only ID column): {}'.format(file_name))

            df_anno.set_index('ID', inplace=True)

            df_list.append(df_anno.reindex(df.index))

        df_anno = pd.concat(df_list, axis=1).reindex(df.index)

        df = pd.concat([df, df_anno], axis=1)

        del(df_anno)

        # Drop columns
        if 'drop_cols' in table_def:
            for col in table_def['drop_cols']:
                if col in df.columns:
                    del(df[col])

        # Rename columns
        if 'rename_cols' in table_def:
            df.columns = [table_def['rename_cols'].get(col, col) for col in df.columns]

        # Move "WIN_*" to end
        tail_cols = [col for col in df.columns if col.startswith('WIN_')]
        head_cols = [col for col in df.columns if col not in tail_cols]

        df = df[head_cols + tail_cols]

        # Shift columns.
        # Config section "shift_col" is a list of dicts. For each dict, take columns in list "shift_col" and move to
        # after "move_cols". Repeat for each element in "shift_col".
        if 'shift_col' in table_def:
            for shift_dict in table_def['shift_col']:
                insert_after = shift_dict['insert_after']
                move_cols = [col for col in shift_dict['move_cols'] if col in df.columns]

                if move_cols:
                    pre_post_cols = [
                        list(),
                        list()
                    ]

                    state = 0

                    for col in df.columns:
                        if col not in move_cols:
                            pre_post_cols[state].append(col)

                        if col == insert_after:
                            state = 1

                    df = df[pre_post_cols[0] + move_cols + pre_post_cols[1]]

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')
#        df.to_excel(output.xlsx, index=True)

# dtab_base_premerge_gt
#
# Get genotypes from pre-merged variants supporting each merged variant.
rule dtab_base_premerge_gt:
    input:
        id_list='sections/{tab_name}/base_table/id_list_{vartype}_{svtype}.txt.gz',
        callable_h1='sections/{tab_name}/base_table/callable/callable_{vartype}_{svtype}_h1.tsv.gz',
        callable_h2='sections/{tab_name}/base_table/callable/callable_{vartype}_{svtype}_h2.tsv.gz',
        tsv=lambda wildcards:
        [
            'temp/{tab_name}/base_table/pre_merge/split/gt_{vartype}_{svtype}_{{svtype_single}}.tsv.gz'.format(
                **wildcards
            ).format(
                svtype_single=svtype_single
            )
                for svtype_single in dtablib.definitions.SVTYPE_EXPAND[wildcards.svtype]
        ]
    output:
        tsv='sections/{tab_name}/base_table/pre_merge/gt_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        id_list = dtablib.dtabutil.get_id_list(input.id_list)

        # Merge GT tables
        df_list = [pd.read_csv(file_name, sep='\t', index_col='ID') for file_name in input.tsv]

        if len(df_list) > 1:
            df = pd.concat(df_list, axis=0)
        else:
            df = df_list[0]

        # Read callable regions
        df_callable_h1 = pd.read_csv(input.callable_h1, sep='\t', index_col='ID')
        df_callable_h2 = pd.read_csv(input.callable_h2, sep='\t', index_col='ID')

        # Check columns and variants
        missing_cols_h1 = set(df.columns) - set(df_callable_h1.columns)
        missing_cols_h2 = set(df.columns) - set(df_callable_h2.columns)

        if missing_cols_h1 or missing_cols_h2:
            missing_cols = missing_cols_h1 | missing_cols_h2

            raise RuntimeError('Missing {} column(s) in callable regions: {}{}'.format(
                len(missing_cols), ', '.join(sorted(missing_cols)[:3]), '...' if len(missing_cols) > 3 else ''
            ))

        missing_id_h1 = set(df.index) - set(df_callable_h1.index)
        missing_id_h2 = set(df.index) - set(df_callable_h2.index)

        if missing_id_h1 or missing_id_h2:
            missing_id = missing_id_h1 | missing_id_h2

            raise RuntimeError('Missing {} variant ID(s) in callable regions: {}{}'.format(
                len(missing_id), ', '.join(sorted(missing_id)[:3]), '...' if len(missing_id) > 3 else ''
            ))

        df_callable_h1 = df_callable_h1.loc[df.index, df.columns]
        df_callable_h2 = df_callable_h2.loc[df.index, df.columns]

        # Get genotype tables
        df_gt_h1 = df.apply(lambda row:
            row[~ pd.isnull(row)].apply(lambda val: val.split('|')[0]).reindex(row.index),
           axis=0
        )

        df_gt_h2 = df.apply(lambda row:
            row[~ pd.isnull(row)].apply(lambda val: val.split('|')[1]).reindex(row.index),
           axis=0
        )

        # Fill 0 or . for missing values
        callable_dict = {
            'h1': df_callable_h1,
            'h2': df_callable_h2
        }

        gt_dict = {
            'h1': df_gt_h1,
            'h2': df_gt_h2
        }

        for col_name in df.columns:
            for hap in ('h1', 'h2'):

                # Get a genotype column and subset to null values
                col = gt_dict[hap][col_name]
                col = col.loc[pd.isnull(col)]

                # Get callable loci and subset to variants in "col"
                col_callable = callable_dict[hap][col_name].loc[col.index]
                col_callable = col_callable.loc[col_callable]  # Subset to True values

                # Set "0" on callable loci in col (was already subset to NA values)
                col[col_callable.index] = '0'

                # Fill remaining values with no-call
                col.fillna('.', inplace=True)

                # Assign missing values back to the haplotype table
                gt_dict[hap][col_name].loc[col.index] = col

        # Make genotype table
        df_gt_list = list()

        for col_name in df.columns:
            # Join h1 and h2 genotypes into a table (one column each) and concatenate.
            gt_col = pd.concat(
                [gt_dict['h1'][col_name], gt_dict['h2'][col_name]],
                axis=1
            ).apply(
                lambda vals: '|'.join(vals), axis=1
            )

            gt_col.name = col_name

            df_gt_list.append(gt_col)

        df = pd.concat(df_gt_list, axis=1)

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_base_premerge_gt_by_svtype
#
# Make genotype table for source (pre-merged) sampleset variants.
rule dtab_base_premerge_gt_by_svtype:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        merge_map='sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}.tsv.gz',
        tsv_sample=lambda wildcards: dtablib.svpop.sampleset_source_dict(
            wildcards.tab_name, dtablib.svpop.VAR_PATTERN, config, wildcards.vartype, wildcards.svtype_single
        ).values()
    output:
        tsv=temp('temp/{tab_name}/base_table/pre_merge/split/gt_{vartype}_{svtype}_{svtype_single}.tsv.gz')
    wildcard_constraints:
        svtype_single='ins|del|inv|snv|dup|sub|rgn'
    run:

        # Get ID list
        df_base = pd.read_csv(input.bed_base, sep='\t', usecols=('ID', 'SVTYPE'))

        id_list = list(df_base.loc[df_base['SVTYPE'] == wildcards.svtype_single.upper(), 'ID'])

        # Get dict of input files
        source_file_dict = dtablib.svpop.sampleset_source_dict(
            wildcards.tab_name, dtablib.svpop.VAR_PATTERN, config, wildcards.vartype, wildcards.svtype_single
        )

        # Read merge map
        df_merge = pd.read_csv(input.merge_map, sep='\t', index_col='ID')

        df_merge = df_merge.loc[id_list]

        # Define a function to translate a column in df_merge from variant IDs to genotypes
        def to_gt(sample_col):

            sample_col = sample_col.dropna()

            # Get sample table
            df_sample = pd.read_csv(source_file_dict[sample_col.name], sep='\t', index_col='ID')

            if 'GT' not in df_sample.columns:
                raise RuntimeError('Missing GT column for sample "{}": {}'.format(
                    sample_col.name,
                    source_file_dict[sample_col.name]
                ))

            # Check for variants missing that the merge set says should come from this sample
            missing_set = set(sample_col) - set(df_sample.index)

            if missing_set:
                raise RuntimeError('Missing {} variants for sample "{}": {}{}: {}'.format(
                    len(missing_set),
                    sample_col.name,
                    ', '.join(sorted(missing_set)[:3]),
                    '...' if len(missing_set) > 3 else '',
                    source_file_dict[sample_col.name]
                ))

            sample_gt = df_sample.loc[list(sample_col), 'GT']
            sample_gt.index = list(sample_col.index)

            return sample_gt

        df_merge_gt = df_merge.apply(to_gt, axis=0)

        df_merge_gt.index.name = 'ID'

        # Write
        df_merge_gt.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_base_callable
#
# Get a table of callable variants by variant ID (rows) and sample (columns). Values are True/False.
rule dtab_base_callable:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        bed_callable=lambda wildcards: dtablib.dtabutil.get_callable_bed_dict(wildcards.tab_name, config, wildcards.hap).values()
    output:
        tsv='sections/{tab_name}/base_table/callable/callable_{vartype}_{svtype}_{hap}.tsv.gz'
    run:

        # Read
        df = pd.read_csv(input.bed_base, sep='\t', usecols=('#CHROM', 'POS', 'END', 'ID', 'SVLEN'), index_col='ID')

        callable_dict = dtablib.dtabutil.get_callable_bed_dict(wildcards.tab_name, config, wildcards.hap)

        # Make callable table
        df_callable_list = list()

        for sample in callable_dict.keys():
            df_sample = pd.read_csv(callable_dict[sample], sep='\t')

            region_tree = dtablib.util.bed_to_tree_dict(df_sample)

            sample_callable = dtablib.util.bed_match_dict_region(df, region_tree)
            sample_callable.name = sample

            df_callable_list.append(sample_callable)

        df_callable = pd.concat(df_callable_list, axis=1)

        # Write
        df_callable.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_base_table_merge_table
#
# Merge base tables.
rule dtab_base_table_merge_table:
    input:
        bed_base=expand(
            'temp/sections/{{tab_name}}/base_table/base_table_{{vartype}}_{{svtype}}/{chrom}.bed.gz',
            chrom=svpoplib.ref.get_df_fai(config['reference'] + '.fai').index
        ),
        id_list=expand(
            'temp/sections/{{tab_name}}/base_table/id_list_{{vartype}}_{{svtype}}/{chrom}.txt.gz',
            chrom=svpoplib.ref.get_df_fai(config['reference'] + '.fai').index
        ),
        sample_pkl=expand(
            'temp/sections/{{tab_name}}/base_table/sample_set_{{vartype}}_{{svtype}}/{chrom}.pkl',
            chrom=svpoplib.ref.get_df_fai(config['reference'] + '.fai').index
        ),
        merge_map=expand(
            'temp/sections/{{tab_name}}/base_table/merge_map_{{vartype}}_{{svtype}}/{chrom}.tsv.gz',
            chrom=svpoplib.ref.get_df_fai(config['reference'] + '.fai').index
        ),
        tsv_id_table=lambda wildcards: dtablib.dtabutil.get_table_def(wildcards.tab_name, config)['id_table'] if 'id_table' in dtablib.dtabutil.get_table_def(wildcards.tab_name, config) else []
    output:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        id_list='sections/{tab_name}/base_table/id_list_{vartype}_{svtype}.txt.gz',
        sample_pkl='sections/{tab_name}/base_table/sample_set_{vartype}_{svtype}.pkl',
        merge_map='sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}.tsv.gz'
    run:

        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

        ### ID list ###
        id_set = set()

        with gzip.open(output.id_list, 'wt') as out_file:
            for file_name in input.id_list:
                if os.stat(file_name).st_size > 0:
                    with gzip.open(file_name, 'rt') as in_file:
                        for line in in_file:
                            out_file.write(line)
                            id_set.add(line.strip())

        ### Sample set PKL ###
        obj_list = list()

        for file_name in input.sample_pkl:
            if os.stat(file_name).st_size > 0:
                with open(file_name, 'rb') as in_file:
                    obj_list.append(pickle.load(in_file))

        if len(obj_list) == 0:
            raise RuntimeError('Found 0 data files to merge (all empty): sample_set')

        sample_set = pd.concat(obj_list, axis=0)

        if set(sample_set.index) != id_set:
            raise RuntimeError('ID set mismatch: sample_set')

        sample_set.to_pickle(output.sample_pkl)

        del sample_set

        ### Merge Map ###
        obj_list = list()

        for file_name in input.merge_map:
            if os.stat(file_name).st_size > 0:
                obj_list.append(pd.read_csv(file_name, sep='\t', low_memory=False))

        if len(obj_list) == 0:
            raise RuntimeError('Found 0 data files to merge (all empty): merge_map')

        df_sample_map = pd.concat(obj_list, axis=0)

        if set(df_sample_map['ID']) != id_set:
            raise RuntimeError('ID set mismatch: sample_map')

        df_sample_map = df_sample_map[['ID'] + [val for val in df_sample_map.columns if val != 'ID']]

        df_sample_map.to_csv(output.merge_map, sep='\t', index=False, compression='gzip')

        ### BED Base ###
        obj_list = list()

        for file_name in input.bed_base:
            if os.stat(file_name).st_size > 0:
                obj_list.append(pd.read_csv(file_name, sep='\t', low_memory=False))

        if len(obj_list) == 0:
            raise RuntimeError('Found 0 data files to merge (all empty): bed_base')

        df = pd.concat(obj_list, axis=0).sort_values(['#CHROM', 'POS'])

        if set(df['ID']) != id_set:
            raise RuntimeError('ID set mismatch: sample_map')

        # Check match against id_table
        if 'id_table' in table_def:

            # Read set of IDs to retain
            id_set = set(pd.read_csv(
                table_def['id_table'].format(**wildcards),
                sep='\t', usecols=('ID', ), squeeze=True
            ))

            # If retained ID set contains IDs not in the unfiltered table, error
            id_set_missing = id_set - set(df['ID'])

            if id_set_missing:
                raise RuntimeError('Filter ID set contains {} IDs not in the original table: {}{})'.format(
                    len(id_set_missing),
                    ','.join(sorted(id_set_missing)[:3]),
                    '...' if len(id_set_missing) > 3 else ''
                ))

            del id_set

        df.to_csv(output.bed_base, sep='\t', index=False, compression='gzip')

# dtab_base_table
#
# Make base table annotation columns will be attached to.
rule dtab_base_table_chrom:
    input:
        bed=lambda wildcards: dtablib.base_table.get_input_bed(wildcards, config)
    output:
        bed_base=temp('temp/sections/{tab_name}/base_table/base_table_{vartype}_{svtype}/{chrom}.bed.gz'),
        id_list=temp('temp/sections/{tab_name}/base_table/id_list_{vartype}_{svtype}/{chrom}.txt.gz'),
        sample_pkl=temp('temp/sections/{tab_name}/base_table/sample_set_{vartype}_{svtype}/{chrom}.pkl'),
        merge_map=temp('temp/sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}/{chrom}.tsv.gz')
    params:
        mem=_dtab_base_table_chrom_mem
    run:

        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

        sourcetype = table_def['sourcetype']

        # Read
        df_list = list()

        df_iter = pd.read_csv(input.bed, sep='\t', low_memory=False, iterator=True, chunksize=20000)

        for df in df_iter:
            df = df.loc[df['#CHROM'] == wildcards.chrom]

            if df.shape[0] > 0:
                df_list.append(df)

        # Write empty files for 0-record chromosomes
        if len(df_list) == 0:
            with open(output.bed_base, 'w'):
                pass
            with open(output.id_list, 'w'):
                pass
            with open(output.sample_pkl, 'w'):
                pass
            with open(output.merge_map, 'w'):
                pass

            return

        # Merge data frames
        df = pd.concat(df_list, axis=0)

        df.set_index('ID', inplace=True, drop=False)
        df.index.name = 'INDEX'

        # Subset variants from a table of variants (with ID column).
        if 'id_table' in table_def:

            # Read set of IDs to retain
            id_set = set(pd.read_csv(
                table_def['id_table'].format(**wildcards),
                sep='\t', usecols=('ID', ), squeeze=True
            ))

            # Subset to IDs
            df = df.loc[df['ID'].apply(lambda val: val in id_set)]

            # Stop if no variants are left
            if df.shape[0] == 0:
                with open(output.bed_base, 'w'):
                    pass
                with open(output.id_list, 'w'):
                    pass
                with open(output.sample_pkl, 'w'):
                    pass
                with open(output.merge_map, 'w'):
                    pass

                return

            del id_set

        # Make a set of samples
        if 'MERGE_SAMPLES' in df.columns:
            sample_set = df['MERGE_SAMPLES'].apply(lambda val: set(val.split(',')))
        else:
            if 'wildcards' in table_def and 'sample' in table_def['wildcards']:
                sample_set = pd.Series([{table_def['wildcards']['sample']}] * df.shape[0], index=df.index)
            else:
                sample_set = pd.Series([{}] * df.shape[0], index=df.index)

        # Make a sample-map table (each merged variant ID to an original sample/ID pair)
        # Result is a table with merged variant IDs (rows), each sample (columns), with the original variant ID from
        # the corresponding sample as the data.
        if 'MERGE_SAMPLES' in df.columns:
            sample_map_list = list()

            for index, row in df.iterrows():
                for sample, org_id in zip(row['MERGE_SAMPLES'].split(','), row['MERGE_VARIANTS'].split(',')):
                    sample_map_list.append(pd.Series([row['ID'], sample, org_id], index=['ID', 'SAMPLE', 'ORG_ID']))

            df_sample_map = pd.concat(
                sample_map_list, axis=1
            ).T.pivot(
                index='ID', columns='SAMPLE', values='ORG_ID'
            )

            # Clear variants if all samples do not match the accepted sample pattern (used for removing replicates)
            if 'sample_pattern' in table_def:
                pattern = re.compile(table_def['sample_pattern'])

                sampleset_filter = sample_set.apply(lambda sample_set:
                        {
                            bool(re.match(pattern, sample)) for sample in sample_set
                        } != {False}
                )

                df = df.loc[sampleset_filter]
                sample_set = sample_set.loc[sampleset_filter]

        else:
            df_sample_map = df['ID']
            df_sample_map.index.name = 'ID'
            df_sample_map.name = table_def['sample']

            df_sample_map = pd.DataFrame(df_sample_map).reset_index()

        # Order variant columns
        df = svpoplib.variant.order_variant_columns(df, allow_missing=True)


        # Write
        with gzip.open(output.id_list, 'wt') as out_file:
            for var_id in df.index:
                out_file.write(var_id)
                out_file.write('\n')

        sample_set.to_pickle(output.sample_pkl)

        df.to_csv(output.bed_base, sep='\t', index=False, compression='gzip')

        df_sample_map.to_csv(output.merge_map, sep='\t', index=True, compression='gzip')
