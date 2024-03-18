"""
Build base variant table.
"""

#
# Definitions
#

CHROM_PARTITIONS = 20

def _dtab_base_table_part_mem(wildcards):
    """
    Get memory requirements for the initial (per chromosome partition) data table resources.

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
#         tsv='table/variants_{tab_name}_{vartype}_{svtype}.tsv.gz',
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
        tsv='table/variants_{tab_name}_{vartype}_{svtype}.tsv.gz'
    output:
        xlsx='table/variants_{tab_name}_{vartype}_{svtype}.xlsx'
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
        tsv='table/variants_{tab_name}_{vartype}_{svtype}.tsv.gz'
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
        callable='sections/{tab_name}/base_table/callable/callable_{vartype}_{svtype}.tsv.gz',
        merge_map='sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}.tsv.gz'
    output:
        tsv_sample='sections/{tab_name}/base_table/pre_merge/gt_sample_{vartype}_{svtype}.tsv.gz',
        tsv_hap='sections/{tab_name}/base_table/pre_merge/gt_hap_{vartype}_{svtype}.tsv.gz'
    run:

        # Get table and sample configs
        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)
        sampleset_config = dtablib.svpop.get_sampleset_config(table_def, config)

        # Determine if haplotype is in sample name
        sample_hap = table_def['sample_hap']

        if sample_hap is None:
            # Auto, determine if all samples end with a haplotype (-h1 or -h2)
            sample_hap = np.all([
                re.search(r'.*-.+$', sample) is not None for sample in sampleset_config['samples']
            ])

        # Read
        id_list = dtablib.dtabutil.get_id_list(input.id_list)

        # Read merge map
        df_map = pd.read_csv(input.merge_map, sep='\t', index_col='ID')

        # Read callable regions
        df_callable = pd.read_csv(input.callable, sep='\t', index_col='ID')

        # Check columns and variants
        missing_cols = set(df_map.columns) - set(df_callable.columns)

        if missing_cols:
            raise RuntimeError('Missing {} column(s) in callable regions: {}{}'.format(
                len(missing_cols), ', '.join(sorted(missing_cols)[:3]), '...' if len(missing_cols) > 3 else ''
            ))

        missing_id = set(id_list) - set(df_map.index)

        if missing_id:
            raise RuntimeError('Missing {} variant ID(s) in the merge map: {}{}'.format(
                len(missing_id), ', '.join(sorted(missing_id)[:3]), '...' if len(missing_id) > 3 else ''
            ))

        missing_id = set(id_list) - set(df_callable.index)

        if missing_id:
            raise RuntimeError('Missing {} variant ID(s) in the callable regions: {}{}'.format(
                len(missing_id), ', '.join(sorted(missing_id)[:3]), '...' if len(missing_id) > 3 else ''
            ))

        df_callable = df_callable.loc[id_list, df_map.columns]
        df_map = df_map.loc[id_list]

        # Transform map
        df_gt = (~ df_map.apply(pd.isnull)).astype(int)  # 1 or 0 genotype

        for sample in df_map.columns:
            # Change genotype to "." if not callable

            df_gt[sample] = pd.concat(
                [df_gt[sample], df_callable[sample]],
                axis=1
            ).apply(
                lambda row: row[0] if row[1] else '.',
                axis=1
            )

        # Check haplotypes
        hap_set = {val.rsplit('-', 1)[1] for val in df_gt.columns}

        # if len(hap_set - {'h1', 'h2'}) > 0:
        #     hap_list = sorted(hap_set - {'h1', 'h2'})
        #     n_hap = len(hap_list)
        #     hap_list_str = ', '.join(hap_list[:3])
        #     elipses = '...' if n_hap > 3 else ''
        #
        #     raise RuntimeError(f'Found {n_hap} unknown haplotypes: Only "h1" and "h2" are currently handled: {hap_list_str}{elipses}')

        hap_list = sorted(hap_set)

        # Get sample order
        sample_list = list()

        for sample_hap in df_gt.columns:
            sample = sample_hap.rsplit('-', 1)[0]
            if sample not in sample_list:
                sample_list.append(sample)

        # Create a table of per-sample genotypes
        df_gt_sample_list = list()

        for sample in sample_list:
            gt_col_list = list()

            for hap in hap_list:
                sample_hap = f'{sample}-{hap}'

                if sample_hap not in df_gt.columns:
                    gt_col_list.append(pd.Series(['.'] * df_gt.shape[0], index=df_gt.index))

                else:
                    gt_col_list.append(df_gt[sample_hap].astype(str))

            row_sample = pd.concat(gt_col_list, axis=1).apply(lambda row: '|'.join(row), axis=1)
            row_sample.name = sample

            df_gt_sample_list.append(row_sample)

        df_gt_sample = pd.concat(df_gt_sample_list, axis=1)

        # Write
        df_gt_sample.to_csv(output.tsv_sample, sep='\t', index=True, compression='gzip')
        df_gt.to_csv(output.tsv_hap, sep='\t', index=True, compression='gzip')

# # dtab_base_premerge_gt_by_svtype
# #
# # Make genotype table for source (pre-merged) sampleset variants.
# rule dtab_base_premerge_gt_by_svtype:
#     input:
#         bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
#         merge_map='sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}.tsv.gz',
#         tsv_sample=lambda wildcards: dtablib.svpop.sampleset_source_dict(
#             wildcards.tab_name, dtablib.svpop.VAR_PATTERN, config, wildcards.vartype, wildcards.svtype_single
#         ).values()
#     output:
#         tsv=temp('temp/{tab_name}/base_table/pre_merge/split/gt_{vartype}_{svtype}_{svtype_single}.tsv.gz')
#     wildcard_constraints:
#         svtype_single='ins|del|inv|snv|dup|sub|rgn'
#     run:
#
#         # Get ID list
#         df_base = pd.read_csv(input.bed_base, sep='\t', usecols=('ID', 'SVTYPE'))
#
#         id_list = list(df_base.loc[df_base['SVTYPE'] == wildcards.svtype_single.upper(), 'ID'])
#
#         # Get dict of input files
#         source_file_dict = dtablib.svpop.sampleset_source_dict(
#             wildcards.tab_name, dtablib.svpop.VAR_PATTERN, config, wildcards.vartype, wildcards.svtype_single
#         )
#
#         # Read merge map
#         df_merge = pd.read_csv(input.merge_map, sep='\t', index_col='ID')
#
#         df_merge = df_merge.loc[id_list]
#
#         # Define a function to translate a column in df_merge from variant IDs to genotypes
#         def to_gt(sample_col):
#
#             sample_col = sample_col.dropna()
#
#             # Get sample table
#             df_sample = pd.read_csv(source_file_dict[sample_col.name], sep='\t', index_col='ID')
#
#             if 'GT' not in df_sample.columns:
#                 raise RuntimeError('Missing GT column for sample "{}": {}'.format(
#                     sample_col.name,
#                     source_file_dict[sample_col.name]
#                 ))
#
#             # Check for variants missing that the merge set says should come from this sample
#             missing_set = set(sample_col) - set(df_sample.index)
#
#             if missing_set:
#                 raise RuntimeError('Missing {} variants for sample "{}": {}{}: {}'.format(
#                     len(missing_set),
#                     sample_col.name,
#                     ', '.join(sorted(missing_set)[:3]),
#                     '...' if len(missing_set) > 3 else '',
#                     source_file_dict[sample_col.name]
#                 ))
#
#             sample_gt = df_sample.loc[list(sample_col), 'GT']
#             sample_gt.index = list(sample_col.index)
#
#             return sample_gt
#
#         df_merge_gt = df_merge.apply(to_gt, axis=0)
#
#         df_merge_gt.index.name = 'ID'
#
#         # Write
#         df_merge_gt.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_base_callable
#
# Get a table of callable variants by variant ID (rows) and sample (columns). Values are True/False.
rule dtab_base_callable:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        bed_callable=lambda wildcards: dtablib.dtabutil.get_callable_bed_dict(wildcards.tab_name, config).values()
    output:
        tsv='sections/{tab_name}/base_table/callable/callable_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df = pd.read_csv(input.bed_base, sep='\t', usecols=('#CHROM', 'POS', 'END', 'ID', 'SVLEN'), index_col='ID')

        callable_dict = dtablib.dtabutil.get_callable_bed_dict(wildcards.tab_name, config)

        # Get table and sample configs
        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)
        sampleset_config = dtablib.svpop.get_sampleset_config(table_def, config)

        # Determine if haplotype is in sample name
        sample_hap = table_def['sample_hap']

        if sample_hap is None:
            # Auto, determine if all samples end with a haplotype (-h1 or -h2)
            sample_hap = np.all([
                '-' in sample for sample in sampleset_config['samples']
                #re.search(r'.*-h\d+', sample) is not None for sample in sampleset_config['samples']
            ])

        # Make callable table
        df_callable_list = list()

        for sample in callable_dict.keys():

            if sample_hap:
                sample_hap_list = [sample]
            else:
                sample_hap_list = ['-'.join([sample, hap]) for hap in ('h1', 'h2')]

            for sample_hap_str in sample_hap_list:
                df_sample = pd.read_csv(callable_dict[sample_hap_str], sep='\t')

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
            'temp/sections/{{tab_name}}/base_table/base_table_{{vartype}}_{{svtype}}/{part}.bed.gz',
            part=range(CHROM_PARTITIONS)
        ),
        id_list=expand(
            'temp/sections/{{tab_name}}/base_table/id_list_{{vartype}}_{{svtype}}/{part}.txt.gz',
            part=range(CHROM_PARTITIONS)
        ),
        sample_pkl=expand(
            'temp/sections/{{tab_name}}/base_table/sample_set_{{vartype}}_{{svtype}}/{part}.pkl',
            part=range(CHROM_PARTITIONS)
        ),
        merge_map=expand(
            'temp/sections/{{tab_name}}/base_table/merge_map_{{vartype}}_{{svtype}}/{part}.tsv.gz',
            part=range(CHROM_PARTITIONS)
        ),
        tsv_id_table=lambda wildcards: dtablib.dtabutil.get_table_def(wildcards.tab_name, config)['id_table'] if 'id_table' in dtablib.dtabutil.get_table_def(wildcards.tab_name, config) else []
    output:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        id_list='sections/{tab_name}/base_table/id_list_{vartype}_{svtype}.txt.gz',
        sample_pkl='sections/{tab_name}/base_table/sample_set_{vartype}_{svtype}.pkl',
        merge_map='sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        vartype='sv|indel|rgn|sub'
    run:
        dtablib.base_table.merge_base_table(input, output, wildcards, config)

# dtab_base_table_merge_table_snv
#
# Merge base tables.
rule dtab_base_table_merge_table_snv:
    input:
        bed_base=expand(
            'temp/sections/{{tab_name}}/base_table/base_table_{{vartype}}_{{svtype}}/{part}.bed.gz',
            part=range(CHROM_PARTITIONS)
        ),
        id_list=expand(
            'temp/sections/{{tab_name}}/base_table/id_list_{{vartype}}_{{svtype}}/{part}.txt.gz',
            part=range(CHROM_PARTITIONS)
        ),
        sample_pkl=expand(
            'temp/sections/{{tab_name}}/base_table/sample_set_{{vartype}}_{{svtype}}/{part}.pkl',
            part=range(CHROM_PARTITIONS)
        ),
        merge_map=expand(
            'temp/sections/{{tab_name}}/base_table/merge_map_{{vartype}}_{{svtype}}/{part}.tsv.gz',
            part=range(CHROM_PARTITIONS)
        ),
        tsv_id_table=lambda wildcards: dtablib.dtabutil.get_table_def(wildcards.tab_name, config)['id_table'] if 'id_table' in dtablib.dtabutil.get_table_def(wildcards.tab_name, config) else []
    output:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        id_list='sections/{tab_name}/base_table/id_list_{vartype}_{svtype}.txt.gz',
        sample_pkl='sections/{tab_name}/base_table/sample_set_{vartype}_{svtype}.pkl',
        merge_map='sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}.tsv.gz'
    wildcard_constraints:
        vartype='snv'
    run:
        dtablib.base_table.merge_base_table(input, output, wildcards, config)

# dtab_base_table_part
#
# Make base table annotation columns will be attached to.
rule dtab_base_table_part:
    input:
        bed=lambda wildcards: dtablib.base_table.get_input_bed(wildcards, config),
        tsv_ref_info='data/ref/contig_info.tsv.gz'
    output:
        bed_base=temp('temp/sections/{tab_name}/base_table/base_table_{vartype}_{svtype}/{part}.bed.gz'),
        id_list=temp('temp/sections/{tab_name}/base_table/id_list_{vartype}_{svtype}/{part}.txt.gz'),
        sample_pkl=temp('temp/sections/{tab_name}/base_table/sample_set_{vartype}_{svtype}/{part}.pkl'),
        merge_map=temp('temp/sections/{tab_name}/base_table/merge_map_{vartype}_{svtype}/{part}.tsv.gz')
    params:
        mem=_dtab_base_table_part_mem
    run:

        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

        sourcetype = table_def['sourcetype']

        # Get chromosomes in this partition
        df_ref = pd.read_csv(input.tsv_ref_info, sep='\t')

        df_ref = df_ref.loc[df_ref['PARTITION'] == int(wildcards.part)]

        chrom_set = set(df_ref['CHROM'])

        # Read
        df_list = list()

        df_iter = pd.read_csv(input.bed, sep='\t', low_memory=False, iterator=True, chunksize=20000)

        for df in df_iter:
            df = df.loc[df['#CHROM'].isin(chrom_set)]

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

# data_ref_contig_table
#
# Contig table.
rule base_reference_contig_table:
    output:
        tsv_ref_info='data/ref/contig_info.tsv.gz'
    params:
        partitions=CHROM_PARTITIONS
    run:

        # Get contig info
        df = svpoplib.ref.get_ref_info(
            config['reference']
        )

        # Add partition info
        partitions = svpoplib.partition.chrom.partition(df['LEN'], params.partitions)

        df['PARTITION'] = -1

        for index in range(len(partitions)):
            for chrom in partitions[index]:
                if df.loc[chrom, 'PARTITION'] != -1:
                    raise RuntimeError(f'Chromosome {chrom} assigned to multiple partitions')
                df.loc[chrom, 'PARTITION'] = index

        missing_list = list(df.loc[df['PARTITION'] == -1].index)

        if missing_list:
            raise RuntimeError(f'Failed assigning {len(missing_list)} chromosomes to partitions: {", ".join(missing_list)}')

        # Write
        df.to_csv(
            output.tsv_ref_info, sep='\t', index=True, compression='gzip'
        )
