"""
Inter-sequence annotations.
"""

INTERSEQ_TRF = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/trf/trf-table_{vartype}_{svtype}.tsv.gz'
INTERSEQ_RMSK = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/rmsk/rmsk-table_{vartype}_{svtype}.tsv.gz'
INTERSEQ_GC = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/gc/gc_content_{vartype}_{svtype}.tsv.gz'


#
# Rules
#

# dtab_repeats_trf
#
# TRF run on variant sequence.
rule dtab_repeats_trf:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(INTERSEQ_TRF, wildcards, config)
    output:
        tsv='sections/{tab_name}/interseq/trf-{params}_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df_svlen = pd.read_csv(input.bed_base, sep='\t', usecols=('ID', 'SVLEN'), index_col='ID', squeeze=True)
        df_trf = pd.read_csv(input.tsv, sep='\t', header=0)

        # Copies
        trf_copies = df_trf.groupby('ID').apply(lambda subdf:
            ', '.join(sorted(set(['{PERIOD} ({COPIES}x)'.format(**row) for index, row in subdf.iterrows()])))
        )

        trf_copies.name = 'TRF_COPIES'

        # Longest perfect repeat
        perfect_n = df_trf.groupby('ID')['LONGEST_PERFECT_N'].apply(max)
        perfect_n.name = 'TRF_PFCT_N'

        perfect_bp = df_trf.groupby('ID')['LONGEST_PERFECT_BP'].apply(max).fillna(0).astype(np.int32)
        perfect_bp.name = 'TRF_PFCT_BP'

        # Percent covered
        df_trf['TRF_LEN'] = df_trf['SV_END'] - df_trf['SV_POS']
        prop_covered = (df_trf.groupby('ID')['TRF_LEN'].max() / df_svlen).dropna()
        prop_covered.name = 'TRF_PROP'

        # Concat and write
        df = pd.concat(
            [trf_copies, perfect_n, perfect_bp, prop_covered],
            axis=1
        )

        df.index.name = 'ID'

        df.to_csv(
            output.tsv, sep='\t', index=True, compression='gzip'
        )

# dtab_repeats_rmsk
#
# RepeatMasker run on variant sequence.
rule dtab_repeats_rmsk:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(INTERSEQ_RMSK, wildcards, config)
    output:
        tsv='sections/{tab_name}/interseq/rmsk_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df_rmsk = pd.read_csv(input.tsv, sep='\t', header=0, usecols=('ID', 'REPEAT_FAMILY', 'REPEAT_TYPE'))

        df_rmsk = df_rmsk.loc[df_rmsk['REPEAT_FAMILY'] != 'Simple_repeat']
        df_rmsk = df_rmsk.loc[df_rmsk['REPEAT_FAMILY'] != 'Low_complexity']

        # Shape - Join REPEAT_TYPE for each ID/REPEAT_FAMILY
        df_rmsk = df_rmsk.groupby(['ID', 'REPEAT_FAMILY']).apply(
            lambda subdf: '{} ({})'.format(
                subdf['REPEAT_FAMILY'].iloc[0],
                ', '.join(sorted(set(subdf['REPEAT_TYPE'])))
            )
        )

        df_rmsk.index = df_rmsk.index.droplevel('REPEAT_FAMILY')

        # Shape - Join REPEAT_TYPE/FAMILY per ID
        df_rmsk = df_rmsk.groupby('ID').apply(lambda subdf:
            ', '.join(sorted(set(subdf)))
        )

        df_rmsk.name = 'RMSK'

        # Add
        df_rmsk.to_csv(output.tsv, sep='\t', index=True, header=True, compression='gzip')

# dtab_interseq_gc
#
# Inter-sequence GC.
rule dtab_interseq_gc:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(INTERSEQ_GC, wildcards, config)
    output:
        tsv='sections/{tab_name}/interseq/gc_{vartype}_{svtype}.tsv.gz'
    run:

        pd.read_csv(
            input.tsv, sep='\t', header=0, usecols=('ID', 'GC_CONTENT'), index_col='ID', squeeze=False
        ).to_csv(
            output.tsv, sep='\t', index=True, header=True, compression='gzip'
        )
