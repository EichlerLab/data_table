"""
SNV-specific annotations.
"""

# dtab_snv_transitions
#
# Annotate transition SNVs
rule dtab_snv_transitions:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_snv_snv.bed.gz'
    output:
        tsv='sections/{tab_name}/snv/transitions_snv_snv.tsv.gz'
    run:

        TI_SET = {
            ('A', 'C'): False,
            ('A', 'G'): True,
            ('A', 'T'): False,

            ('C', 'A'): False,
            ('C', 'G'): False,
            ('C', 'T'): True,

            ('G', 'A'): True,
            ('G', 'C'): False,
            ('G', 'T'): False,

            ('T', 'A'): False,
            ('T', 'C'): True,
            ('T', 'G'): False,
        }

        df = pd.read_csv(input.bed_base, sep='\t', usecols=('ID', 'REF', 'ALT'), index_col='ID')

        df['REF'] = df['REF'].apply(str.upper)
        df['ALT'] = df['ALT'].apply(str.upper)

        df['SNV_TRANS'] = df.apply(lambda row: TI_SET.get((row['REF'], row['ALT']), np.nan), axis=1)

        del(df['REF'])
        del(df['ALT'])

        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')
