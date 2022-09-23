"""
Mapping-based annotations.
"""

MAPPING_MINIMAP2 = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/altmap/altmap-minimap2-distance-2-any-90_{vartype}_ins.bed.gz'  # INS only

#
# Rules
#

# rmsk_mapping_tandem_dup
#
# Tandem duplication by INS mapping to reference.
rule rmsk_mapping_tandem_dup:
    input:
        tab=lambda wildcards: dtablib.svpop.resolve_rel_path(MAPPING_MINIMAP2, wildcards, config)
    output:
        tsv='sections/{tab_name}/mapping/tandem-dup_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df_tandem = pd.read_csv(input.tab, sep='\t', header=0)

        # Shape
        df_tandem['POS1'] = df_tandem['POS'] + 1
        df_tandem['TANDEM_DUP'] = df_tandem.apply(lambda row: '{#CHROM}:{POS1}-{END}'.format(**row), axis=1)

        # Write
        df_tandem[['ID', 'TANDEM_DUP']].to_csv(output.tsv, sep='\t', index=False, compression='gzip')
