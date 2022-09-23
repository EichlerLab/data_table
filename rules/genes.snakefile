"""
Annotations from intersecting variants with regions.
"""

REFSEQ_COUNT_PATTERN = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/refseq/refseq-count_{vartype}_{svtype}.tsv.gz'

REFSEQ_PROX_UP = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/refseq-prox/refseq-prox-up-{flank}_{vartype}_{svtype}.tsv.gz'
REFSEQ_PROX_DN = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/refseq-prox/refseq-prox-dn-{flank}_{vartype}_{svtype}.tsv.gz'


#
# Rules
#

# dtab_genes_refseq_intersect
#
# RefSeq intersect gene annotations.
rule dtab_genes_refseq_intersect:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(REFSEQ_COUNT_PATTERN, wildcards, config),
        id_list='sections/{tab_name}/base_table/id_list_{vartype}_{svtype}.txt.gz'
    output:
        tsv='sections/{tab_name}/genes/refseq-intersect_{vartype}_{svtype}.tsv.gz',
        geneset='sections/{tab_name}/geneset/genes-refseq-intersect_{vartype}_{svtype}.pkl'
    run:

        # Read
        df_refseq = pd.read_csv(input.tsv, sep='\t')

        id_list = dtablib.dtabutil.get_id_list(input.id_list)

        df_refseq_cds = df_refseq.loc[df_refseq['CDS'] > 0].groupby('ID')['GENE'].apply(lambda vals: ','.join(sorted(set(vals)))).reindex(id_list)
        df_refseq_cds.name = 'REFSEQ_CDS'

        df_refseq_ncrna = df_refseq.loc[df_refseq['NCRNA'] > 0].groupby('ID')['GENE'].apply(lambda vals: ','.join(sorted(set(vals)))).reindex(id_list)
        df_refseq_ncrna.name = 'REFSEQ_NCRNA'

        df_refseq_utr3 = df_refseq.loc[df_refseq['UTR3'] > 0].groupby('ID')['GENE'].apply(lambda vals: ','.join(sorted(set(vals)))).reindex(id_list)
        df_refseq_utr3.name = 'REFSEQ_UTR3'

        df_refseq_utr5 = df_refseq.loc[df_refseq['UTR5'] > 0].groupby('ID')['GENE'].apply(lambda vals: ','.join(sorted(set(vals)))).reindex(id_list)
        df_refseq_utr5.name = 'REFSEQ_UTR5'

        df_refseq_intron = df_refseq.loc[df_refseq['INTRON'] > 0].groupby('ID')['GENE'].apply(lambda vals: ','.join(sorted(set(vals)))).reindex(id_list)
        df_refseq_intron.name = 'REFSEQ_INTRON'

        # Get RefSeq hits (gene names keyed by variant ID)
        df = pd.concat(
            [
                df_refseq_cds,
                df_refseq_ncrna,
                df_refseq_utr3,
                df_refseq_utr5,
                df_refseq_intron
            ],
            axis=1
        )

        # Get gene set
        df_gene_set = df_refseq.groupby('ID')['GENE'].apply(set)

        with open(output.geneset, 'wb') as out_file:
            pickle.dump(df_gene_set, out_file)

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_genes_refseq_prox
#
# Refseq proximal hits.
rule dtab_genes_refseq_prox:
    input:
        tsv_up=lambda wildcards: dtablib.svpop.resolve_rel_path(REFSEQ_PROX_UP, wildcards, config),
        tsv_dn=lambda wildcards: dtablib.svpop.resolve_rel_path(REFSEQ_PROX_DN, wildcards, config),
        id_list='sections/{tab_name}/base_table/id_list_{vartype}_{svtype}.txt.gz'
    output:
        tsv='sections/{tab_name}/genes/refseq-prox-{flank}_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df_prox_up = pd.read_csv(input.tsv_up, sep='\t', header=0, usecols=('ID', 'GENE', 'DISTANCE'))
        df_prox_dn = pd.read_csv(input.tsv_dn, sep='\t', header=0, usecols=('ID', 'GENE', 'DISTANCE'))

        id_list = dtablib.dtabutil.get_id_list(input.id_list)

        # Shape
        df_prox_up['GENE_DIST'] = df_prox_up.apply(lambda row: '{GENE}:{DISTANCE}'.format(**row), axis=1)
        df_prox_dn['GENE_DIST'] = df_prox_dn.apply(lambda row: '{GENE}:{DISTANCE}'.format(**row), axis=1)

        # Merge all intersections into a list with each element GENE:DISTANCE, key by SV ID

        prox_up = df_prox_up.groupby('ID')['GENE_DIST'].apply(lambda val: ','.join(sorted(set(val)))).reindex(id_list)
        prox_up.name = f'REFSEQ_UP_{wildcards.flank}'

        prox_dn = df_prox_dn.groupby('ID')['GENE_DIST'].apply(lambda val: ','.join(sorted(set(val)))).reindex(id_list)
        prox_dn.name = f'REFSEQ_DN_{wildcards.flank}'

        df = pd.concat([prox_up, prox_dn], axis=1)

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')
