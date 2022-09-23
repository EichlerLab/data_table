"""
Get gene annotations (pLI, LOEUF) for gene hits.
"""

# dtab_geneset_pli
#
# Get pLI table for each variant.
rule dtab_geneset_pli:
    input:
        tsv=lambda wildcards: dtablib.dtabutil.section_input_files(wildcards, config, max_tier=1),
        pkl_list=lambda wildcards: dtablib.geneset.find_geneset_pkl(
            path='sections/{tab_name}/geneset'.format(**wildcards),
            path_re='.*_{vartype}_{svtype}.pkl'.format(**wildcards)
        )
    output:
        tsv='sections/{tab_name}/geneset/pli_{vartype}_{svtype}.tsv.gz'
    run:

        if len(input.pkl_list) == 0:
            raise RuntimeError('No genesets were output by previous rules')

        # Check for pLI table
        if 'pli_table' not in config:
            raise RuntimeError('Missing config element "pli_table" (Table of pLI values for each gene)')

        if not os.path.isfile(config['pli_table']):
            raise RuntimeError('Configured "pli_table" element points to a non-existent file: {}'.format(config['pli_table']))

        # Read
        df_pli_val = pd.read_csv(
            config['pli_table'],
            sep='\t',
            usecols=('GENE', 'PLI'), index_col='GENE',
            squeeze=True
        )

        # Read geneset tables. Each is a DataFrame of "ID" and "GENE" where "GENE" is a set of genes that ID affects.
        gene_dict = collections.defaultdict(set)

        for pkl_file in input.pkl_list:
            with open(pkl_file, 'rb') as in_file:
                df_gene_set = pickle.load(in_file)

            for id, gene_set in df_gene_set.items():
                gene_dict[id] |= gene_set

        df_var_gene = pd.Series(gene_dict)

        # Assign pLI
        df_pli_tuple = df_var_gene.apply(lambda gene_set:
            sorted([
                (
                    df_pli_val[gene] if gene in df_pli_val.index else -1000,
                    gene
                ) for gene in gene_set
            ], reverse=True)
        )

        # Get string of pLI values
        df_pli_str = df_pli_tuple.apply(lambda pli_list:
            ', '.join(
                [
                    '{}:{}'.format(pli_tuple[1], '{:.2f}'.format(pli_tuple[0]) if pli_tuple[0] >= 0 else 'NA')
                        for pli_tuple in pli_list
                ]
            )
        )

        df_pli_str.name = 'PLI_ALL'

        # Get max pLI
        df_pli_max = df_pli_tuple.apply(lambda pli_list:
            np.max([pli_tuple[0] for pli_tuple in pli_list])
        )

        df_pli_max = df_pli_max.apply(lambda val: val if val >= 0 else np.nan)
        df_pli_max = df_pli_max.apply(lambda val: '{:.2f}'.format(val) if not np.isnan(val) else '')

        df_pli_max.name = 'PLI_MAX'

        # Make dataframe and write
        df_pli = pd.concat([df_pli_max, df_pli_str], axis=1)
        df_pli.index.name = 'ID'

        # Write
        df_pli.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_geneset_loeuf
#
# Get LOEUF table for each variant.
rule dtab_geneset_loeuf:
    input:
        tsv=lambda wildcards: dtablib.dtabutil.section_input_files(wildcards, config, max_tier=1),
        pkl_list=lambda wildcards: dtablib.geneset.find_geneset_pkl(
            path='sections/{tab_name}/geneset'.format(**wildcards),
            path_re='.*_{vartype}_{svtype}.pkl'.format(**wildcards)
        )
    output:
        tsv='sections/{tab_name}/geneset/loeuf_{vartype}_{svtype}.tsv.gz'
    run:

        if len(input.pkl_list) == 0:
            raise RuntimeError('No genesets were output by previous rules')

        # Check for LOEUF table
        if 'loeuf_table' not in config:
            raise RuntimeError('Missing config element "loeuf_table" (Table of LOEUF values for each gene)')

        if not os.path.isfile(config['loeuf_table']):
            raise RuntimeError('Configured "loeuf_table" element points to a non-existent file: {}'.format(config['loeuf_table']))

        # Read
        df_loeuf_val = pd.read_csv(
            config['loeuf_table'],
            sep='\t'
        )

        df_loeuf_val = df_loeuf_val.loc[df_loeuf_val['CANONICAL']]

        df_loeuf_val = df_loeuf_val.groupby('GENE')['LOEUF'].apply(np.min)

        # Read geneset tables. Each is a DataFrame of "ID" and "GENE" where "GENE" is a set of genes that ID affects.
        gene_dict = collections.defaultdict(set)

        for pkl_file in input.pkl_list:
            with open(pkl_file, 'rb') as in_file:
                df_gene_set = pickle.load(in_file)

            for id, gene_set in df_gene_set.items():
                gene_dict[id] |= gene_set

        df_var_gene = pd.Series(gene_dict)

        # Assign LOEUF
        df_loeuf_tuple = df_var_gene.apply(lambda gene_set:
            sorted([
                (
                    df_loeuf_val[gene] if gene in df_loeuf_val.index else 10000,
                    gene
                ) for gene in gene_set
            ])
        )

        # Get string of LOEUF values
        df_loeuf_str = df_loeuf_tuple.apply(lambda loeuf_list:
            ', '.join(
                [
                    '{}:{}'.format(loeuf_tuple[1], '{:.2f}'.format(loeuf_tuple[0]) if loeuf_tuple[0] < 10000 else 'NA')
                        for loeuf_tuple in loeuf_list
                ]
            )
        )

        df_loeuf_str.name = 'LOEUF_ALL'

        # Get max LOEUF
        df_loeuf_min = df_loeuf_tuple.apply(lambda loeuf_list:
            np.min([loeuf_tuple[0] for loeuf_tuple in loeuf_list])
        )

        df_loeuf_min = df_loeuf_min.apply(lambda val: val if val < 10000 else np.nan)
        df_loeuf_min = df_loeuf_min.apply(lambda val: '{:.2f}'.format(val) if not np.isnan(val) else '')

        df_loeuf_min.name = 'LOEUF_MIN'

        # Make dataframe and write
        df_loeuf = pd.concat([df_loeuf_min, df_loeuf_str], axis=1)
        df_loeuf.index.name = 'ID'

        # Write
        df_loeuf.to_csv(output.tsv, sep='\t', index=True, compression='gzip')


# dtab_geneset_loeuf
#
# Get LOEUF table for each variant.
rule dtab_geneset_acmg59:
    input:
        tsv=lambda wildcards: dtablib.dtabutil.section_input_files(wildcards, config, max_tier=1),
        pkl_list=lambda wildcards:
        dtablib.geneset.find_geneset_pkl(
            path='sections/{tab_name}/geneset'.format(**wildcards),
            path_re='.*_{vartype}_{svtype}.pkl'.format(**wildcards)
        )
    output:
        tsv='sections/{tab_name}/geneset/acmg59_{vartype}_{svtype}.tsv.gz'
    run:

        if len(input.pkl_list) == 0:
            raise RuntimeError('No genesets were output by previous rules')

        # Check for LOEUF table
        if 'acmg59_table' not in config:
            raise RuntimeError('Missing config element "acmg59_table" (Table of ACMG genes)')

        if not os.path.isfile(config['acmg59_table']):
            raise RuntimeError('Configured "acmg59_table" element points to a non-existent file: {}'.format(config['acmg59_table']))

        # Read
        df_acmg = pd.read_csv(
            config['acmg59_table'],
            sep='\t'
        )

        df_acmg = df_acmg.groupby('Gene')['Disease'].apply(lambda vals: ', '.join(vals))

        # Read geneset tables. Each is a DataFrame of "ID" and "GENE" where "GENE" is a set of genes that ID affects.
        # Collet geneset information from all files into one gene_dict
        gene_dict = collections.defaultdict(set)

        for pkl_file in input.pkl_list:
            with open(pkl_file, 'rb') as in_file:
                df_gene_set = pickle.load(in_file)

            for id, gene_set in df_gene_set.items():
                gene_dict[id] |= gene_set

        df_var_gene = pd.Series(gene_dict)

        # Get disease list for each gene
        df_acmg_disease = df_var_gene.apply(lambda gene_set:
            ';'.join(sorted([
                (
                    '{} ({})'.format(df_acmg[gene], gene)
                ) for gene in gene_set if gene in df_acmg.index
            ])
        ))

        df_acmg_disease.name = 'ACMG59'
        df_acmg_disease.index.name = 'ID'

        # Write
        df_acmg_disease.to_csv(output.tsv, sep='\t', index=True, header=True, compression='gzip')
