"""
Reference region annotations.
"""

REGIONS_SD = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/sd/sd-max-frac_{vartype}_{svtype}.tsv.gz'
REGIONS_TRF = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/trf/trf_regions_{params}_{vartype}_{svtype}.tsv.gz'
REGIONS_ENCODE = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/allreg/allreg_{vartype}_{svtype}.tsv.gz'
REGIONS_DHS2020 = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/dhs/dhs2020-1000_{vartype}_{svtype}.tsv.gz'
REGIONS_CCRE2020 = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/ccre/ccre2020-all_{vartype}_{svtype}.tsv.gz'
REGIONS_OREGANNO = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/oreganno/oreganno_{vartype}_{svtype}.tsv.gz'
REGIONS_HOMOP = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/homopolymer/homopolymer_intersect_all_any_{vartype}_{svtype}.tsv.gz'
REGIONS_HOMOP_UP = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/homopolymer/homopolymer_nearest_all_up_{vartype}_{svtype}.tsv.gz'
REGIONS_HOMOP_DN = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/homopolymer/homopolymer_nearest_all_dn_{vartype}_{svtype}.tsv.gz'

REGIONS_DINUCL = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/homopolymer/dinucleotide_intersect_all_any_{vartype}_{svtype}.tsv.gz'
REGIONS_DINUCL_UP = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/homopolymer/dinucleotide_nearest_all_up_{vartype}_{svtype}.tsv.gz'
REGIONS_DINUCL_DN = 'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/anno/homopolymer/dinucleotide_nearest_all_dn_{vartype}_{svtype}.tsv.gz'

# dtab_regions_win
#
# Get a region string with a flanking windown around each variant.
rule dtab_regions_win:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz'
    output:
        tsv='sections/{tab_name}/regions/ref_win_{win_spec}_{vartype}_{svtype}.tsv.gz'
    run:

        # Get win flank
        win_spec = wildcards.win_spec.lower()

        win_mult = 1

        if win_spec.endswith('k'):
            win_mult = 1000
            win_spec = win_spec[:-1]

        elif win_spec.endswith('m'):
            win_mult = 1000000
            win_spec = win_spec[:-1]

        try:
            win_size = np.int64(win_spec) * win_mult

        except ValueError:
            raise RuntimeError('Window spec is not an integer: {win_spec}'.format(**wildcards))

        # Get column name
        col_name = 'WIN_{}'.format(wildcards.win_spec.upper())

        # Read
        df = pd.read_csv(input.bed_base, sep='\t', usecols=('#CHROM', 'POS', 'END', 'ID'))
        df_fai = svpoplib.ref.get_df_fai(REF_FAI)

        fai_dict = dict(df_fai)

        # Get window
        df['POS'] = df['POS'].astype(np.int64)
        df['END'] = df['END'].astype(np.int64)

        df['POS'] = df['POS'] - win_size + 1
        df.loc[df['POS'] < 1, 'POS'] = 0

        df['END'] += win_size
        df['END'] = df.apply(lambda row: np.min([row['END'], fai_dict.get(row['#CHROM'], row['END'])]), axis=1)

        #df['END'] = df.apply(lambda row: nrow['END'] if row['#CHROM'] not in df_fai.index or row['END'] > df_fai[row['#CHROM']] else df_fai[row['#CHROM']], axis=1)

        df[col_name] = df.apply(lambda row: '{#CHROM}:{POS}-{END}'.format(**row), axis=1)

        # Write
        df[['ID', col_name]].to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# dtab_regions_homop
#
# Homopolymer intersections.
rule dtab_regions_homop:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_HOMOP, wildcards, config)
    output:
        tsv='sections/{tab_name}/regions/homop-intersect_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df = pd.read_csv(input.tsv, sep='\t')

        # Drop duplicates
        df.sort_values('LEN', inplace=True)
        df.drop_duplicates('ID', keep='first', inplace=True)

        # Shape
        df = df[['ID', 'LEN', 'BASE']]
        df.columns = ['ID', 'HOMOP_LEN', 'HOMOP_BASE']

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

rule dtab_regions_homop_prox:
    input:
        tsv_up=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_HOMOP_UP, wildcards, config),
        tsv_dn=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_HOMOP_DN, wildcards, config)
    output:
        tsv='sections/{tab_name}/regions/homop-nearest_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df_up = pd.read_csv(input.tsv_up, sep='\t', index_col='ID')
        df_dn = pd.read_csv(input.tsv_dn, sep='\t', index_col='ID')

        # Shape
        df_up['HOMOP_UP_BASE'] = df_up.apply(lambda row: '{BASE}x{LEN:d}'.format(**row), axis=1)
        df_dn['HOMOP_DN_BASE'] = df_dn.apply(lambda row: '{BASE}x{LEN:d}'.format(**row), axis=1)

        df_up = df_up[['DISTANCE', 'HOMOP_UP_BASE']]
        df_up.columns = ['HOMOP_UP_DIST', 'HOMOP_UP_BASE']

        df_dn = df_dn[['DISTANCE', 'HOMOP_DN_BASE']]
        df_dn.columns = ['HOMOP_DN_DIST', 'HOMOP_DN_BASE']

        # Merge
        df = pd.concat([df_up, df_dn], axis=1)

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_regions_dinucl
#
# Dinucleotide intersections.
rule dtab_regions_dinucl:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_DINUCL, wildcards, config)
    output:
        tsv='sections/{tab_name}/regions/dinucl-intersect_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df = pd.read_csv(input.tsv, sep='\t')

        # Drop duplicates
        df.sort_values('LEN', inplace=True)
        df.drop_duplicates('ID', keep='first', inplace=True)

        # Shape
        df = df[['ID', 'LEN', 'BASE']]
        df.columns = ['ID', 'DINUCL_LEN', 'DINUCL_BASE']

        # Write
        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# dtab_regions_dinucl_prox
#
# Proximal dinucleotide repeat.
rule dtab_regions_dinucl_prox:
    input:
        tsv_up=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_DINUCL_UP, wildcards, config),
        tsv_dn=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_DINUCL_DN, wildcards, config)
    output:
        tsv='sections/{tab_name}/regions/dinucl-nearest_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df_up = pd.read_csv(input.tsv_up, sep='\t', index_col='ID')
        df_dn = pd.read_csv(input.tsv_dn, sep='\t', index_col='ID')

        # Shape
        df_up['DINUCL_UP_BASE'] = df_up.apply(lambda row: '{BASE}x{LEN:d}'.format(**row), axis=1)
        df_dn['DINUCL_DN_BASE'] = df_dn.apply(lambda row: '{BASE}x{LEN:d}'.format(**row), axis=1)

        df_up = df_up[['DISTANCE', 'DINUCL_UP_BASE']]
        df_up.columns = ['DINUCL_UP_DIST', 'DINUCL_UP_BASE']

        df_dn = df_dn[['DISTANCE', 'DINUCL_DN_BASE']]
        df_dn.columns = ['DINUCL_DN_DIST', 'DINUCL_DN_BASE']

        # Merge
        df = pd.concat([df_up, df_dn], axis=1)

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_regions_sd
#
# SD region intersect table.
rule dtab_regions_sd:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_SD, wildcards, config)
    output:
        tsv='sections/{tab_name}/regions/sd_{vartype}_{svtype}.tsv.gz'
    run:

        df = pd.read_csv(input.tsv, sep='\t')

        df = df.loc[df['MATCH_MAX'] != '.', ['ID', 'MATCH_MAX']]
        df.columns = ['ID', 'REF_SD']

        df.to_csv(output.tsv, sep='\t', index=False, compression='gzip')

# dtab_regions_trf
#
# TRF region intersect table.
rule dtab_regions_trf:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(
            REGIONS_TRF, wildcards, config,
            param_sub=lambda param: f'200_0_{param}' if param == 'any' or re.search(r'^\d\d?', param) else param
        ),
        id_list='sections/{tab_name}/base_table/id_list_{vartype}_{svtype}.txt.gz'
    output:
        tsv='sections/{tab_name}/regions/trf_{vartype}_{svtype}.tsv.gz'
    run:
        df = pd.read_csv(input.tsv, sep='\t', squeeze=False)
        id_list = dtablib.dtabutil.get_id_list(input.id_list)

        df['REF_TRF'] = True
        df.set_index('ID', inplace=True)
        df = df.reindex(id_list, fill_value=False)

        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_regions_ccre2020
#
# Candidate cis-regulatory (cCRE) regions (ENCODE 2020)
rule dtab_regions_ccre2020:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_CCRE2020, wildcards, config)
    output:
        tsv='sections/{tab_name}/regions/ccre2020_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df = pd.read_csv(
            input.tsv,
            sep='\t',
            dtype=str
        )

        df = df.groupby('ID').apply(lambda subdf: pd.Series(
            [
                ','.join(sorted(
                    {
                        val for val_list in subdf['CCRE2020'] for val in (
                            val_list.split(',') if not pd.isnull(val_list) else {}
                        )
                    }
                )),
                ','.join(['{EH38D}:{EH38E}'.format(**row) for index, row in subdf.iterrows()])
            ],
            index=['CCRE2020', 'CCRE2020_ACC']
        ))

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_regions_dhs2020
#
# Regulatory regions (Vierstra 2020)
rule dtab_regions_dhs2020:
    input:
        bed_base='sections/{tab_name}/base_table/base_table_{vartype}_{svtype}.bed.gz',
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_DHS2020, wildcards, config)
    output:
        tsv='sections/{tab_name}/regions/dhs2020_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df = pd.read_csv(
            input.tsv,
            sep='\t',
            usecols=('ID', 'DHS_ID'),
            dtype=str
        )

        df = df.groupby('ID')['DHS_ID'].apply(lambda vals: ';'.join(vals))
        df.name = 'DHS2020'

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, header=True, compression='gzip')

# dtab_regions_encode
#
# Regulatory regions.
rule dtab_regions_encode:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_ENCODE, wildcards, config)
    output:
        tsv='sections/{tab_name}/regions/encode_{vartype}_{svtype}.tsv.gz'
    run:

        # Read
        df_all_reg = pd.read_csv(
            input.tsv,
            sep='\t',
            usecols=('ID', 'CPG_ISLAND', 'ENCODE_H3K27AC', 'ENCODE_H3K4ME3', 'ENCODE_H3K4ME1', 'DHS'),
            index_col='ID'
        )

        # Shape
        set_cpg = set(df_all_reg.loc[df_all_reg['CPG_ISLAND']].index)

        reg_col_names = ('H3K27Ac', 'H3K4Me3', 'H3K4Me1', 'DHS')
        df_all_reg = df_all_reg.loc[:, ('ENCODE_H3K27AC', 'ENCODE_H3K4ME3', 'ENCODE_H3K4ME1', 'DHS')]
        df_all_reg.columns = reg_col_names

        df_all_reg = df_all_reg.loc[df_all_reg.apply(any, axis=1)].copy()

        df = pd.DataFrame(df_all_reg.apply(lambda row: ','.join(row.loc[row].index), axis=1))
        df.columns = ['ENCODE']

        # Set CpG
        df['CPG'] = [val in set_cpg for val in df.index]

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

# dtab_regions_oreganno
#
# Annotate ORegAnno hits.
rule dtab_regions_oreganno:
    input:
        tsv=lambda wildcards: dtablib.svpop.resolve_rel_path(REGIONS_OREGANNO, wildcards, config),
        tsv_oreg=config['oreganno_table'] if 'oreganno_table' in config else [],
        id_list='sections/{tab_name}/base_table/id_list_{vartype}_{svtype}.txt.gz'
    output:
        tsv='sections/{tab_name}/regions/oreganno_{vartype}_{svtype}.tsv.gz',
        geneset='sections/{tab_name}/geneset/regions-oreganno_{vartype}_{svtype}.pkl'
    run:

        # Check for ORegAnno table
        if 'oreganno_table' not in config:
            raise RuntimeError('Missing config element "oreganno_table" (Table of ORegAnno IDs and their annotations)')

        if not os.path.isfile(config['oreganno_table']):
            raise RuntimeError('Configured "oreganno_table" element points to a non-existent file: {}'.format(config['oreganno_table']))

        # Read
        id_list = dtablib.dtabutil.get_id_list(input.id_list)

        df_oreg = pd.read_csv(input.tsv, sep='\t', header=0)
        df_oreg = df_oreg.loc[~ pd.isnull(df_oreg['OREG_ID'])]

        oreg_id_set = set(df_oreg['OREG_ID'].dropna())

        df_oreg_table = pd.read_csv(
            config['oreganno_table'],
            sep='\t', header=0,
            usecols=('ORegAnno_ID', 'Type', 'Gene_Symbol', 'Regulatory_Element_Symbol'),
            low_memory=False
        )

        df_oreg_table = df_oreg_table.loc[[val in oreg_id_set for val in df_oreg_table['ORegAnno_ID']]]

        # Normalize ORegAnno data table
        df_oreg_table['Type'] = df_oreg_table['Type'].apply(lambda val: val.upper())

        oreg_type_trans = {
            'REGULATORY POLYMORPHISM': 'REG_POLY',
            'REGULATORY HAPLOTYPE': 'REG_HAPLO',
            'TRANSCRIPTION FACTOR BINDING SITE': 'TFBS',
            'REGULATORY REGION': 'REG_REGION',
            'OPERON': 'OPERON',
            'SRNA BINDING SITE': 'sRNA_BS',
            'MIRNA BINDING SITE': 'miRNA_BS',
        }

        df_oreg_table['Type'] = df_oreg_table['Type'].apply(lambda val: oreg_type_trans[val])

        df_oreg_table['Gene_Symbol'].fillna('Unknown', inplace=True)
        df_oreg_table['Regulatory_Element_Symbol'].fillna('Unknown', inplace=True)

        # Aggregate by ORegAnno ID
        df_oreg_gene = df_oreg_table.groupby('ORegAnno_ID')['Gene_Symbol'].apply(lambda vals: set(vals))
        df_oreg_element = df_oreg_table.groupby('ORegAnno_ID')['Regulatory_Element_Symbol'].apply(lambda vals: set(vals))

        # Assign genes and elements to df_oreg
        df_oreg['GENE'] = list(df_oreg_gene.loc[list(df_oreg['OREG_ID'])])
        df_oreg['REG_ELEMENT'] = list(df_oreg_element.loc[list(df_oreg['OREG_ID'])])

        # Aggregate df_oreg by ID
        oreg_gene = df_oreg.groupby('ID')['GENE'].apply(lambda valsetlist: ','.join(sorted({val for valset in valsetlist for val in valset})))
        oreg_gene.name = 'OREG_GENE'

        oreg_regelement = df_oreg.groupby('ID')['REG_ELEMENT'].apply(lambda valsetlist: ','.join(sorted({val for valset in valsetlist for val in valset})))
        oreg_regelement.name = 'OREG_ELEMENT'

        oreg_ctcf = df_oreg.groupby('ID')['REG_ELEMENT'].apply(lambda valsetlist: np.any(list({val == 'CTCF' for valset in valsetlist for val in valset})))
        oreg_ctcf.name = 'OREG_CTCF'

        # Merge
        df = pd.concat([oreg_gene, oreg_regelement, oreg_ctcf], axis=1).reindex(id_list)

        df['OREG_CTCF'] = df['OREG_CTCF'].fillna(False).astype(bool)

        # Write
        df.to_csv(output.tsv, sep='\t', index=True, compression='gzip')

        # Create gene set table
        df_gene_set = oreg_gene.apply(lambda gene_list: {val for val in gene_list.split(',') if val != 'Unknown'})
        df_gene_set = df_gene_set.loc[df_gene_set.apply(len) > 0]

        with open(output.geneset, 'wb') as out_file:
            pickle.dump(df_gene_set, out_file)
