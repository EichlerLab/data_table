"""
Write VCFs.
"""

# vcf_write_vcf
#
# Make VCF headers.
rule vcf_write_vcf:
    input:
        tsv='table/variants_{tab_name}_{vartype}_{svtype}.tsv.gz',
        tsv_gt='sections/{tab_name}/base_table/pre_merge/gt_sample_{vartype}_{svtype}.tsv.gz',
        ref_tsv='data/ref/contig_info.tsv.gz'
    output:
        vcf='vcf/variants_{tab_name}_{vartype}_{svtype}_{alt_fmt}.vcf.gz'
    wildcard_constraints:
        alt_fmt='alt|sym'
    run:

        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

        symbolic_alt = wildcards.alt_fmt == 'sym'

        if symbolic_alt and wildcards.svtype == 'snv':
            raise RuntimeError('Cannot generate symbolic ALT VCF for SNV')

        # Read
        df = pd.read_csv(input.tsv, sep='\t')
        df.set_index('ID', inplace=True, drop=False)
        df.index.name = 'INDEX'

        # Check

        if wildcards.svtype == 'snv':
            if 'REF' not in df.columns:
                raise RuntimeError('SNV VCF missing column: REF')

            if 'ALT' not in df.columns:
                raise RuntimeError('SNV VCF missing column: ALT')

            df['POS'] += 1

        # Support renaming '#CHROM' to 'CHROM'

        if 'CHROM' in df.columns and '#CHROM' not in df.columns:
            df['#CHROM'] = df['CHROM']


        ### Set format fields ###

        info_header_list = list()
        info_element = collections.defaultdict(list)

        # ID
        info_header_list.append(('ID', 'A', 'String', 'ID of merged variant set (max one variant per sample)'))

        for val in df['ID']:
            info_element[val].append(f'ID={val}')

        # VARTYPE (Needed for svpoplib.vcf.ref_base)
        if 'VARTYPE' not in df.columns:
            if wildcards.vartype not in {'sv', 'indel', 'snv'}:
                raise RuntimeError('VARTYPE field is missing and wildcards.vartype is not "sv", "indel", or "snv"')

            df['VARTYPE'] = wildcards.vartype.upper()

        info_header_list.append(('VARTYPE', 'A', 'String', 'Variant class'))

        for index, element in df['VARTYPE'].apply(lambda val: f'VARTYPE={val}').iteritems():
            info_element[index].append(element)

        # SVTYPE
        if 'SVTYPE' in df.columns:
            info_header_list.append(('SVTYPE', 'A', 'String', 'Variant type'))

            for index, element in df['SVTYPE'].apply(lambda val: f'SVTYPE={val}').iteritems():
                info_element[index].append(element)

        # SVLEN
        if 'SVLEN' in df.columns:
            info_header_list.append(('SVLEN', '.', 'Integer', 'Difference in length of ref and alt alleles'))

            for index, element in \
            (
                df.apply(lambda row: 'SVLEN={}'.format(
                    (-1 if row['SVTYPE'] == 'DEL' else 1) * np.abs(row['SVLEN'])
                ), axis=1)
            ).iteritems():
                info_element[index].append(element)

        # Lead sample
        if 'MERGE_SAMPLES' in df.columns:
            info_header_list.append(('SAMPLE', '.', 'String', 'Lead sample variant was called on. SV sequence, breakpoints, and contig locations come from the SV call in this sample.'))

            for index, element in df['MERGE_SAMPLES'].apply(lambda val: 'SAMPLE=' + val.split(',')[0]).iteritems():
                info_element[index].append(element)

        # TIG_REGION and QUERY_STRAND
        if 'TIG_REGION' in df.columns:
            info_header_list.append(('TIG_REGION', '.', 'String', 'Contig region of variant on an assembly (lead sample and haplotype only, may be present on other assemblies)'))

            for index, element in df['TIG_REGION'].apply(lambda val: 'TIG_REGION=' + val.split(';')[0]).iteritems():
                info_element[index].append(element)

        if 'QUERY_STRAND' in df.columns:
            info_header_list.append(('QUERY_STRAND', '.', 'String', 'Reference strand lead contig mapped to (+ or -)'))

            for index, element in df['QUERY_STRAND'].apply(lambda val: 'QUERY_STRAND=' + val.split(';')[0]).iteritems():
                info_element[index].append(element)

        # REF_SD and REF_TRF
        if 'REF_SD' in df.columns:
            info_header_list.append(('REF_SD', '.', 'Float', 'Max segmental duplication (SD) identity variant intersects.'))

            for index, element in df.loc[~ pd.isnull(df['REF_SD']), 'REF_SD'].iteritems():
                info_element[index].append(f'REF_SD={element:.2f}')

        if 'REF_TRF' in df.columns:
            info_header_list.append(('REF_TRF', '0', 'Flag', 'Variant intersects a reference TRF region'))

            for val in df.loc[df['REF_TRF']].index:
                info_element[val].append('REF_TRF')

        # HOMOP_LEN and DINUCL_LEN
        if 'HOMOP_LEN' in df.columns:
            info_header_list.append(('HOMOP_LEN', '.', 'Integer', 'Largest a reference homopolymer this variant intersects'))

            for index, element in df.loc[~ pd.isnull(df['HOMOP_LEN']), 'HOMOP_LEN'].iteritems():
                info_element[index].append(f'HOMOP_LEN={element:.2f}')

        if 'DINUCL_LEN' in df.columns:
            info_header_list.append(('DINUCL_LEN', '.', 'Integer', 'Largest a reference dinucleotide repeat this variant intersects'))

            for index, element in df.loc[~ pd.isnull(df['DINUCL_LEN']), 'DINUCL_LEN'].iteritems():
                info_element[index].append(f'DINUCL_LEN={element:.2f}')

        # Make INFO string
        df['INFO'] = pd.Series({
            index: ';'.join(vals) for index, vals in info_element.items()
        })

        ### REF ###
        if 'REF' not in df.columns:
            df['REF'] = list(svpoplib.vcf.ref_base(df, config['reference']))

        ### SEQ ###

        # Get SEQ field (Needed to set ALT if not symbolic, and needed as a SEQ info field for symbolic ALT)
        if 'SEQ' not in df.columns and ('ALT' not in df.columns or symbolic_alt):

            # Find FA file by first parsing in config elements, then wildcards
            seq_fa_file = os.path.join(
                config['svpop_dir'],
                'results/variant/{sourcetype}/{sourcename}/{sample}/{filter}/{svset}/bed/fa/{{vartype}}_{{svtype}}.fa.gz'.format(
                    **table_def
                ).format(**wildcards)
            )

            seq_fai_file = seq_fa_file + '.fai'

            if not os.path.isfile(seq_fa_file):
                raise RuntimeError(f'Sequnece FASTA file is missing: {seq_fa_file}')

            if not os.path.isfile(seq_fai_file):
                raise RuntimeError(f'Sequnece FASTA FAI file is missing: {seq_fai_file}')

            # Read
            id_set = set(df.index)

            with gzip.open(seq_fa_file, 'rt') as in_file:
                df_seq = pd.Series(
                    {record.id: str(record.seq).upper() for record in SeqIO.parse(in_file, 'fasta') if record.id in id_set}
                )

            missing_id_set = id_set - set(df_seq.index)

            if missing_id_set:
                raise RuntimeError('Missing {} IDs in FASTA: {}{}: {}'.format(
                    len(missing_id_set),
                    ', '.join(sorted(missing_id_set)[:3]),
                    '...' if len(missing_id_set) > 3 else '',
                    seq_fa_file
                ))

            df['SEQ'] = df_seq

            del df_seq

        ### ALT ###
        
        alt_header_list = list()

        if symbolic_alt:

            svtype_set = set(df['SVTYPE'])

            if 'INS' in svtype_set:
                alt_header_list.append(('INS', 'Sequence insertion'))
                svtype_set.remove('INS')

            if 'DEL' in svtype_set:
                alt_header_list.append(('DEL', 'Sequence deletion'))
                svtype_set.remove('DEL')

            if 'INV' in svtype_set:
                alt_header_list.append(('INV', 'Inversion'))
                svtype_set.remove('INV')

            if svtype_set:
                raise RuntimeError('Unknown SVTYPEs: '.format(', '.join(sorted(svtype_set))))

            df['ALT'] = df['SVTYPE'].apply(lambda val: f'<{val}>')

        elif wildcards.vartype != 'snv':

            # Check for sequence types that cannot be placed in ALT (may need symbolic ALTs)
            svtype_set_uk = set(df['SVTYPE']) - {'INS', 'DEL'}

            if svtype_set_uk:
                raise RuntimeError('Unknown or illegal SVTYPEs for sequence in ALT: {}'.format(','.join(sorted(svtype_set_uk))))

            df['REF'] = df.apply(lambda row: row['REF'] + row['SEQ'] if row['SVTYPE'] == 'DEL' else row['REF'], axis=1)
            df['ALT'] = df.apply(lambda row: row['REF'] + row['SEQ'] if row['SVTYPE'] == 'INS' else row['REF'][0], axis=1)

            del df['SEQ']


        ### QUAL, FILTER, FORMAT ###

        if 'QUAL' not in df.columns:
            df['QUAL'] = '.'

        if 'FILTER' in df.columns:
            raise RuntimeError('FILTER is defined in dataframe, but FILTER headers are not yet implemented')

        filter_header_list = list()

        df['FILTER'] = '.'

        df['FORMAT'] = 'GT'

        format_header_list = [
            ('GT', '1', 'String', 'Genotype')
        ]

        ### Rearrange to VCF order ###

        df = df[['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']]

        ### Add genotypes ###
        df_gt = pd.read_csv(input.tsv_gt, sep='\t', index_col='ID')

        gt_missing_id = set(df.index) - set(df_gt.index)

        if gt_missing_id:
            raise RuntimeError('GT dataframe missing {} IDs: {}{}: {}'.format(
                len(gt_missing_id),
                ', '.join(sorted(gt_missing_id)),
                '...' if len(gt_missing_id) > 3 else '',
                input.tsv_gt
            ))

        df = pd.concat([df, df_gt], axis=1)

        ### Write ###
        df_ref = pd.read_csv(input.ref_tsv, sep='\t')

        with Bio.bgzf.open(output.vcf, 'wt') as out_file:
            for line in svpoplib.vcf.header_list(
                df_ref,
                info_header_list,
                format_header_list,
                alt_header_list,
                filter_header_list,
                variant_source='PAV (HGSVC2)',
                ref_file_name=config.get('vcf_reference', config['reference'])
            ):
                out_file.write(line)

            out_file.write('\t'.join(df.columns))
            out_file.write('\n')

            for index, row in df.iterrows():
                out_file.write('\t'.join(row.astype(str)))
                out_file.write('\n')

        # Index
        shell("""tabix {output.vcf}""")
