"""
Generate browser tracks.
"""

# dtab_tracks_variant_bb
#
# BED to BigBed.
rule dtab_tracks_variant_bb:
    input:
        bed='temp/tracks/variant/{tab_name}/{vartype}_{svtype}.bed',
        asfile='temp/tracks/variant/{tab_name}/{vartype}_{svtype}.as',
        fai=config['reference'] + '.fai'
    output:
        bb='tracks/variant/{tab_name}/{tab_name}_{vartype}_{svtype}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""


# dtab_tracks_variant_bed
#
# Make track from variant BED.
rule dtab_tracks_variant_bed:
    input:
        tsv='table/variants_{tab_name}_{vartype}_{svtype}.tsv.gz'
    output:
        bed=temp('temp/tracks/variant/{tab_name}/{vartype}_{svtype}.bed'),
        asfile=temp('temp/tracks/variant/{tab_name}/{vartype}_{svtype}.as'),
        track_tsv=temp('temp/tracks/variant/{tab_name}/{vartype}_{svtype}.as.tsv')
    run:

        table_def = dtablib.dtabutil.get_table_def(wildcards.tab_name, config)

        # Read
        df = pd.read_csv(input.tsv, sep='\t', header=0)
        df_fai = svpoplib.ref.get_df_fai(config['reference'] + '.fai')

        # Set #CHROM
        if '#CHROM' not in df.columns and 'CHROM' in df.columns:
            df.columns = [col if col != 'CHROM' else '#CHROM' for col in df.columns]

        # Read default fields
        field_table_default = os.path.join(DATA_TABLE_DIR, 'files/tracks/variant_fields.tsv')

        if not os.path.isfile(field_table_default):
            raise RuntimeError(f'Cannot locate pipeline default field table: {field_table_default}: Missing or not a regular file')

        df_field_default = pd.read_csv(field_table_default, sep='\t')

        # Read custom field table
        field_table = table_def.get('track_tsv', None)

        if field_table is not None:

            if not os.path.isfile(field_table):
                raise RuntimeError(f'Cannot locate configured field table ("track_tsv" configuration option): {field_table}: Missing or not a regular file')

            if field_table.lower().endswith('.tsv') or field_table.lower().endswith('.tsv.gz'):
                df_field = pd.read_csv(field_table, sep='\t')
            elif field_table.lower().endswith('.xlsx') or field_table.lower().endswith('.xls'):
                df_field = pd.read_excel(field_table)
            else:
                raise RuntimeError(f'Unrecognized variant table file type (expected TSV or XLSX): {field_table}')

            # Add missing default fields
            set_missing = set(df_field_default['FIELD']) - set(df_field['FIELD'])

            if set_missing:
                df_field = pd.concat([
                    df_field,
                    df_field_default.loc[df_field_default['FIELD'].isin(set_missing)]
                ])

        else:
            df_field = df_field_default

        # Write field table
        df_field.to_csv(output.track_tsv, sep='\t', index=False)

        # Filter columns that have track annotations
        field_set = set(df_field['FIELD'])

        df = df.loc[:, [col for col in df.columns if col in field_set]]

        # Get track name and description for AS file
        track_name = 'VariantTable'
        track_description = '{tab_name} - {svtype}'.format(**wildcards)

        # Make track
        svpoplib.tracks.variant.make_bb_track(
            df, df_fai, output.bed, output.asfile,
            track_name, track_description,
            field_table_file_name=output.track_tsv
        )
