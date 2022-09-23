"""
Generate browser tracks.
"""

# dtab_tracks_variant_bb
#
# BED to BigBed.
rule dtab_tracks_variant_bb:
    input:
        bed='temp/tracks/variant/{tab_name}_{vartype}_{svtype}.bed',
        asfile='temp/tracks/variant/{tab_name}_{vartype}_{svtype}.as',
        fai=config['reference'] + '.fai'
    output:
        bb='tracks/variant/{tab_name}_{vartype}_{svtype}.bb'
    shell:
        """bedToBigBed -tab -as={input.asfile} -type=bed9+ {input.bed} {input.fai} {output.bb}"""


# dtab_tracks_variant_bed
#
# Make track from variant BED.
rule dtab_tracks_variant_bed:
    input:
        tsv='tsv/variants_{tab_name}_{vartype}_{svtype}.tsv.gz'
    output:
        bed=temp('temp/tracks/variant/{tab_name}_{vartype}_{svtype}.bed'),
        asfile=temp('temp/tracks/variant/{tab_name}_{vartype}_{svtype}.as')
    run:

        # Read
        df = pd.read_csv(input.tsv, sep='\t', header=0)
        df_fai = svpoplib.ref.get_df_fai(config['reference'] + '.fai')

        field_table_file_name = os.path.join(DATA_TABLE_DIR, 'files/tracks/variant_fields.tsv')

        # Filter columns that have track annotations
        field_set = set(
            pd.read_csv(
                field_table_file_name,
                sep='\t', header=0
            )['FIELD']
        )

        df = df.loc[:, [col for col in df.columns if col in field_set]]

        # Get track name and description for AS file
        track_name = 'VariantTable'
        track_description = '{tab_name} - {svtype}'.format(**wildcards)

        # Make track
        svpoplib.tracks.variant.make_bb_track(
            df, df_fai, output.bed, output.asfile,
            track_name, track_description,
            field_table_file_name=field_table_file_name
        )
