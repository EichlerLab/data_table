# Callset support configuration #

In the final variant table (one row per variant), several columns can be added to annotate support from other
sources including pre-formatted tables or SV-Pop intersections. For merged callsets, the pipeline will pull
support annotations from individual samples.

## Configuration section ##

The configuration section is called "support". It is a dictionary with keys for each variant type (vartype_svtype).
The next level is the name of the support column in the final table, and it contains elements telling this pipeline
where to find files and how to process them.

For each support section, `path` and `type` are required attributes.

* path: Path to input files with wildcards for sample and variant types.
* type: Type of support.

### Example ###
```
{
    "svpop_dir": "/path/to/svpop",
    "reference": "/path/to/ref.fa.gz",
    "table_def": {
        "section": {
            "sourcetype": "sampleset",
            "sourcename": "pav",
            "sample": "batch1:hap",
            "filter": "all",
            "svset": "all",
            "callable_filter": "/path/to/pav/results/{sample}/callable/callable_regions_{hap}_500.bed.gz",
            "support": {
                "sv_insdel": {
                    "PAVLRA": {
                        "path": "/path/to/svpop/results/variant/intersect/caller+pav+{sample}/caller+pav-lra+{sample}/szro-50-200/all/all/{vartype}_{svtype}/intersect.tsv.gz",
                        "type": "svpopinter"
                    },
                    "PBSV": {
                        "path": "/path/to/svpop/results/variant/intersect/caller+pav+{sample_l}/caller+pbsv+{sample_r}/szro-50-200/all/all/{vartype}_{svtype}/intersect.tsv.gz",
                        "type": "svpopinter-striphap"
                    }
                }
            }
        }
    }
}
```

This defines two support sources, PAVLRA and PBSV.

### Using support_sections ###

If the same support section appears multiple times (e.g. inside "sv_insdel" and "indel_insdel" support sections), the
definition can be separated into a section called "support_sections" and referenced by name in the "support" section.

Example:
```
{
    "svpop_dir": "/path/to/svpop",
    "reference": "/path/to/ref.fa.gz",
    "table_def": {
        "section": {
            "sourcetype": "sampleset",
            "sourcename": "pav",
            "sample": "batch1:hap",
            "filter": "all",
            "svset": "all",
            "callable_filter": "/path/to/pav/results/{sample}/callable/callable_regions_{hap}_500.bed.gz",
            "support": {
                "sv_insdel": {
                    "PAVLRA": "PAVLRA",
                    "PBSV": "PBSV"
                },
                "indel_insdel": {
                    "PAVLRA": "PAVLRA",
                    "PBSV": "PBSV"
                }
            },
            "support_sections": {
                "PAVLRA": {
                    "path": "/path/to/svpop/results/variant/intersect/caller+pav+{sample}/caller+pav-lra+{sample}/szro-50-200/all/all/{vartype}_{svtype}/intersect.tsv.gz",
                    "type": "svpopinter"
                },
                "PBSV": {
                    "path": "/path/to/svpop/results/variant/intersect/caller+pav+{sample_l}/caller+pbsv+{sample_r}/szro-50-200/all/all/{vartype}_{svtype}/intersect.tsv.gz",
                    "type": "svpopinter-striphap"
                }
            }
        }
    }
}
```

The support section names in "support" have two parts. The first is the name of the section (before the colon) and the
second is the name of the element in "support_sections" to pull the configuration from (after the colon). In most cases,
these will be the same, but it allows for some flexibility.


### Supported types ###

* svpopinter: Pull from an SV-Pop intersects where sample names match on both sides.
  * Path: Path to the SV-Pop "intersect.tsv.gz" file. 
  * Path wildcards: sample, vartype, svtype
* svpopinter-striphap: Pull from an SV-Pop where the left sample (A) has haplotype appended ("-h1" or "-h2"), but the
  right sample (B) does not.
  * Path: Path to the SV-Pop "intersect.tsv.gz" file.
  * Path wildcards: sample_l (left sample name), sample_r (right sample name), vartype, svtype 
* subseq: Read output from the subseq validation pipeline
* table: Pull annotations from individual samples using the lead sample in the merge.
* table_bool: Same as table, but expect a table of IDs and write True/False values if the variant ID is present
  or absent (respectively).
* preformat: Preformatted table using the IDs in the merge.


### Other options ###

* allow-missing: Allow missing samples. Useful if intersects were not done for all samples, for example, checking
  against another callset where all samples were not present. The value of this option should be "true" or "false"
  (not case sensitive, but other values such as "t", "f", "0", "1", "yes", and "no" are recognized).
* column-name: When pulling from a table, all columns are extracted by default. This option specifies the column name
  to keep and discards others.
* sample-pattern: A regular expression that will match valid samples. For example, `"sample-pattern": ".*-(h1|h2)$"`
  would restrict the support section for sample names ending with "h1" or "h2", but not "unassigned". All samples not
  matching the pattern are ignored and left empty. Sample pattern `"sample-pattern": "^sample1|sample2|sample3$"` would
  restrict analysis to just three samples, "sample1", "sample2", and "sample3".