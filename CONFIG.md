# Data Table Configuration

Configuration is in a JSON file `config/config.json`.

Basic example:
```

```

## Top level configuration Elements ##

* svpop_dir: Directory where SV-Pop is installed.
* svpop_run_dir: Directory where SV-Pop was run for this project. Variants and annotations are read from this location.
* reference: Path to the reference FASTA file.
* table_def: A dictionary of table definitions. One entry per variant set. Most configurations have just one entry
  in `table_def`.

## Table Definitions ##

The name of a table definition may not contain a dash character ("-"), this is reserved for table definitions with
wildcards (see below).

* sourcetype
* sourcename
* sample
* filter
* svset
* callable_filter
* shift_col
* sections
* support

### Table definitions with wildcards ###

Wildcard values are in the table definition name after a single dash character ("-"). This facility allows for creating
a table for a number of inputs, for example, one for each sample using the same table definition.

For example, if "freeze1-SAMPLE1" is requested, the table definition for "freeze1" is searched for in the table
definitions, and "SAMPLE1" is used as a replacement value for wildcards inside the definition (e.g. likely "{sample}"
in this case).

In the examlpe above, "freeze1-SAMPLE1" is split into "table_name" ("freeze1") and a wildcard string ("SAMPLE1"). More
than one wildcard may be present, for example "freeze1-pav-SAMPLE1" is split into "table_name" ("freeze1") and a
wildcard string ("pav-SAMPLE1"), which by default is split into two separate wildcard values on "-" ("pav" and
"SAMPLE").

Additional options to the table definition for splitting the wildcard string and matching them to keys:
* wildcard_delim [string, default = "-"]: Delimiter for splitting the wildcards into a list. 
* wildcard_limit [int, default = -1]: Split the wildcard string
* wildcard_list [list(string)]: A list of wildcard keys (e.g. ["sample"] or ["callsource", "sample"]). The length
  must match the length of the wildcard list split from the table definition name.

Changing the delimiter and the limit helps for values that might contain the delimiter, for example "freeze1-SAMPLE-1"
would be split into two wildcards by default ("SAMPLE" and "1"), however, if the sample was intended to be "SAMPLE-1",
then either the delimiter could be changed to something not found in the sample name (wildcard_delim = ":") or the
limit could be set (e.g. wildcard_limit = 1).

The wildcard_list contains the names of the keys for the wildcards while everything following the first "-" in the table
definition is the values those keys replace. For examlpe, if the table definition is "freeze1-pav-SAMPLE1" and
wildcard_list is ["callsounce", "sample"], then inside the table definition, "{callsource}" is replaced with "pav" and
"{sample}" is replaced with "SAMPLE1".

### sections ###

Each entry in `sections` is a dictionary with key `{vartype}_{svtype}` and value is a list of keywords. Most are
annotations pulled from SV-Pop and must be run in SV-Pop before running the data table.

Annotations include:
* ref_trf
* refseq
* refseq_prox-2000
* ccre2020
* gc
* homop
* homop_nearest
* dinucl
* dinucl_nearest
* win-LEN: Coordinate of a window of LEN bases around the variant. LEN may be number of bases (e.g. "200") and has
  optional multipliers including "k" and "m" (e.g. "2k").

### shift_col ###

Rearrange column ordering.

* insert_after [string]: Columns are moved after this column.
* move_cols [list(string)]: List of columns to move after the `insert_after` column in the order they should appear
  after moving.


