"""
Build a table of variants and annotations.
"""

import collections
import gzip
import numpy as np
import os
import pandas as pd
import pickle
import re
import sys

from Bio import SeqIO
import Bio.bgzf

DATA_TABLE_DIR = os.path.dirname(workflow.snakefile)


#
# Config
#

SITE_CONFIG_FILE_NAME = os.path.join(DATA_TABLE_DIR, 'local/config/config.json')
RUN_CONFIG_FILE_NAME = 'config/config.json'

if os.path.isfile(SITE_CONFIG_FILE_NAME):
    configfile: SITE_CONFIG_FILE_NAME

if os.path.isfile(RUN_CONFIG_FILE_NAME):
    configfile: RUN_CONFIG_FILE_NAME


# Check SVPOP_DIR

if 'svpop_dir' not in config:
    raise RuntimeError('Missing "svpop_dir" in config (path to SV-Pop pipeline installation)')

if 'svpop_run_dir' not in config:
    raise RuntimeError('Missing "svpop_run_dir" in config (path to SV-Pop run directory where data are found)')

if 'table_def' not in config:
    raise RuntimeError('Missing "table_def" in config (output table definitions)')

if 'reference' not in config:
    raise RuntimeError('Missing "reference" in config (output table definitions)')


# Get reference (for convenience)

REF_FA = config['reference']
REF_FAI = config['reference'] + '.fai'

# Get run directories

SVPOP_RUN_DIR = config['svpop_run_dir']


#
# Library imports
#

sys.path.append(config['svpop_dir'])
sys.path.append(os.path.join(config['svpop_dir'], 'dep'))

import svpoplib
import dtablib


#
# Rules
#

include: 'rules/base_table.snakefile'
include: 'rules/extern.snakefile'
include: 'rules/genes.snakefile'
include: 'rules/geneset.snakefile'
include: 'rules/pop.snakefile'
include: 'rules/regions.snakefile'
include: 'rules/snv.snakefile'
include: 'rules/support.snakefile'
include: 'rules/tracks.snakefile'
include: 'rules/interseq.snakefile'
include: 'rules/mapping.snakefile'
include: 'rules/vcf.snakefile'
