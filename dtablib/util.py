"""
General utilities.
"""

import collections
import numpy as np
import intervaltree


def bed_to_tree_dict(df_bed):
    """
    Turn a BED file into a dictionary keyed by chromosome names containing an interval tree for each chromosome.

    :param df_bed: BED file as a Pandas DataFrame.

    :return: collections.defaultdict of intervals.
    """

    bed_tree = collections.defaultdict(intervaltree.IntervalTree)

    for index, row in df_bed.iterrows():
        bed_tree[row['#CHROM']][row['POS']:row['END']] = index

    return bed_tree


def bed_match_dict_region(df, bed_tree):

    return df.apply(lambda row:
        len(bed_tree[row['#CHROM']][row['POS']:row['END']]) > 0,
        axis=1
    )

    # return df.apply(lambda row:
    #     np.any(
    #         [
    #             (
    #                 (
    #                     min(row['END'], bed_region.end) - max(row['POS'], bed_region.begin)
    #                 ) / row['SVLEN'] > overlap
    #             ) for bed_region in bed_tree[row['#CHROM']][row['POS']:row['END']]
    #         ]
    #     ),
    #     axis=1
    # )


def get_bool(val):
    """
    Get a boolean from a value.

    :param val: Value.

    :return: `True` or `False`.
    """

    val_lower = val.lower()

    if val in {'t', 'true', 'yes', 1, '1', True}:
        return True

    if val in {'f', 'false', 'no', 0, '0', False}:
        return False

    raise RuntimeError('Unrecognized value for bool conversion: {}'.format(val))
