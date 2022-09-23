"""
Methods to support gene sets.
"""

import os
import re

def find_geneset_pkl(path, path_re='.*\.pkl'):
    """
    Locate PKL files in a given path.

    :param path: Path to search.
    :param path_re: Path pattern to expect.

    :return: A list of files in `path` that match pattern `path_re`.
    """

    path_list = list()

    # Find PKL files
    for root, dir_list, file_list in os.walk(path):
        for pkl_file in file_list:
            pkl_file = os.path.join(root, pkl_file)

            if re.match(path_re, pkl_file):
                path_list.append(pkl_file)

    return path_list
