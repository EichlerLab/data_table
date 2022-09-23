"""
Fst function definitions.
"""

import numpy as np

def fst_wc(ac, an):
    """
    Weir and Cockerham's Fst.
    
    Code from adpted from Benson.

    :param ac: A 1d array of allele counts for each population. This is the total number of alleles (e.g. ALT alleles).
    :param an: A 1d array of allele number for each population. This is the total number of alleles (e.g. 2 * sample_count - no_call_alleles).
    """
    
    # Check ac and an
    if not np.issubdtype(ac.dtype, np.integer):
        raise RuntimeError('ac must be an integer')
    
    if not np.issubdtype(an.dtype, np.integer):
        raise RuntimeError('an must be an integer')
    
    if ac.shape[0] != an.shape[0]:
        raise RuntimeError('ac ({}) an ({}) length mismatch'.format(ac.shape[0], an.shape[0]))
    
    if np.all(ac == an):  # No Fst calculation for AF == 1 for all populations
        return np.float64(0.0)
    
    # Allele frequency
    af = ac / an
    
    # calcuate Ht = 2 * p_bar * q_bar; p refer to the sum of ingroup allele freqs. q is for the same quality for outgroup.
    # calculate Hs = (1/n) * sum of 2 * pi * qi, for i  in range(1,n+1) and n =numPop
    
    # Get number of unique alleles (r) and total number of alleles (n)
    r = an.shape[0]
    n = np.sum(an)
    
    p_bar = np.sum((an / n) * af)
    
    n_bar = n / r

    S_sq_part1 = 1.0 / ((r - 1) * n_bar)
    S_sq_part2 = np.sum(an * (af - p_bar)**2)
    
    sum_n_i_sq = np.sum(an**2)

    S_sq = S_sq_part1 * S_sq_part2

    T1 = S_sq - (1.0 / (2*n_bar - 1)) * (p_bar * (1-p_bar) - (r-1)/r * S_sq)

    nc_part1 = 1.0 / (r - 1)
    nc_part2 = (n - (sum_n_i_sq / n))

    nc = nc_part1 * nc_part2

    T2_part1 = (2 * nc - 1) / (2 * n_bar - 1)
    T2_part2 = p_bar * (1- p_bar)
    T2_part3 = (1 + (2*(r-1)*(n_bar - nc)) / (2 * n_bar - 1))
    T2_part4 = S_sq / r

    T2 = T2_part1 * T2_part2 + T2_part3 * T2_part4

    locus_fst = T1 / T2
    
    # Return value (0 if negative)
    return np.max([locus_fst, 0])


def fst_wright(ac, an):
    """
    Wright Fst.
    
    :param ac: A 1d array of allele counts for each population. This is the total number of alleles (e.g. ALT alleles).
    :param an: A 1d array of allele number for each population. This is the total number of alleles (e.g. 2 * sample_count - no_call_alleles).
    """
    
    # Check ac and an
    if not np.issubdtype(ac.dtype, np.integer):
        raise RuntimeError('ac must be an integer')
    
    if not np.issubdtype(an.dtype, np.integer):
        raise RuntimeError('an must be an integer')
    
    if ac.shape[0] != an.shape[0]:
        raise RuntimeError('ac ({}) an ({}) length mismatch'.format(ac.shape[0], an.shape[0]))
    
    if np.all(ac == an):  # No Fst calculation for AF == 1 for all populations
        return np.float64(0.0)
    
    af = ac / an

    an_t = np.sum(an)

    p = np.sum(ac) / np.sum(an)

    h_t = p * (1 - p)

    h_s = np.sum(an / an_t * af * (1 - af))

    return (h_t - h_s) / h_t

def fst_wc(ac, an):
    """
    Weir and Cockerham's Fst.
    
    Code from adapted from Benson.

    :param ac: A 1d array of allele counts for each population. This is the total number of alleles (e.g. ALT alleles).
    :param an: A 1d array of allele number for each population. This is the total number of alleles (e.g. 2 * sample_count - no_call_alleles).
    """
    
    # Check ac and an
    if not np.issubdtype(ac.dtype, np.integer):
        raise RuntimeError('ac must be an integer')
    
    if not np.issubdtype(an.dtype, np.integer):
        raise RuntimeError('an must be an integer')
    
    if ac.shape[0] != an.shape[0]:
        raise RuntimeError('ac ({}) an ({}) length mismatch'.format(ac.shape[0], an.shape[0]))
    
    if np.all(ac == an):  # No Fst calculation for AF == 1 for all populations
        return np.float64(0.0)
    
    # Allele frequency
    af = ac / an
    
    # calcuate Ht = 2 * p_bar * q_bar; p refer to the sum of ingroup allele freqs. q is for the same quality for outgroup.
    # calculate Hs = (1/n) * sum of 2 * pi * qi, for i  in range(1,n+1) and n =numPop
    
    # Get number of unique alleles (r) and total number of alleles (n)
    r = an.shape[0]
    n = np.sum(an)
    
    p_bar = np.sum((an / n) * af)
    
    n_bar = n / r

    S_sq_part1 = 1.0 / ((r - 1) * n_bar)
    S_sq_part2 = np.sum(an * (af - p_bar)**2)
    
    sum_n_i_sq = np.sum(an**2)

    S_sq = S_sq_part1 * S_sq_part2

    T1 = S_sq - (1.0 / (2*n_bar - 1)) * (p_bar * (1-p_bar) - (r-1)/r * S_sq)

    nc_part1 = 1.0 / (r - 1)
    nc_part2 = (n - (sum_n_i_sq / n))

    nc = nc_part1 * nc_part2

    T2_part1 = (2 * nc - 1) / (2 * n_bar - 1)
    T2_part2 = p_bar * (1- p_bar)
    T2_part3 = (1 + (2*(r-1)*(n_bar - nc)) / (2 * n_bar - 1))
    T2_part4 = S_sq / r

    T2 = T2_part1 * T2_part2 + T2_part3 * T2_part4

    locus_fst = T1 / T2
    
    # Return value (0 if negative)
    return np.max([locus_fst, 0])


# Original code from Benson
def calc_SinglePair_Fst(list_popsizes, list_freqs):
    # calcuate Ht = 2 * p_bar * q_bar; p refer to the sum of ingroup allele freqs. q is for the same quality for outgroup.
    # calculate Hs = (1/n) * sum of 2 * pi * qi, for i  in range(1,n+1) and n =numPop
    
    if np.sum(list_popsizes) == 0:
        return 0

    list_popsizes = [float(x) for x in list_popsizes]
    list_freqs = [float(x) for x in list_freqs]
#   print (n1,n2,p1,p2)

    r = float(len(list_popsizes))
    n = float(sum(list_popsizes))

    p_bar = 0.0 
    for i, v in enumerate(list_freqs):
        p_bar += (list_popsizes[i] / n) * list_freqs[i]

    n_bar = sum(list_popsizes) / r 

    S_sq_part1 = 1.0 / ((r - 1) * n_bar)

    S_sq_part2 = 0.0 
    sum_n_i = 0.0 
    sum_n_i_sq = 0.0 
    for i, v in enumerate(list_freqs):
        S_sq_part2 += list_popsizes[i] * ((list_freqs[i] - p_bar)**2)
        sum_n_i += list_popsizes[i]
        sum_n_i_sq += (list_popsizes[i])**2


    S_sq = S_sq_part1 * S_sq_part2

    T1 = S_sq - (1.0 / (2*n_bar - 1)) * (p_bar * (1-p_bar) - (r-1)/r * S_sq)

    nc_part1 = 1.0 / (r - 1)
    nc_part2 = (sum_n_i - (sum_n_i_sq / sum_n_i))
    nc = nc_part1 * nc_part2

    T2_part1 = (2 * nc - 1) / (2 * n_bar - 1)
    T2_part2 = p_bar * (1- p_bar)
    T2_part3 = (1 + (2*(r-1)*(n_bar - nc)) / (2 * n_bar - 1)) 
    T2_part4 = S_sq / r 
    T2 = T2_part1 * T2_part2 + T2_part3 * T2_part4

    try:
        locus_fst = T1 / T2
    except ZeroDivisionError:
        locus_fst = 0.0
        
    return locus_fst
