# import pytest
import allel
# from allel.test.tools import assert_array_equal as aeq
import numpy as np
from archaic_match.funcmodule import get_informative_sites
from archaic_match.__main__ import calc_match_pct


genotype_array = allel.GenotypeArray([
    [[0, 0], [0, 0], [0, 0]],
    [[0, 1], [0, 0], [0, 0]],
    [[1, 0], [0, 0], [0, 0]],
    [[1, 1], [0, 0], [0, 0]],
    [[0, 0], [0, 1], [0, 0]],
    [[0, 1], [0, 1], [0, 0]],
    [[1, 0], [0, 1], [0, 0]],
    [[1, 1], [0, 1], [0, 0]],
    [[0, 0], [1, 0], [0, 0]],
    [[0, 1], [1, 0], [0, 0]],
    [[1, 0], [1, 0], [0, 0]],
    [[1, 1], [1, 0], [0, 0]],
    [[0, 0], [1, 1], [0, 0]],
    [[0, 1], [1, 1], [0, 0]],
    [[1, 0], [1, 1], [0, 0]],
    [[1, 1], [1, 1], [0, 0]],

    [[0, 0], [0, 0], [0, 1]],
    [[0, 1], [0, 0], [0, 1]],
    [[1, 0], [0, 0], [0, 1]],
    [[1, 1], [0, 0], [0, 1]],
    [[0, 0], [0, 1], [0, 1]],
    [[0, 1], [0, 1], [0, 1]],
    [[1, 0], [0, 1], [0, 1]],
    [[1, 1], [0, 1], [0, 1]],
    [[0, 0], [1, 0], [0, 1]],
    [[0, 1], [1, 0], [0, 1]],
    [[1, 0], [1, 0], [0, 1]],
    [[1, 1], [1, 0], [0, 1]],
    [[0, 0], [1, 1], [0, 1]],
    [[0, 1], [1, 1], [0, 1]],
    [[1, 0], [1, 1], [0, 1]],
    [[1, 1], [1, 1], [0, 1]],

    [[0, 0], [0, 0], [1, 0]],
    [[0, 1], [0, 0], [1, 0]],
    [[1, 0], [0, 0], [1, 0]],
    [[1, 1], [0, 0], [1, 0]],
    [[0, 0], [0, 1], [1, 0]],
    [[0, 1], [0, 1], [1, 0]],
    [[1, 0], [0, 1], [1, 0]],
    [[1, 1], [0, 1], [1, 0]],
    [[0, 0], [1, 0], [1, 0]],
    [[0, 1], [1, 0], [1, 0]],
    [[1, 0], [1, 0], [1, 0]],
    [[1, 1], [1, 0], [1, 0]],
    [[0, 0], [1, 1], [1, 0]],
    [[0, 1], [1, 1], [1, 0]],
    [[1, 0], [1, 1], [1, 0]],
    [[1, 1], [1, 1], [1, 0]],

    [[0, 0], [0, 0], [1, 1]],
    [[0, 1], [0, 0], [1, 1]],
    [[1, 0], [0, 0], [1, 1]],
    [[1, 1], [0, 0], [1, 1]],
    [[0, 0], [0, 1], [1, 1]],
    [[0, 1], [0, 1], [1, 1]],
    [[1, 0], [0, 1], [1, 1]],
    [[1, 1], [0, 1], [1, 1]],
    [[0, 0], [1, 0], [1, 1]],
    [[0, 1], [1, 0], [1, 1]],
    [[1, 0], [1, 0], [1, 1]],
    [[1, 1], [1, 0], [1, 1]],
    [[0, 0], [1, 1], [1, 1]],
    [[0, 1], [1, 1], [1, 1]],
    [[1, 0], [1, 1], [1, 1]],
    [[1, 1], [1, 1], [1, 1]],
])

archaic_sample_idx = [0]
modern_sample_idx = [1, 2]

subpops = {'archaic': archaic_sample_idx,
           'modern': modern_sample_idx}

allele_counts = genotype_array.count_alleles_subpops(subpops=subpops)

modern_haplotypes = allel.HaplotypeArray(
    genotype_array[:, modern_sample_idx]
    .flatten().reshape(-1, len(modern_sample_idx) * 2))
archic_haplotypes = allel.HaplotypeArray(
    genotype_array[:, archaic_sample_idx]
    .flatten().reshape(-1, len(archaic_sample_idx) * 2))

archaic_allele_counts = allel.AlleleCountsArray([
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2],
    [2, 0],
    [1, 1],
    [1, 1],
    [0, 2]
])


modern_allele_counts = allel.AlleleCountsArray([
    [4, 0],
    [4, 0],
    [4, 0],
    [4, 0],
    [3, 1],
    [3, 1],
    [3, 1],
    [3, 1],
    [3, 1],
    [3, 1],
    [3, 1],
    [3, 1],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],

    [3, 1],
    [3, 1],
    [3, 1],
    [3, 1],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [1, 3],
    [1, 3],
    [1, 3],
    [1, 3],

    [3, 1],
    [3, 1],
    [3, 1],
    [3, 1],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [1, 3],
    [1, 3],
    [1, 3],
    [1, 3],

    [2, 2],
    [2, 2],
    [2, 2],
    [2, 2],
    [1, 3],
    [1, 3],
    [1, 3],
    [1, 3],
    [1, 3],
    [1, 3],
    [1, 3],
    [1, 3],
    [0, 4],
    [0, 4],
    [0, 4],
    [0, 4]
])


informative_sites_derived_in_archaic = [
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True,
    False, True, True, True
]

informative_sites_derived_in_archaic_or_modern = [
    False, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True,
    True, True, True, True
]


g_test = allel.GenotypeArray([
    [[0, 0], [0, 1]],
    [[0, 2], [1, 1]],
    [[2, 2], [-1, -1]]
])

c_data = [
    [3, 1, 0],
    [1, 2, 1],
    [0, 0, 2]
]
c_test = allel.AlleleCountsArray([
    [3, 1, 0],
    [1, 2, 1],
    [0, 0, 2]
])


def test_count_allele():
    print(c_test)
    print(g_test.count_alleles())
    assert np.array_equal(np.asarray(g_test.count_alleles()),
                          np.asarray(c_test))


def test_count_allele_subpops():
    assert np.array_equal(allele_counts['archaic'],
                          np.asarray(archaic_allele_counts))
    assert np.array_equal(allele_counts['modern'],
                          np.asarray(modern_allele_counts))


def test_get_informative_sites_derived_in_archaic():
    assert np.array_equal(
        get_informative_sites(allele_counts, "derived_in_archaic"),
        informative_sites_derived_in_archaic)


def test_get_informative_sites_derived_in_archaic_or_modern():
    assert np.array_equal(
        get_informative_sites(allele_counts, "derived_in_archaic_or_modern"),
        informative_sites_derived_in_archaic_or_modern)


def test_calc_match_pct_derived_in_archaic():
    assert (calc_match_pct(
            informative_sites_derived_in_archaic,
            archic_haplotypes,
            modern_haplotypes.T[1, ]
            ) == (24 / 48))


def test_calc_match_pct_derived_in_archaic_or_modern():
    assert (calc_match_pct(
            informative_sites_derived_in_archaic_or_modern,
            archic_haplotypes,
            modern_haplotypes.T[1, ]
            ) == (24 / 63))
