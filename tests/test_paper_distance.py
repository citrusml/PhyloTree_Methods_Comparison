#!/usr/bin/env python3
"""
Unit tests for src/paper_distance.py, bin/run_pwa_nj.py, and bin/run_msa_nj.py.

Verifies the mathematical properties of:
1. WAG + Gamma4 spectral decomposition and transition probability matrix P(t).
2. Maximum Likelihood pairwise evolutionary distance estimation under WAG + Gamma4.
3. Complete gap deletion (NJ_MSA_c / exclude='all').
4. MMseqs2 symmetrical relative sequence similarity score transformation to distance matrix.
"""

from __future__ import annotations

import os
import sys
import math
import numpy as np
import pytest

_PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _PROJECT_ROOT not in sys.path:
    sys.path.insert(0, _PROJECT_ROOT)
_BIN_DIR = os.path.join(_PROJECT_ROOT, "bin")
if _BIN_DIR not in sys.path:
    sys.path.insert(0, _BIN_DIR)

from src.paper_distance import (
    WAG_PI,
    WAGGammaModel,
    apply_complete_gap_deletion,
    calculate_wag_gamma_ml_distance,
    similarity_to_distance_matrix,
)
from run_pwa_nj import calculate_distance


def test_discrete_gamma_rates_mean_is_one() -> None:
    """
    Verifies that the 4-category discrete Gamma rates have mean 1.0 and positive ascending values.

    Returns
    -------
    None
    """
    rates = WAGGammaModel._compute_discrete_gamma_rates(alpha=1.0, k=4)
    assert rates.shape == (4,)
    assert np.all(rates > 0.0)
    assert np.all(np.diff(rates) > 0.0)
    assert math.isclose(float(np.mean(rates)), 1.0, rel_tol=1e-9)


def test_wag_rate_matrix_normalization_and_stationarity() -> None:
    """
    Verifies that the WAG instantaneous rate matrix Q satisfies:
    1. Row sums equal 0.
    2. Stationary distribution pi * Q = 0.
    3. Expected substitution rate -sum_i pi_i Q_ii = 1.0 substitutions/site.
    4. Transition matrix P(t) has row sums equal to 1.0.

    Returns
    -------
    None
    """
    model = WAGGammaModel(alpha=1.0, num_categories=4)
    q_mat = model.q_matrix

    # Row sums of Q must be 0
    assert np.allclose(np.sum(q_mat, axis=1), 0.0, atol=1e-12)
    # Detailed balance / stationarity: pi @ Q = 0
    assert np.allclose(model.pi @ q_mat, 0.0, atol=1e-12)
    # Normalized rate: -sum_i pi_i Q_ii = 1.0
    mean_rate = -float(np.sum(model.pi * np.diag(q_mat)))
    assert math.isclose(mean_rate, 1.0, rel_tol=1e-12)

    # Transition matrix at t = 0.5 substitutions/site must be row-stochastic
    p_mat = model.transition_matrix(0.5)
    assert p_mat.shape == (20, 20)
    assert np.all(p_mat > 0.0)
    assert np.allclose(np.sum(p_mat, axis=1), 1.0, atol=1e-12)


def test_wag_gamma_ml_distance_monotonicity() -> None:
    """
    Verifies that calculate_wag_gamma_ml_distance returns 0.0 for identical sequences
    and increases monotonically as more amino acid substitutions are introduced.

    Returns
    -------
    None
    """
    seq_base = "ACDEFGHIKLMNPQRSTVWY" * 10  # 200 aa
    d_zero = calculate_wag_gamma_ml_distance(seq_base, seq_base, alpha=1.0)
    assert d_zero == 0.0

    # 10% substitutions (Conservative I -> V)
    seq_mut1 = ("ACDEFGHVKLMNPQRSTVWY" * 10)
    d_mut1 = calculate_wag_gamma_ml_distance(seq_base, seq_mut1, alpha=1.0)

    # 20% substitutions (I -> W and A -> Y)
    seq_mut2 = ("YCDEFGHWKLMNPQRSTVWY" * 10)
    d_mut2 = calculate_wag_gamma_ml_distance(seq_base, seq_mut2, alpha=1.0)

    assert 0.0 < d_mut1 < d_mut2
    # Also verify dispatch via run_pwa_nj.calculate_distance
    assert math.isclose(
        calculate_distance(seq_base, seq_mut1, dist_model="wag_gamma", alpha=1.0),
        d_mut1,
        rel_tol=1e-9,
    )


def test_complete_gap_deletion() -> None:
    """
    Verifies that apply_complete_gap_deletion removes any MSA column containing '-'
    in at least one sequence (matching NJ_MSA_c / exclude='all').

    Returns
    -------
    None
    """
    msa = [
        "AC-DEFG",
        "ACDDE-G",
        "ACDDEGG",
    ]
    filtered = apply_complete_gap_deletion(msa)
    assert filtered == ["ACDEG", "ACDEG", "ACDEG"]


def test_similarity_to_distance_matrix() -> None:
    """
    Verifies linear (1 - w_ij) and log (-ln(w_ij)) transformations of MMseqs2
    similarity matrix w_ij.

    Returns
    -------
    None
    """
    w_mat = np.array([
        [1.0, 0.8, 0.0],
        [0.8, 1.0, 0.25],
        [0.0, 0.25, 1.0],
    ], dtype=np.float64)

    d_lin = similarity_to_distance_matrix(w_mat, dist_model="mmseqs_linear")
    assert np.allclose(np.diag(d_lin), 0.0)
    assert math.isclose(float(d_lin[0, 1]), 0.2, rel_tol=1e-9)
    assert math.isclose(float(d_lin[0, 2]), 1.0, rel_tol=1e-9)

    d_log = similarity_to_distance_matrix(w_mat, dist_model="mmseqs_log", epsilon=1e-4)
    assert np.allclose(np.diag(d_log), 0.0)
    assert math.isclose(float(d_log[0, 1]), -math.log(0.8), rel_tol=1e-9)
    assert math.isclose(float(d_log[0, 2]), -math.log(1e-4), rel_tol=1e-9)
