#!/usr/bin/env python3
"""
Paper-Exact Evolutionary Distance and Sequence Similarity Calculation Module.

Implements the exact pairwise and multiple-sequence-alignment distance matrix
calculations specified in Matsui & Iwasaki (2020, Systematic Biology, syz049):

1. NJ_PSA (Pairwise Sequence Alignment + Neighbor-Joining):
   - All-to-all pairwise local alignment via MMseqs2 (-s 7.5, -e 10).
   - Symmetrical relative sequence similarity score (Dufour et al. 2010):
       w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))
   - Distance matrix transformations from w_ij in [0, 1]:
     * 'mmseqs_linear' (or 'mmseqs_sss'): d_ij = 1.0 - w_ij (linear dissimilarity)
     * 'mmseqs_log':                      d_ij = -ln(max(w_ij, epsilon))

2. NJ_MSA_p & NJ_MSA_c / Pairwise WAG + Gamma4 Maximum Likelihood Distance:
   - Exact Whelan & Goldman (2001) empirical substitution rate matrix (WAG) +
     4-category discrete Gamma rate heterogeneity (Gamma4, alpha=1.0), matching
     R phangorn::dist.ml(model="WAG", exclude="pairwise" | "all") and INDELible v1.03.
   - Supports both gap exclusion modes:
     * 'wag_gamma' (or 'wag_gamma_pairwise'): Pairwise gap deletion (NJ_MSA_p / PSA)
     * 'wag_gamma_complete':                  Complete column gap deletion (NJ_MSA_c)

3. PhyPA-style Needleman-Wunsch Pairwise Alignment + Poisson / Gamma-Poisson Distance:
   - Global pairwise alignment using BLOSUM62 (Gap Open = 20, Gap Extend = 2)
     followed by Poisson ('poisson'), Gamma-Poisson ('gamma_poisson'), or
     WAG+Gamma4 ML ('wag_gamma') distance estimation.
"""

from __future__ import annotations

import math
import os
import shutil
import subprocess
import sys
import tempfile
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from Bio import SeqIO
from scipy.optimize import minimize_scalar
from scipy.special import gammainc
from scipy.stats import gamma as scipy_gamma

# Standard PAML / INDELible / phangorn 20 amino acid ordering
AA_ORDER: str = "ARNDCQEGHILKMFPSTWYV"
AA_TO_IDX: Dict[str, int] = {aa: idx for idx, aa in enumerate(AA_ORDER)}

# Lookup table mapping ASCII byte (0..255) to amino acid index (0..19), or -1 for gap/unknown
_ASCII_TO_AA_IDX: np.ndarray = np.full(256, -1, dtype=np.int8)
for _aa, _idx in AA_TO_IDX.items():
    _ASCII_TO_AA_IDX[ord(_aa)] = _idx
    _ASCII_TO_AA_IDX[ord(_aa.lower())] = _idx

# Official Whelan & Goldman (2001) equilibrium amino acid frequencies (unitless probabilities)
WAG_PI: np.ndarray = np.array(
    [
        0.08662790,  # A
        0.04397200,  # R
        0.03908940,  # N
        0.05704510,  # D
        0.01930780,  # C
        0.03672810,  # Q
        0.05805890,  # E
        0.08325180,  # G
        0.02443130,  # H
        0.04846600,  # I
        0.08620970,  # L
        0.06202860,  # K
        0.01950273,  # M
        0.03843190,  # F
        0.04576310,  # P
        0.06951790,  # S
        0.06101270,  # T
        0.01438590,  # W
        0.03527420,  # Y
        0.07089560,  # V
    ],
    dtype=np.float64,
)
WAG_PI = WAG_PI / np.sum(WAG_PI)

# Official Whelan & Goldman (2001) symmetric exchangeability matrix S (20x20, relative rate units)
WAG_EXCHANGEABILITY: np.ndarray = np.array(
    [
        [0.0000000, 0.5515710, 0.5098480, 0.7389980, 1.0270400, 0.9085980, 1.5828500, 1.4167200, 0.3169540, 0.1933350, 0.3979150, 0.9062650, 0.8934960, 0.2104940, 1.4385500, 3.3707900, 2.1211100, 0.1131330, 0.2407350, 2.0060100],
        [0.5515710, 0.0000000, 0.6353460, 0.1473040, 0.5281910, 3.0355000, 0.4391570, 0.5846650, 2.1371500, 0.1869790, 0.4976710, 5.3514200, 0.6831620, 0.1027110, 0.6794890, 1.2241900, 0.5544130, 1.1639200, 0.3815330, 0.2518490],
        [0.5098480, 0.6353460, 0.0000000, 5.4294200, 0.2652560, 1.5436400, 0.9471980, 1.1255600, 3.9562900, 0.5542360, 0.1315280, 3.0120100, 0.1982210, 0.0961621, 0.1950810, 3.9742300, 2.0300600, 0.0719167, 1.0860000, 0.1962460],
        [0.7389980, 0.1473040, 5.4294200, 0.0000000, 0.0302949, 0.6167830, 6.1741600, 0.8655840, 0.9306760, 0.0394370, 0.0848047, 0.4798550, 0.1037540, 0.0467304, 0.4239840, 1.0717600, 0.3748660, 0.1297670, 0.3257110, 0.1523350],
        [1.0270400, 0.5281910, 0.2652560, 0.0302949, 0.0000000, 0.0988179, 0.0213520, 0.3066740, 0.2489720, 0.1701350, 0.3842870, 0.0740339, 0.3904820, 0.3980200, 0.1094040, 1.4076600, 0.5129840, 0.7170700, 0.5438330, 1.0021400],
        [0.9085980, 3.0355000, 1.5436400, 0.6167830, 0.0988179, 0.0000000, 5.4694700, 0.3300520, 4.2941100, 0.1139170, 0.8694890, 3.8949000, 1.5452600, 0.0999208, 0.9333720, 1.0288700, 0.8579280, 0.2157370, 0.2277100, 0.3012810],
        [1.5828500, 0.4391570, 0.9471980, 6.1741600, 0.0213520, 5.4694700, 0.0000000, 0.5677170, 0.5700250, 0.1273950, 0.1542630, 2.5844300, 0.3151240, 0.0811339, 0.6823550, 0.7049390, 0.8227650, 0.1565570, 0.1963030, 0.5887310],
        [1.4167200, 0.5846650, 1.1255600, 0.8655840, 0.3066740, 0.3300520, 0.5677170, 0.0000000, 0.2494100, 0.0304501, 0.0613037, 0.3735580, 0.1741000, 0.0499310, 0.2435700, 1.3418200, 0.2258330, 0.3369830, 0.1036040, 0.1872470],
        [0.3169540, 2.1371500, 3.9562900, 0.9306760, 0.2489720, 4.2941100, 0.5700250, 0.2494100, 0.0000000, 0.1381900, 0.4994620, 0.8904320, 0.4041410, 0.6793710, 0.6961980, 0.7401690, 0.4733070, 0.2625690, 3.8734400, 0.1183580],
        [0.1933350, 0.1869790, 0.5542360, 0.0394370, 0.1701350, 0.1139170, 0.1273950, 0.0304501, 0.1381900, 0.0000000, 3.1709700, 0.3238320, 4.2574600, 1.0594700, 0.0999288, 0.3194400, 1.4581600, 0.2124830, 0.4201700, 7.8213000],
        [0.3979150, 0.4976710, 0.1315280, 0.0848047, 0.3842870, 0.8694890, 0.1542630, 0.0613037, 0.4994620, 3.1709700, 0.0000000, 0.2575550, 4.8540200, 2.1151700, 0.4158440, 0.3447390, 0.3266220, 0.6653090, 0.3986180, 1.8003400],
        [0.9062650, 5.3514200, 3.0120100, 0.4798550, 0.0740339, 3.8949000, 2.5844300, 0.3735580, 0.8904320, 0.3238320, 0.2575550, 0.0000000, 0.9342760, 0.0888360, 0.5568960, 0.9671300, 1.3869800, 0.1375050, 0.1332640, 0.3054340],
        [0.8934960, 0.6831620, 0.1982210, 0.1037540, 0.3904820, 1.5452600, 0.3151240, 0.1741000, 0.4041410, 4.2574600, 4.8540200, 0.9342760, 0.0000000, 1.1906300, 0.1713290, 0.4939050, 1.5161200, 0.5157060, 0.4284370, 2.0584500],
        [0.2104940, 0.1027110, 0.0961621, 0.0467304, 0.3980200, 0.0999208, 0.0811339, 0.0499310, 0.6793710, 1.0594700, 2.1151700, 0.0888360, 1.1906300, 0.0000000, 0.1614440, 0.5459310, 0.1719030, 1.5296400, 6.4542800, 0.6498920],
        [1.4385500, 0.6794890, 0.1950810, 0.4239840, 0.1094040, 0.9333720, 0.6823550, 0.2435700, 0.6961980, 0.0999288, 0.4158440, 0.5568960, 0.1713290, 0.1614440, 0.0000000, 1.6132800, 0.7953840, 0.1394050, 0.2160460, 0.3148870],
        [3.3707900, 1.2241900, 3.9742300, 1.0717600, 1.4076600, 1.0288700, 0.7049390, 1.3418200, 0.7401690, 0.3194400, 0.3447390, 0.9671300, 0.4939050, 0.5459310, 1.6132800, 0.0000000, 4.3780200, 0.5237420, 0.7869930, 0.2327390],
        [2.1211100, 0.5544130, 2.0300600, 0.3748660, 0.5129840, 0.8579280, 0.8227650, 0.2258330, 0.4733070, 1.4581600, 0.3266220, 1.3869800, 1.5161200, 0.1719030, 0.7953840, 4.3780200, 0.0000000, 0.1108640, 0.2911480, 1.3882300],
        [0.1131330, 1.1639200, 0.0719167, 0.1297670, 0.7170700, 0.2157370, 0.1565570, 0.3369830, 0.2625690, 0.2124830, 0.6653090, 0.1375050, 0.5157060, 1.5296400, 0.1394050, 0.5237420, 0.1108640, 0.0000000, 2.4853900, 0.3653690],
        [0.2407350, 0.3815330, 1.0860000, 0.3257110, 0.5438330, 0.2277100, 0.1963030, 0.1036040, 3.8734400, 0.4201700, 0.3986180, 0.1332640, 0.4284370, 6.4542800, 0.2160460, 0.7869930, 0.2911480, 2.4853900, 0.0000000, 0.3147300],
        [2.0060100, 0.2518490, 0.1962460, 0.1523350, 1.0021400, 0.3012810, 0.5887310, 0.1872470, 0.1183580, 7.8213000, 1.8003400, 0.3054340, 2.0584500, 0.6498920, 0.3148870, 0.2327390, 1.3882300, 0.3653690, 0.3147300, 0.0000000],
    ],
    dtype=np.float64,
)


class WAGGammaModel:
    """
    Precomputed spectral decomposition of the WAG + discrete Gamma substitution model.

    Matches R phangorn::dist.ml(model="WAG") and INDELible v1.03 WAG + Gamma4 model.

    Parameters
    ----------
    alpha : float, default=1.0
        Shape parameter of the discrete Gamma distribution (dimensionless).
        If alpha <= 0 or None, homogeneous rates (1 category, rate=1.0) are used.
    num_categories : int, default=4
        Number of discrete Gamma rate categories (count, default: 4).
    """

    def __init__(self, alpha: Optional[float] = 1.0, num_categories: int = 4) -> None:
        self.alpha: Optional[float] = alpha
        self.num_categories: int = num_categories
        self.pi: np.ndarray = WAG_PI.copy()

        # Construct normalized instantaneous rate matrix Q (substitutions per site)
        q_mat = WAG_EXCHANGEABILITY * self.pi[np.newaxis, :]
        np.fill_diagonal(q_mat, 0.0)
        row_sums = np.sum(q_mat, axis=1)
        np.fill_diagonal(q_mat, -row_sums)

        mean_rate = float(np.sum(self.pi * row_sums))
        q_mat /= mean_rate
        self.q_matrix: np.ndarray = q_mat

        # Symmetrize via S_sym = Pi^{1/2} Q Pi^{-1/2} for exact real eigendecomposition
        sqrt_pi = np.sqrt(self.pi)
        inv_sqrt_pi = 1.0 / sqrt_pi
        sym_mat = (sqrt_pi[:, np.newaxis] * q_mat) * inv_sqrt_pi[np.newaxis, :]
        sym_mat = 0.5 * (sym_mat + sym_mat.T)

        eigvals, eigvecs = np.linalg.eigh(sym_mat)
        # Ensure largest eigenvalue is exactly 0.0
        eigvals[-1] = 0.0
        self.eigvals: np.ndarray = eigvals
        # Right eigenvectors V = Pi^{-1/2} U, Left eigenvectors V^{-1} = U^T Pi^{1/2}
        self.v_right: np.ndarray = inv_sqrt_pi[:, np.newaxis] * eigvecs
        self.v_inv: np.ndarray = eigvecs.T * sqrt_pi[np.newaxis, :]

        # Discrete Gamma category rates (Yang 1994 category means)
        self.rates: np.ndarray = self._compute_discrete_gamma_rates(alpha, num_categories)

    @staticmethod
    def _compute_discrete_gamma_rates(alpha: Optional[float], k: int) -> np.ndarray:
        """
        Computes k discrete Gamma category mean rates (unitless, mean = 1.0)
        following Yang (1994) as used in INDELible, RAxML, IQ-TREE, and phangorn.

        Parameters
        ----------
        alpha : float or None
            Gamma shape parameter (dimensionless).
        k : int
            Number of discrete Gamma categories (count).

        Returns
        -------
        np.ndarray
            Array of k relative substitution rates (unitless).
        """
        if alpha is None or alpha <= 0.0 or k <= 1:
            return np.ones(1, dtype=np.float64)

        quantiles = np.linspace(0.0, 1.0, k + 1)
        # Cutoffs of Gamma(alpha, scale=1/alpha)
        cutoffs = scipy_gamma.ppf(quantiles[1:-1], a=alpha, scale=1.0 / alpha)
        # Incomplete gamma integrals I(x * alpha, alpha + 1)
        inc_vals = np.concatenate(
            [[0.0], gammainc(alpha + 1.0, cutoffs * alpha), [1.0]]
        )
        rates = np.diff(inc_vals) * k
        # Normalize to exact mean 1.0
        return rates / np.mean(rates)

    def transition_matrix(self, distance: float) -> np.ndarray:
        """
        Computes the site-averaged 20x20 transition probability matrix P(d)
        at evolutionary distance d (substitutions per site).

        Parameters
        ----------
        distance : float
            Evolutionary branch distance d (expected substitutions per site).

        Returns
        -------
        np.ndarray
            20x20 transition probability matrix P_ab(d) (unitless probabilities).
        """
        d_val = max(0.0, float(distance))
        # exp(lambda_m * r_c * d): shape (num_categories, 20)
        exp_terms = np.exp(
            np.outer(self.rates * d_val, self.eigvals)
        )  # (K, 20)
        mean_exp = np.mean(exp_terms, axis=0)  # (20,)
        p_mat = (self.v_right * mean_exp[np.newaxis, :]) @ self.v_inv
        return np.clip(p_mat, 1e-15, 1.0)


_CACHED_WAG_MODELS: Dict[Tuple[float, int], WAGGammaModel] = {}


def get_wag_gamma_model(alpha: Optional[float] = 1.0, num_categories: int = 4) -> WAGGammaModel:
    """
    Retrieves or creates a cached WAGGammaModel instance.

    Parameters
    ----------
    alpha : float or None, default=1.0
        Gamma shape parameter alpha (dimensionless).
    num_categories : int, default=4
        Number of discrete Gamma categories (count).

    Returns
    -------
    WAGGammaModel
        Cached WAG+Gamma model instance.
    """
    key_alpha = float(alpha) if (alpha is not None and alpha > 0.0) else 0.0
    key = (key_alpha, int(num_categories))
    if key not in _CACHED_WAG_MODELS:
        _CACHED_WAG_MODELS[key] = WAGGammaModel(alpha=key_alpha, num_categories=num_categories)
    return _CACHED_WAG_MODELS[key]


def encode_sequence_to_aa_indices(seq: str) -> np.ndarray:
    """
    Converts an amino acid sequence string into an array of integer indices (0..19),
    with gaps ('-') or ambiguous characters mapped to -1.

    Parameters
    ----------
    seq : str
        Amino acid sequence string.

    Returns
    -------
    np.ndarray
        1D int8 array of length L containing values in {-1, 0, ..., 19}.
    """
    raw_bytes = np.frombuffer(seq.encode("ascii", errors="replace"), dtype=np.uint8)
    return _ASCII_TO_AA_IDX[raw_bytes]


def compute_aa_pair_counts(idx1: np.ndarray, idx2: np.ndarray) -> Tuple[np.ndarray, int, int]:
    """
    Computes the 20x20 amino acid co-occurrence count matrix C_ab between two
    encoded sequences at valid (non-gap) overlapping positions.

    Parameters
    ----------
    idx1 : np.ndarray
        Encoded amino acid indices for sequence 1 (-1 for gap, 0..19 for AA).
    idx2 : np.ndarray
        Encoded amino acid indices for sequence 2 (-1 for gap, 0..19 for AA).

    Returns
    -------
    counts : np.ndarray
        20x20 matrix of amino acid pair counts (count).
    valid_sites : int
        Total number of valid non-gap overlapping sites (sites).
    mismatches : int
        Number of non-identical valid overlapping sites (sites).
    """
    valid_mask = (idx1 >= 0) & (idx2 >= 0)
    valid_sites = int(np.count_nonzero(valid_mask))
    if valid_sites == 0:
        return np.zeros((20, 20), dtype=np.float64), 0, 0

    a1 = idx1[valid_mask].astype(np.int32)
    a2 = idx2[valid_mask].astype(np.int32)
    mismatches = int(np.count_nonzero(a1 != a2))
    flat_indices = a1 * 20 + a2
    counts = np.bincount(flat_indices, minlength=400).reshape((20, 20)).astype(np.float64)
    return counts, valid_sites, mismatches


def calculate_wag_gamma_ml_distance(
    seq1: Union[str, np.ndarray],
    seq2: Union[str, np.ndarray],
    alpha: float = 1.0,
    min_overlap: int = 1,
    max_dist: float = 10.0,
) -> float:
    """
    Computes the Maximum Likelihood evolutionary distance d between two aligned
    amino acid sequences under the WAG + Gamma4 substitution model, matching
    R phangorn::dist.ml(model="WAG").

    Parameters
    ----------
    seq1 : str or np.ndarray
        First aligned sequence (string or pre-encoded int8 index array).
    seq2 : str or np.ndarray
        Second aligned sequence (string or pre-encoded int8 index array).
    alpha : float, default=1.0
        Gamma shape parameter alpha (dimensionless).
    min_overlap : int, default=1
        Minimum required valid non-gap overlapping sites (sites).
        In phangorn::dist.ml, any pair with >= 1 overlapping site is evaluated,
        while 0 overlapping sites yield max_dist fallback.
    max_dist : float, default=10.0
        Maximum ceiling evolutionary distance (substitutions per site).

    Returns
    -------
    float
        Estimated ML evolutionary distance d (substitutions per site).
    """
    idx1 = encode_sequence_to_aa_indices(seq1) if isinstance(seq1, str) else seq1
    idx2 = encode_sequence_to_aa_indices(seq2) if isinstance(seq2, str) else seq2

    counts, valid_sites, mismatches = compute_aa_pair_counts(idx1, idx2)
    if valid_sites < min_overlap:
        return float(max_dist)
    if mismatches == 0:
        return 0.0

    # Symmetrize count matrix for time-reversible likelihood optimization
    sym_counts = 0.5 * (counts + counts.T)
    nz_rows, nz_cols = np.nonzero(sym_counts > 0.0)
    nz_weights = sym_counts[nz_rows, nz_cols]

    model = get_wag_gamma_model(alpha=alpha, num_categories=4)

    def neg_log_likelihood(d_val: float) -> float:
        p_mat = model.transition_matrix(d_val)
        log_p = np.log(p_mat[nz_rows, nz_cols])
        return -float(np.dot(nz_weights, log_p))

    res = minimize_scalar(
        neg_log_likelihood,
        bounds=(1e-6, max_dist),
        method="bounded",
        options={"xatol": 1e-5},
    )
    if not res.success or np.isnan(res.x):
        return float(max_dist)
    return float(np.clip(res.x, 0.0, max_dist))


def apply_complete_gap_deletion(aligned_seqs: List[str]) -> List[str]:
    """
    Removes every column in an MSA that contains any gap or non-standard character
    across any sequence, strictly implementing Matsui & Iwasaki (2020) NJ_MSA_c
    (phangorn::dist.ml exclude="all").

    Parameters
    ----------
    aligned_seqs : List[str]
        List of N aligned amino acid sequences of equal length L (residues).

    Returns
    -------
    List[str]
        List of N gap-free aligned sequences containing only universally present columns.
    """
    if not aligned_seqs:
        return []
    encoded_rows = np.vstack([encode_sequence_to_aa_indices(s) for s in aligned_seqs])
    # Keep columns where all sequences have valid amino acid indices (>= 0)
    keep_col_mask = np.all(encoded_rows >= 0, axis=0)
    if not np.any(keep_col_mask):
        return ["" for _ in aligned_seqs]

    char_matrix = np.array([list(s.upper()) for s in aligned_seqs], dtype="U1")
    filtered_matrix = char_matrix[:, keep_col_mask]
    return ["".join(row) for row in filtered_matrix]


def find_mmseqs_binary() -> str:
    """
    Locates the MMseqs2 executable in PATH or conda/micromamba environments.

    Returns
    -------
    str
        Absolute path to the 'mmseqs' executable.

    Raises
    ------
    RuntimeError
        If 'mmseqs' cannot be located.
    """
    cand = shutil.which("mmseqs")
    if cand and os.path.isfile(cand) and os.access(cand, os.X_OK):
        return cand

    py_dir = os.path.dirname(sys.executable)
    local_cand = os.path.join(py_dir, "mmseqs")
    if os.path.isfile(local_cand) and os.access(local_cand, os.X_OK):
        return local_cand

    env_dir = os.environ.get("CONDA_PREFIX")
    if env_dir:
        env_cand = os.path.join(env_dir, "bin", "mmseqs")
        if os.path.isfile(env_cand) and os.access(env_cand, os.X_OK):
            return env_cand

    user_home = os.path.expanduser("~")
    home_bin = os.path.join(user_home, "bin", "mmseqs")
    if os.path.isfile(home_bin) and os.access(home_bin, os.X_OK):
        return home_bin

    mamba_root = os.environ.get("MAMBA_ROOT_PREFIX")
    roots = [os.path.join(user_home, sub) for sub in [".micromamba", "micromamba", ".conda", "miniforge3", "mambaforge", "miniconda3", "anaconda3"]]
    if mamba_root:
        roots.insert(0, mamba_root)
    if env_dir:
        roots.insert(0, os.path.dirname(os.path.dirname(env_dir)))

    for root in roots:
        for env_name in ["phylomethod_env", "mmseqs_env"]:
            m_cand = os.path.join(root, "envs", env_name, "bin", "mmseqs")
            if os.path.isfile(m_cand) and os.access(m_cand, os.X_OK):
                return m_cand

    raise RuntimeError("MMseqs2 binary 'mmseqs' was not found in PATH. Please ensure mmseqs2 is installed.")


def compute_mmseqs_similarity_matrix(
    fasta_path: str,
    sensitivity: float = 7.5,
    evalue: float = 10.0,
    threads: int = 1,
) -> Tuple[List[str], np.ndarray, Dict[str, float]]:
    """
    Computes the all-to-all MMseqs2 symmetrical relative sequence similarity score matrix w_ij,
    strictly following Matsui & Iwasaki (2020, Systematic Biology, syz049) and gs2 (src/format.cpp):
        w_ij = (S_bits(i, j) + S_bits(j, i)) / (S_bits(i, i) + S_bits(j, j))

    Parameters
    ----------
    fasta_path : str
        Path to input FASTA sequence file (unaligned or aligned; gaps are stripped).
    sensitivity : float, default=7.5
        MMseqs2 sensitivity parameter -s (dimensionless, paper default: 7.5).
    evalue : float, default=10.0
        MMseqs2 E-value threshold -e (dimensionless, paper default: 10.0).
    threads : int, default=1
        Number of CPU threads for MMseqs2 search (count).

    Returns
    -------
    taxa_names : List[str]
        Ordered list of taxon identifiers.
    w_mat : np.ndarray
        Symmetric (N, N) relative similarity matrix with values in [0.0, 1.0] (unitless).
    summary_stats : Dict[str, float]
        Summary statistics of off-diagonal SSS values (sss_mean, sss_median, etc.).
    """
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        raise ValueError(f"No FASTA records found in {fasta_path}")

    taxa_names: List[str] = [rec.id for rec in records]
    clean_seqs: List[str] = [str(rec.seq).replace("-", "").upper() for rec in records]
    if any(len(s) == 0 for s in clean_seqs):
        raise ValueError(f"Zero-length sequence detected in {fasta_path}")

    n = len(taxa_names)
    id_map = {name: i for i, name in enumerate(taxa_names)}

    mmseqs_bin = find_mmseqs_binary()
    s_bits = np.zeros((n, n), dtype=np.float64)

    with tempfile.TemporaryDirectory(prefix="mmseqs_psa_") as tmpdir:
        clean_fa = os.path.join(tmpdir, "input_clean.fasta")
        with open(clean_fa, "w") as wf:
            for name, seq in zip(taxa_names, clean_seqs):
                wf.write(f">{name}\n{seq}\n")

        db = os.path.join(tmpdir, "db")
        res = os.path.join(tmpdir, "res")
        tmp = os.path.join(tmpdir, "tmp")
        out_m8 = os.path.join(tmpdir, "out.m8")
        os.makedirs(tmp, exist_ok=True)

        subprocess.run(
            [mmseqs_bin, "createdb", clean_fa, db],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        subprocess.run(
            [
                mmseqs_bin,
                "search",
                db,
                db,
                res,
                tmp,
                "--threads",
                str(threads),
                "-e",
                str(evalue),
                "-s",
                str(sensitivity),
            ],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )
        subprocess.run(
            [mmseqs_bin, "convertalis", db, db, res, out_m8, "--format-mode", "0"],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            text=True,
        )

        if os.path.exists(out_m8):
            with open(out_m8, "r") as rf:
                for line in rf:
                    parts = line.strip().split("\t")
                    if len(parts) >= 12:
                        q, t = parts[0], parts[1]
                        if q in id_map and t in id_map:
                            try:
                                bits = float(parts[11])
                            except ValueError:
                                bits = 0.0
                            qi, ti = id_map[q], id_map[t]
                            if bits > s_bits[qi, ti]:
                                s_bits[qi, ti] = bits

    w_mat = np.zeros((n, n), dtype=np.float64)
    off_diag_scores: List[float] = []

    for i in range(n):
        w_mat[i, i] = 1.0
        for j in range(i + 1, n):
            self_score = s_bits[i, i] + s_bits[j, j]
            comp_score = s_bits[i, j] + s_bits[j, i]
            val = (comp_score / self_score) if self_score > 0.0 else 0.0
            val = min(1.0, max(0.0, float(val)))
            w_mat[i, j] = val
            w_mat[j, i] = val
            off_diag_scores.append(val)

    scores_arr = np.array(off_diag_scores, dtype=np.float64) if off_diag_scores else np.array([0.0])
    summary_stats: Dict[str, float] = {
        "sss_mean": float(np.mean(scores_arr)),
        "sss_median": float(np.median(scores_arr)),
        "sss_std": float(np.std(scores_arr)),
        "sss_min": float(np.min(scores_arr)),
        "sss_max": float(np.max(scores_arr)),
        "sss_zero_frac": float(np.mean(scores_arr <= 1e-6)),
        "num_taxa": float(n),
    }

    return taxa_names, w_mat, summary_stats


def similarity_to_distance_matrix(
    w_mat: np.ndarray,
    dist_model: str = "linear",
    epsilon: float = 1e-4,
) -> np.ndarray:
    """
    Converts a symmetrical relative sequence similarity matrix w_ij in [0, 1]
    into a pairwise distance matrix d_ij.

    Parameters
    ----------
    w_mat : np.ndarray
        (N, N) symmetrical relative sequence similarity matrix in [0.0, 1.0] (unitless).
    dist_model : str, default="linear"
        Transformation formula:
        - 'linear' / 'mmseqs_linear' / 'mmseqs_sss': d_ij = 1.0 - w_ij (unitless in [0, 1])
        - 'log' / 'mmseqs_log':                      d_ij = -ln(max(w_ij, epsilon)) (unitless)
    epsilon : float, default=1e-4
        Small positive floor for logarithmic transformation when w_ij == 0 (unitless).

    Returns
    -------
    np.ndarray
        (N, N) symmetrical pairwise distance matrix with zero diagonal (unitless).
    """
    n = w_mat.shape[0]
    d_mat = np.zeros((n, n), dtype=np.float64)
    mode = dist_model.lower().strip()
    is_linear = mode in ("linear", "mmseqs_linear", "mmseqs_sss", "sss")

    for i in range(n):
        for j in range(i + 1, n):
            w = float(w_mat[i, j])
            if is_linear:
                d_val = max(0.0, 1.0 - w)
            else:
                eff_w = max(w, epsilon)
                d_val = max(0.0, -math.log(eff_w))
            d_mat[i, j] = d_val
            d_mat[j, i] = d_val

    return d_mat

