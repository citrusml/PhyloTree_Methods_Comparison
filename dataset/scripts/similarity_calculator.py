"""
similarity_calculator.py
=========================
Core module for computing Sequence Similarity Score (SSS) and pairwise sequence identity
following the definition in Matsui & Iwasaki (Syst. Biol., 2020) and standard PSA metrics.
"""

from typing import Dict, List, Tuple, Optional, Any
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices


# Load standard BLOSUM62 substitution matrix
_BLOSUM62 = substitution_matrices.load("BLOSUM62")


def get_pairwise_aligner(
    gap_open: float = 10.0,
    gap_extend: float = 0.5,
    mode: str = "global"
) -> PairwiseAligner:
    """
    Constructs and returns a Bio.Align PairwiseAligner configured with BLOSUM62.

    Parameters
    ----------
    gap_open : float, optional
        Penalty for opening a gap (positive float, default: 10.0).
    gap_extend : float, optional
        Penalty for extending a gap (positive float, default: 0.5).
    mode : str, optional
        Alignment mode, either 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman).
        Default is 'global'.

    Returns
    -------
    PairwiseAligner
        Configured Biopython pairwise aligner object.
    """
    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = _BLOSUM62
    aligner.open_gap_score = -abs(gap_open)
    aligner.extend_gap_score = -abs(gap_extend)
    return aligner


def compute_pairwise_sss_and_identity(
    sequences: List[str],
    seq_ids: Optional[List[str]] = None,
    aligner: Optional[PairwiseAligner] = None
) -> Dict[str, Any]:
    """
    Computes pairwise Sequence Similarity Scores (SSS) w_ij and pairwise identities
    for all unique pairs (i, j) among a list of unaligned protein sequences.

    The Relative Sequence Similarity Score (SSS) is defined as:
        w_ij = max(0, S(i, j)) / mean{S(i, i), S(j, j)}
    where S(i, j) is the pairwise global alignment score with BLOSUM62,
    and S(i, i), S(j, j) are self-alignment scores.

    Parameters
    ----------
    sequences : List[str]
        List of unaligned amino acid sequences (characters '-' are removed if present).
    seq_ids : Optional[List[str]], optional
        List of sequence identifiers corresponding to `sequences`.
        If None, sequential string IDs 'seq_0', 'seq_1', ... are generated.
    aligner : Optional[PairwiseAligner], optional
        Configured aligner. If None, default global BLOSUM62 aligner is used.

    Returns
    -------
    Dict[str, Any]
        Dictionary containing:
        - 'num_sequences': int, number of sequences N
        - 'num_pairs': int, number of pairs N*(N-1)/2
        - 'mean_sss': float, average SSS w_ij across all pairs [unitless, 0.0 - 1.0]
        - 'median_sss': float, median SSS w_ij [unitless, 0.0 - 1.0]
        - 'min_sss': float, minimum SSS w_ij [unitless, 0.0 - 1.0]
        - 'max_sss': float, maximum SSS w_ij [unitless, 0.0 - 1.0]
        - 'mean_identity': float, average pairwise identity [ratio, 0.0 - 1.0]
        - 'median_identity': float, median pairwise identity [ratio, 0.0 - 1.0]
        - 'pairs': List[Dict[str, Any]], per-pair details with keys:
            ['id_i', 'id_j', 'len_i', 'len_j', 'score_raw', 'sss', 'identity', 'aligned_len']
    """
    if aligner is None:
        aligner = get_pairwise_aligner()

    n = len(sequences)
    if seq_ids is None:
        seq_ids = [f"seq_{k}" for k in range(n)]

    # Clean sequences: strip '-' gaps and whitespace, convert to upper
    clean_seqs = [s.replace("-", "").strip().upper() for s in sequences]

    if n < 2:
        return {
            "num_sequences": n,
            "num_pairs": 0,
            "mean_sss": 1.0 if n == 1 else 0.0,
            "median_sss": 1.0 if n == 1 else 0.0,
            "min_sss": 1.0 if n == 1 else 0.0,
            "max_sss": 1.0 if n == 1 else 0.0,
            "std_sss": 0.0,
            "mean_identity": 1.0 if n == 1 else 0.0,
            "median_identity": 1.0 if n == 1 else 0.0,
            "min_identity": 1.0 if n == 1 else 0.0,
            "max_identity": 1.0 if n == 1 else 0.0,
            "std_identity": 0.0,
            "pairs": []
        }

    # Pre-calculate self-alignment scores S(i, i)
    self_scores = [float(aligner.score(s, s)) for s in clean_seqs]

    pairs_data: List[Dict[str, Any]] = []
    sss_list: List[float] = []
    identities_list: List[float] = []

    for i in range(n):
        s1 = clean_seqs[i]
        len1 = len(s1)
        self_i = self_scores[i]

        for j in range(i + 1, n):
            s2 = clean_seqs[j]
            len2 = len(s2)
            self_j = self_scores[j]

            # Raw alignment score S(i, j)
            raw_score = float(aligner.score(s1, s2))

            # Symmetric relative score w_ij
            mean_self = (self_i + self_j) / 2.0
            w_ij = max(0.0, raw_score) / mean_self if mean_self > 0 else 0.0
            w_ij = min(1.0, w_ij)
            sss_list.append(w_ij)

            # Extract alignment for pairwise sequence identity
            alns = aligner.align(s1, s2)
            top_aln = alns[0]
            a1, a2 = top_aln[0], top_aln[1]

            # Calculate identity on non-gap aligned positions
            matches = 0
            valid_len = 0
            for c1, c2 in zip(a1, a2):
                if c1 != "-" and c2 != "-":
                    valid_len += 1
                    if c1 == c2:
                        matches += 1

            identity = matches / valid_len if valid_len > 0 else 0.0
            identities_list.append(identity)

            pairs_data.append({
                "id_i": seq_ids[i],
                "id_j": seq_ids[j],
                "len_i": len1,
                "len_j": len2,
                "score_raw": raw_score,
                "sss": w_ij,
                "identity": identity,
                "aligned_len": valid_len
            })

    return {
        "num_sequences": n,
        "num_pairs": len(pairs_data),
        "mean_sss": float(np.mean(sss_list)),
        "median_sss": float(np.median(sss_list)),
        "min_sss": float(np.min(sss_list)),
        "max_sss": float(np.max(sss_list)),
        "std_sss": float(np.std(sss_list)),
        "mean_identity": float(np.mean(identities_list)),
        "median_identity": float(np.median(identities_list)),
        "min_identity": float(np.min(identities_list)),
        "max_identity": float(np.max(identities_list)),
        "std_identity": float(np.std(identities_list)),
        "pairs": pairs_data
    }
