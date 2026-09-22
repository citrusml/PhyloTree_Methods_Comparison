#!/usr/bin/env python3
"""
plot_dual_trees.py
===================
Generates publication-quality dual-tree (Tanglegram / Cophylo) comparison plots:
  - Left: Estimated phylogenetic tree with evolutionary branch lengths (Phylogram)
  - Right: NCBI Taxonomy ground-truth topology tree (Cladogram)
  - Center: Connecting bezier links matching homologous taxa/genes

Optimized for small file sizes (clean vector paths, tight layout, lightweight PNG/PDF).
"""

import re
import sys
import argparse
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any

import matplotlib.pyplot as plt
from matplotlib.path import Path as MPath
import matplotlib.patches as mpatches

SCRIPT_DIR = Path(__file__).resolve().parent
DATASET_DIR = SCRIPT_DIR.parent
RESULTS_DIR = DATASET_DIR / "results"
TREES_DIR = RESULTS_DIR / "trees"
DUAL_PLOTS_DIR = RESULTS_DIR / "dual_plots"

plt.rcParams.update({
    "font.sans-serif": ["Arial", "DejaVu Sans", "Helvetica", "sans-serif"],
    "font.size": 10,
    "axes.titlesize": 12,
    "lines.linewidth": 1.2,
})


class PhyloNode:
    def __init__(self, name: str = "", length: float = 0.0):
        self.name = name.strip("'\"")
        self.length = max(0.0, float(length))
        self.children: List["PhyloNode"] = []
        self.x: float = 0.0
        self.y: float = 0.0
        self.depth: int = 0

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def get_leaves(self) -> List["PhyloNode"]:
        if self.is_leaf:
            return [self]
        leaves = []
        for c in self.children:
            leaves.extend(c.get_leaves())
        return leaves


def parse_newick(nwk: str) -> Optional[PhyloNode]:
    """Robust, standalone Newick format parser."""
    nwk = nwk.strip()
    if not nwk:
        return None
    if nwk.endswith(";"):
        nwk = nwk[:-1].strip()

    tokens = []
    curr = []
    for char in nwk:
        if char in "(),:":
            if curr:
                tokens.append("".join(curr).strip())
                curr = []
            tokens.append(char)
        else:
            curr.append(char)
    if curr:
        tokens.append("".join(curr).strip())

    stack = [PhyloNode()]
    i = 0
    while i < len(tokens):
        token = tokens[i]
        if token == "(":
            new_node = PhyloNode()
            stack[-1].children.append(new_node)
            stack.append(new_node)
        elif token == ",":
            stack.pop()
            new_node = PhyloNode()
            stack[-1].children.append(new_node)
            stack.append(new_node)
        elif token == ")":
            stack.pop()
        elif token == ":":
            i += 1
            if i < len(tokens):
                try:
                    stack[-1].length = float(tokens[i])
                except ValueError:
                    stack[-1].length = 0.0
        else:
            stack[-1].name = token
        i += 1

    return stack[0]


def layout_tree(root: PhyloNode, is_right: bool = False, is_cladogram: bool = False):
    """
    Computes (x, y) coordinates for rectangular phylogram or cladogram.
    If is_right is True, tree faces left (root at right, leaves at left).
    If is_right is False, tree faces right (root at left, leaves at right).
    """
    leaves = root.get_leaves()
    n_leaves = len(leaves)

    # Assign y coordinates sequentially to leaves
    for idx, leaf in enumerate(leaves):
        leaf.y = float(idx)

    # Internal node y is mean of children's y
    def calc_y(node: PhyloNode):
        if not node.is_leaf:
            for c in node.children:
                calc_y(c)
            node.y = sum(c.y for c in node.children) / len(node.children)

    calc_y(root)

    # Assign depths
    def calc_depth(node: PhyloNode, current_depth: int = 0):
        node.depth = current_depth
        for c in node.children:
            calc_depth(c, current_depth + 1)

    calc_depth(root)
    max_depth = max(n.depth for n in root.get_leaves()) if n_leaves > 0 else 1

    # Calculate x coordinates
    if is_cladogram:
        # Uniform branch lengths based on depth
        def calc_x_clad(node: PhyloNode):
            node.x = float(node.depth)
            for c in node.children:
                calc_x_clad(c)
        calc_x_clad(root)
        # Normalize to [0, 1]
        for node in [root] + [n for n in root.get_leaves()]:
            pass  # will normalize below
        all_nodes = []
        def collect_nodes(n):
            all_nodes.append(n)
            for c in n.children: collect_nodes(c)
        collect_nodes(root)
        max_x = max(n.x for n in all_nodes) if all_nodes else 1.0
        if max_x == 0: max_x = 1.0
        for n in all_nodes:
            n.x = n.x / max_x
    else:
        # Distance from root
        def calc_x_dist(node: PhyloNode, current_x: float = 0.0):
            node.x = current_x + (node.length if node != root else 0.0)
            for c in node.children:
                calc_x_dist(c, node.x)
        calc_x_dist(root)
        all_nodes = []
        def collect_nodes(n):
            all_nodes.append(n)
            for c in n.children: collect_nodes(c)
        collect_nodes(root)
        max_x = max(n.x for n in all_nodes) if all_nodes else 1.0
        if max_x == 0: max_x = 1.0
        for n in all_nodes:
            n.x = n.x / max_x

    # Orient tree
    if is_right:
        # Root on right (x=1), leaves facing left (x=0)
        for n in all_nodes:
            n.x = 1.0 - n.x


def draw_tree_lines(ax, node: PhyloNode, color: str = "#334155"):
    """Draws rectangular tree branches with PhyloWeaver-style internal node dots."""
    if node.is_leaf:
        return
    # Vertical connecting line across children
    min_y = min(c.y for c in node.children)
    max_y = max(c.y for c in node.children)
    ax.plot([node.x, node.x], [min_y, max_y], color=color, lw=1.6, solid_capstyle="round")

    # Internal node circle
    ax.plot(node.x, node.y, "o", markersize=4, markerfacecolor="white", markeredgecolor=color, markeredgewidth=1.4, zorder=4)

    # Horizontal branches to children
    for c in node.children:
        ax.plot([node.x, c.x], [c.y, c.y], color=color, lw=1.6, solid_capstyle="round")
        draw_tree_lines(ax, c, color)


def plot_dual_tree(
    left_nwk: str,
    right_nwk: str,
    left_title: str = "Estimated Tree (Branch Lengths)",
    right_title: str = "NCBI Taxonomy Ground Truth (Topology)",
    hog_title: str = "HOG Comparison",
    meta_text: str = "",
    output_path: Optional[Path] = None
):
    """
    Renders the Tanglegram with Left Tree, Right Tree, and Connecting Links.
    """
    left_root = parse_newick(left_nwk)
    right_root = parse_newick(right_nwk)

    # Check if left tree actually has branch lengths
    all_left_nodes = []
    def col_l_check(n):
        all_left_nodes.append(n)
        for c in n.children: col_l_check(c)
    col_l_check(left_root)
    left_has_lengths = any(n.length > 1e-6 for n in all_left_nodes)

    layout_tree(left_root, is_right=False, is_cladogram=(not left_has_lengths))
    layout_tree(right_root, is_right=True, is_cladogram=True)

    left_leaves = {l.name: l for l in left_root.get_leaves()}
    right_leaves = {l.name: l for l in right_root.get_leaves()}

    n_leaves = max(len(left_leaves), len(right_leaves))
    fig_height = max(5.0, n_leaves * 0.22 + 1.8)
    fig, ax = plt.subplots(figsize=(12, fig_height), dpi=150)

    # Canvas mapping:
    # Left tree: X in [0.05, 0.40]
    # Right tree: X in [0.60, 0.95]
    # Center link zone: X in [0.40, 0.60]

    LEFT_START, LEFT_END = 0.05, 0.38
    RIGHT_START, RIGHT_END = 0.62, 0.95

    # Scale left tree nodes
    all_left = []
    def col_l(n):
        all_left.append(n)
        for c in n.children: col_l(c)
    col_l(left_root)
    for n in all_left:
        n.x = LEFT_START + n.x * (LEFT_END - LEFT_START)

    # Scale right tree nodes
    all_right = []
    def col_r(n):
        all_right.append(n)
        for c in n.children: col_r(c)
    col_r(right_root)
    for n in all_right:
        # Remember for right tree, 0 is leaves (left), 1 is root (right)
        n.x = RIGHT_START + (1.0 - n.x) * (RIGHT_END - RIGHT_START)

    # Draw left tree (PhyloWeaver Slate Blue)
    draw_tree_lines(ax, left_root, color="#286699")
    # Draw right tree (PhyloWeaver Dark Slate)
    draw_tree_lines(ax, right_root, color="#334155")

    # Draw Scale Bar on Left Tree (if phylogram)
    if left_has_lengths:
        max_branch = max((n.x - LEFT_START) / (LEFT_END - LEFT_START) for n in all_left)
        # Find raw max distance
        all_raw_dists = []
        def col_dists(n, d=0):
            all_raw_dists.append(d + n.length)
            for c in n.children: col_dists(c, d + n.length)
        col_dists(left_root)
        max_dist = max(all_raw_dists) if all_raw_dists else 1.0
        if max_dist <= 0: max_dist = 1.0

        scale_val = 0.1
        if max_dist < 0.15: scale_val = 0.05
        elif max_dist > 1.5: scale_val = 0.5
        scale_px_w = (scale_val / max_dist) * (LEFT_END - LEFT_START)
        scale_px_w = min(scale_px_w, (LEFT_END - LEFT_START) * 0.9)

        sb_y = -0.45
        ax.plot([LEFT_START, LEFT_START + scale_px_w], [sb_y, sb_y], color="#334155", lw=1.5)
        ax.plot([LEFT_START, LEFT_START], [sb_y - 0.12, sb_y + 0.12], color="#334155", lw=1.5)
        ax.plot([LEFT_START + scale_px_w, LEFT_START + scale_px_w], [sb_y - 0.12, sb_y + 0.12], color="#334155", lw=1.5)
        ax.text(LEFT_START + scale_px_w / 2, sb_y - 0.25, f"{scale_val} substitutions / site",
                ha="center", va="top", fontsize=8, color="#64748b", fontweight="bold")

    # Draw left leaf labels for all left leaves
    for l_name, l_node in left_leaves.items():
        ax.text(LEFT_END - 0.006, l_node.y, l_name, ha="right", va="center", fontsize=7.5, color="#1e293b", fontweight="bold")

    # Draw right leaf labels for all right leaves
    for r_name, r_node in right_leaves.items():
        ax.text(RIGHT_START + 0.006, r_node.y, r_name, ha="left", va="center", fontsize=7.5, color="#1e293b", fontweight="bold")

    # Match leaves and draw connecting links
    used_right = set()
    for l_name, l_node in left_leaves.items():
        # Clean lookup
        r_node = None
        if l_name in right_leaves and right_leaves[l_name] not in used_right:
            r_node = right_leaves[l_name]
        else:
            # Fallback prefix match (first 5 chars species code)
            for r_cand_name, cand_node in right_leaves.items():
                if cand_node not in used_right and r_cand_name[:5] == l_name[:5]:
                    r_node = cand_node
                    break

        if r_node:
            used_right.add(r_node)
            y1 = l_node.y
            y2 = r_node.y
            is_parallel = abs(y1 - y2) < 0.2
            link_color = "#3498db" if is_parallel else "#e74c3c"
            alpha = 0.6 if is_parallel else 0.85
            lw = 1.0 if is_parallel else 1.3

            # Cubic Bezier curve across middle gap
            x1 = LEFT_END
            x2 = RIGHT_START
            cx1 = x1 + (x2 - x1) * 0.45
            cx2 = x1 + (x2 - x1) * 0.55

            verts = [(x1, y1), (cx1, y1), (cx2, y2), (x2, y2)]
            codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
            path = MPath(verts, codes)
            patch = mpatches.PathPatch(path, facecolor="none", edgecolor=link_color, lw=lw, alpha=alpha)
            ax.add_patch(patch)

    # Titles & Meta Header
    fig.suptitle(hog_title, fontsize=13, fontweight="bold", y=0.98)
    if meta_text:
        ax.text(0.50, 1.01, meta_text, transform=ax.transAxes, ha="center", va="bottom",
                fontsize=9.5, color="#555555", style="italic")

    ax.text((LEFT_START + LEFT_END) / 2, -0.02, left_title, transform=ax.transAxes, ha="center", va="top",
            fontsize=11, fontweight="bold", color="#2980b9")
    ax.text((RIGHT_START + RIGHT_END) / 2, -0.02, right_title, transform=ax.transAxes, ha="center", va="top",
            fontsize=11, fontweight="bold", color="#27ae60")

    # Set bounds
    max_y = max(left_root.y, right_root.y, n_leaves)
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(-0.8, max_y + 0.8)
    ax.axis("off")

    plt.tight_layout(rect=[0, 0.03, 1, 0.95])

    if output_path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        file_size_kb = output_path.stat().st_size / 1024
        print(f"Saved: {output_path} ({file_size_kb:.1f} KB)")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot dual tree comparison (Estimated vs Ground Truth)")
    parser.add_argument("--group", type=str, default="amn", choices=["amn", "euk", "act", "bac"], help="Taxonomic group")
    parser.add_argument("--hog", type=str, default="E0732656", help="HOG identifier (e.g. E0732656)")
    parser.add_argument("--method", type=str, default="msa_ml", choices=["msa_ml", "msa_nj", "psa_nj", "gs"], help="Estimation method")
    parser.add_argument("--all", action="store_true", help="Batch generate plots for all HOGs in the group")
    parser.add_argument("--output-dir", type=str, default=str(DUAL_PLOTS_DIR), help="Output directory")
    args = parser.parse_args()

    import json
    data_json = RESULTS_DIR / "dual_tree_data.json"
    if not data_json.exists():
        print(f"Error: {data_json} not found. Run build_taxonomy_trees.py first.")
        sys.exit(1)

    with open(data_json, "r", encoding="utf-8") as f:
        data = json.load(f)

    out_dir = Path(args.output_dir)

    def plot_single(g_key: str, h_id: str, m_key: str):
        key = f"{g_key}_{h_id}"
        if key not in data["hogs"]:
            print(f"Warning: HOG {key} not found.")
            return

        hog_info = data["hogs"][key]
        gt_nwk = hog_info["ground_truth_tree"]
        est_nwk = hog_info["trees"].get(m_key, "")

        if not est_nwk:
            print(f"Warning: Tree for {m_key} on {key} is empty.")
            return

        method_name = data["methods"].get(m_key, m_key)
        group_name = data["groups"].get(g_key, g_key)

        sss = hog_info["mean_sss"]
        tcs = hog_info["tcs"].get(m_key, 0.0)
        meta_str = f"Taxa: {hog_info['num_leaves']} | Mean SSS: {sss:.4f} | Identity: {hog_info['mean_identity']}% | TCS Score: {tcs:.2f}"

        out_file = out_dir / g_key / f"{g_key}_{h_id}_{m_key}.png"
        plot_dual_tree(
            left_nwk=est_nwk,
            right_nwk=gt_nwk,
            left_title=f"{method_name} (Branch Lengths)",
            right_title="NCBI Taxonomy Ground Truth (Topology)",
            hog_title=f"{group_name} — HOG {h_id}",
            meta_text=meta_str,
            output_path=out_file
        )

    if args.all:
        print(f"Batch generating plots for group: {args.group}...")
        for key, info in data["hogs"].items():
            if info["group"] == args.group:
                plot_single(args.group, info["hog_id"], args.method)
    else:
        plot_single(args.group, args.hog, args.method)


if __name__ == "__main__":
    main()
