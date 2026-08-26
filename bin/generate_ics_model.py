#!/usr/bin/env python3
"""
Generate Invariant Category Sites (ICS) NEXUS Model File
Location: bin/generate_ics_model.py

Extracts the empirical LG (Le and Gascuel, 2008) substitution rate matrix from IQ-TREE,
applies Dayhoff 6-category invariance masking (zeroing out substitutions between categories),
and saves the resulting model definition as a NEXUS file for use with AliSim.
"""

import sys
import os
import shutil
import argparse
import subprocess

# Standard PAML amino acid order
AA_ORDER = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]

# Dayhoff 6 amino acid physicochemical categories
DAYHOFF_GROUPS = [
    {"A", "G", "P", "S", "T"},  # Group 1: Small / Tiny
    {"D", "E", "N", "Q"},        # Group 2: Acidic / Amide
    {"H", "K", "R"},             # Group 3: Basic
    {"I", "L", "M", "V"},        # Group 4: Aliphatic / Hydrophobic
    {"F", "W", "Y"},             # Group 5: Aromatic
    {"C"}                        # Group 6: Sulfhydryl
]

def find_iqtree_cmd():
    """Locates iqtree executable (iqtree3, iqtree2, or iqtree)."""
    candidates = ["iqtree3", "iqtree2", "iqtree"]
    py_dir = os.path.dirname(sys.executable)
    for bin_name in candidates:
        full_path = os.path.join(py_dir, bin_name)
        if os.path.isfile(full_path) and os.access(full_path, os.X_OK):
            return full_path

    for bin_name in candidates:
        found = shutil.which(bin_name)
        if found:
            return found
        try:
            res = subprocess.run([bin_name, "--version"], capture_output=True, text=True)
            if res.returncode == 0:
                return bin_name
        except FileNotFoundError:
            pass
    return None

def extract_lg_from_iqtree(iqtree_cmd=None):
    """
    Extracts the empirical LG rate matrix and equilibrium frequencies from IQ-TREE binary.
    """
    if not iqtree_cmd:
        iqtree_cmd = find_iqtree_cmd()
    if not iqtree_cmd:
        raise RuntimeError("IQ-TREE executable not found. Please ensure IQ-TREE is installed.")

    res = subprocess.run(["strings", iqtree_cmd], capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Failed to inspect strings from {iqtree_cmd}")

    lines = res.stdout.splitlines()
    for i, line in enumerate(lines):
        if line.strip() == "model LG=":
            matrix_lines = lines[i+1:i+20]
            freq_line = lines[i+20].strip().rstrip(";")
            return matrix_lines, freq_line

    raise RuntimeError("Could not find built-in LG model definition inside IQ-TREE binary.")

def build_ics_nexus_content(matrix_lines, freq_line, model_name="ICS"):
    """
    Builds NEXUS model definition for ICS by zeroing out cross-category exchangeabilities.
    All amino acids within the same Dayhoff category maintain their empirical LG exchangeability.
    Cross-category substitutions are strictly set to 0.0.
    For singleton categories (Group 6: {C}), C cannot mutate to any other amino acid,
    making C an invariant site when chosen at the root of an ICS site.
    """
    aa_to_group = {aa: gid for gid, g in enumerate(DAYHOFF_GROUPS) for aa in g}

    # Parse 190 lower-triangular rates
    rates = []
    for r, l in enumerate(matrix_lines):
        vals = [float(x) for x in l.split()]
        if len(vals) != r + 1:
            raise ValueError(f"Expected {r+1} values in row {r+1}, got {len(vals)}")
        rates.append(vals)

    ics_rows = []
    for i in range(1, 20):
        row_vals = []
        for j in range(i):
            aa1, aa2 = AA_ORDER[i], AA_ORDER[j]
            g1, g2 = aa_to_group[aa1], aa_to_group[aa2]
            if g1 == g2:
                row_vals.append(f"{rates[i-1][j]:.6f}")
            else:
                row_vals.append("0.000000")
        ics_rows.append(" ".join(row_vals))

    model_body = "\n".join(ics_rows) + "\n" + freq_line + ";"

    nexus_content = f"""#nexus
begin models;
    model {model_name} = 
{model_body}
end;
"""
    return nexus_content

def main():
    parser = argparse.ArgumentParser(description="Generate Invariant Category Sites (ICS) NEXUS model file from IQ-TREE LG matrix")
    parser.add_argument("--out", default="models/ics_model.nex", help="Output NEXUS model file path (default: models/ics_model.nex)")
    parser.add_argument("--model_name", default="ICS", help="Model name in NEXUS file (default: ICS)")
    parser.add_argument("--iqtree", help="Optional path to IQ-TREE executable")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    print("Extracting empirical LG matrix from IQ-TREE...")
    m_lines, f_line = extract_lg_from_iqtree(args.iqtree)

    print(f"Applying Dayhoff 6-category invariance masking...")
    nexus_content = build_ics_nexus_content(m_lines, f_line, model_name=args.model_name)

    with open(args.out, "w") as f:
        f.write(nexus_content)

    print(f"ICS NEXUS model file successfully generated: {args.out}")

if __name__ == "__main__":
    main()
