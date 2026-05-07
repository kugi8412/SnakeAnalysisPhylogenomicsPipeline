#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Protein Structure-Aware Phylogenetics via ESMFold API

This module adds functionality to:
1. Predict 3D protein structures using ESMFold API
2. Compute structural distances (TM-score, RMSD) between proteins
3. Build structure-based phylogenetic trees
4. Combine sequence-based and structure-based phylogenetic signals

The structural criterion for tree construction is based on:
- TM-score: global structural similarity (0-1, >0.5 = same fold)
- RMSD: root-mean-square deviation of C_alpha atoms
- GDT-TS: Global Distance Test - Total Score

Usage:
    # Predict structures
    python protein_structure_phylogeny.py predict \
        --fasta results/sco/families/family_001.fasta \
        --output_dir results/structures/

    # Build structure-based tree
    python protein_structure_phylogeny.py tree \
        --structures_dir results/structures/ \
        --output results/trees/structure_tree.newick \
        --method upgma --criterion tm_score
"""


import os
import sys
import time
import json
import argparse

import numpy as np

from Bio import Phylo
from io import StringIO
from pathlib import Path
from typing import Dict, List, Tuple, Optional


class ESMFoldClient:
    """Client for ESMFold structure prediction API."""

    API_URL = "https://api.esmatlas.com/foldSequence/v1/pdb/"
    MAX_SEQ_LENGTH = 400
    RATE_LIMIT_DELAY = 1.0  # seconds between requests

    def __init__(self, output_dir: str, max_length: int = 400):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.max_length = max_length

    def predict_structure(self, sequence: str, name: str) -> Optional[str]:
        """
        Predict 3D structure using ESMFold API.
        Returns path to PDB file, or None on failure.
        """
        import requests

        pdb_path = self.output_dir / f"{name}.pdb"

        # Skip if already predicted
        if pdb_path.exists() and pdb_path.stat().st_size > 100:
            print(f"  [{name}] Using cached structure")
            return str(pdb_path)

        # Truncate long sequences
        seq = sequence.replace('\n', '').replace(' ', '').strip()
        if len(seq) > self.max_length:
            print(f"  [{name}] Truncating sequence from {len(seq)} to {self.max_length} residues")
            seq = seq[:self.max_length]

        if len(seq) < 10:
            print(f"  [{name}] Sequence too short ({len(seq)} residues), skipping")
            return None

        try:
            response = requests.post(
                self.API_URL,
                data=seq,
                headers={"Content-Type": "text/plain"},
                timeout=120,
                verify=True
            )

            if response.status_code == 200:
                with open(pdb_path, 'w') as f:
                    f.write(response.text)
                print(f"  [{name}] Structure saved to {pdb_path}")
                time.sleep(self.RATE_LIMIT_DELAY)
                return str(pdb_path)
            else:
                print(f"  [{name}] API error: {response.status_code} - {response.text[:200]}")
                return None

        except Exception as e:
            print(f"  [{name}] Request failed: {e}")
            return None

    def predict_batch(self, sequences: Dict[str, str]) -> Dict[str, str]:
        """Predict structures for multiple sequences. Returns {name: pdb_path}."""
        results = {}
        total = len(sequences)

        for i, (name, seq) in enumerate(sequences.items(), 1):
            print(f"[{i}/{total}] Processing {name}")
            pdb_path = self.predict_structure(seq, name)
            if pdb_path:
                results[name] = pdb_path

        print(f"\nPredicted {len(results)}/{total} structures successfully")
        return results


def parse_pdb_ca_coords(pdb_path: str) -> np.ndarray:
    """Extract C_alpha coordinates from a PDB file."""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
    return np.array(coords)


def parse_pdb_plddt(pdb_path: str) -> List[float]:
    """Extract per-residue pLDDT scores from ESMFold PDB B-factor column."""
    plddts = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    bfactor = float(line[60:66])
                    plddts.append(bfactor)
                except ValueError:
                    pass
    return plddts


def _kabsch_superpose(coords1: np.ndarray, coords2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Center and superpose coords1 onto coords2 using Kabsch algorithm.
    Returns (aligned_coords1, centered_coords2)."""
    c1 = coords1 - coords1.mean(axis=0)
    c2 = coords2 - coords2.mean(axis=0)

    H = c1.T @ c2
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.eye(3)
    sign_matrix[2, 2] = np.sign(d)
    R = Vt.T @ sign_matrix @ U.T

    c1_aligned = (R @ c1.T).T
    return c1_aligned, c2


def compute_rmsd(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Compute RMSD between two sets of coordinates after optimal superposition."""
    assert coords1.shape == coords2.shape, "Coordinate arrays must have same shape"
    c1_aligned, c2 = _kabsch_superpose(coords1, coords2)
    rmsd = np.sqrt(np.mean(np.sum((c1_aligned - c2) ** 2, axis=1)))
    return float(rmsd)


def compute_tm_score(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Compute TM-score between two structures.
    TM-score is length-normalized and ranges from 0 to 1.
    >0.5 generally indicates same fold, >0.17 indicates same superfold.
    """
    L_target = len(coords2)
    if L_target < 5:
        return 0.0

    d0 = 1.24 * np.cbrt(L_target - 15) - 1.8
    d0 = max(d0, 0.5)

    L_common = min(len(coords1), len(coords2))
    c1_aligned, c2_c = _kabsch_superpose(coords1[:L_common], coords2[:L_common])

    distances = np.sqrt(np.sum((c1_aligned - c2_c) ** 2, axis=1))
    tm = np.sum(1.0 / (1.0 + (distances / d0) ** 2)) / L_target

    return float(tm)


def compute_gdt_ts(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """
    Compute GDT-TS (Global Distance Test - Total Score).
    Average percentage of residues within 1, 2, 4, 8 Angstroms
    after optimal superposition.
    """
    L_common = min(len(coords1), len(coords2))
    c1_aligned, c2_c = _kabsch_superpose(coords1[:L_common], coords2[:L_common])

    distances = np.sqrt(np.sum((c1_aligned - c2_c) ** 2, axis=1))
    n = len(distances)

    gdt = 0.0
    for threshold in [1.0, 2.0, 4.0, 8.0]:
        gdt += np.sum(distances < threshold) / n

    return float(gdt / 4.0)


CRITERIA = {
    "tm_score": {
        "func": compute_tm_score,
        "to_distance": lambda sim: max(0.0, 1.0 - sim),  # convert similarity to distance
        "description": "TM-score (higher = more similar)"
    },
    "rmsd": {
        "func": compute_rmsd,
        "to_distance": lambda d: d,  # already a distance
        "description": "RMSD in Angstroms (lower = more similar)"
    },
    "gdt_ts": {
        "func": compute_gdt_ts,
        "to_distance": lambda sim: max(0.0, 1.0 - sim),  # convert similarity to distance
        "description": "GDT-TS (higher = more similar)"
    },
    "combined": {
        "description": "Weighted combination: 0.5*TM + 0.3*GDT + 0.2*(1/(1+RMSD))"
    }
}


def compute_all_pairwise(pdb_paths: Dict[str, str],
                         criterion: str = "tm_score") -> Tuple[List[str], np.ndarray]:
    """
    Compute pairwise structural distance matrix.

    Args:
        pdb_paths: {name: pdb_file_path}
        criterion: one of "tm_score", "rmsd", "gdt_ts", "combined"

    Returns:
        (names, distance_matrix)
    """
    names = sorted(pdb_paths.keys())
    n = len(names)
    dist_matrix = np.zeros((n, n))

    # Parse all coordinates
    coords_cache = {}
    for name in names:
        try:
            coords_cache[name] = parse_pdb_ca_coords(pdb_paths[name])
        except Exception as e:
            print(f"  Warning: Could not parse {name}: {e}")

    print(f"Computing pairwise distances ({criterion}) for {n} structures...")
    total_pairs = n * (n - 1) // 2
    done = 0

    for i in range(n):
        for j in range(i + 1, n):
            n1, n2 = names[i], names[j]
            if n1 not in coords_cache or n2 not in coords_cache:
                dist_matrix[i, j] = dist_matrix[j, i] = float('inf')
                continue

            c1, c2 = coords_cache[n1], coords_cache[n2]

            if len(c1) == 0 or len(c2) == 0:
                dist_matrix[i, j] = dist_matrix[j, i] = float('inf')
                continue

            try:
                if criterion == "combined":
                    tm = compute_tm_score(c1, c2)
                    gdt = compute_gdt_ts(c1, c2)
                    rmsd = compute_rmsd(c1[:min(len(c1), len(c2))],
                                        c2[:min(len(c1), len(c2))])
                    combined_sim = 0.5 * tm + 0.3 * gdt + 0.2 * (1.0 / (1.0 + rmsd))
                    d = max(0.0, 1.0 - combined_sim)
                else:
                    raw_value = CRITERIA[criterion]["func"](c1, c2)
                    d = CRITERIA[criterion]["to_distance"](raw_value)

                dist_matrix[i, j] = d
                dist_matrix[j, i] = d
            except Exception as e:
                print(f"  Warning: {n1} vs {n2} failed: {e}")
                dist_matrix[i, j] = dist_matrix[j, i] = float('inf')

            done += 1
            if done % 50 == 0:
                print(f"  {done}/{total_pairs} pairs computed")

    return names, dist_matrix


def build_upgma_tree(names: List[str], dist_matrix: np.ndarray) -> str:
    """Build a UPGMA tree from a distance matrix. Returns Newick string."""
    n = len(names)
    clusters = {i: names[i] for i in range(n)}
    heights = {i: 0.0 for i in range(n)}
    sizes = {i: 1 for i in range(n)}

    # Use dict-of-dicts for sparse distance storage (avoids O(n²) matrix realloc)
    D = {}
    for i in range(n):
        for j in range(i + 1, n):
            D[(i, j)] = dist_matrix[i, j]

    active = set(range(n))
    next_idx = n

    while len(active) > 1:
        # Find minimum distance pair
        min_d = float('inf')
        mi, mj = -1, -1
        for (i, j), d in D.items():
            if i in active and j in active and d < min_d:
                min_d = d
                mi, mj = i, j

        if mi == -1:
            break

        new_height = min_d / 2.0
        bl_i = max(0.0, new_height - heights[mi])
        bl_j = max(0.0, new_height - heights[mj])

        new_name = f"({clusters[mi]}:{bl_i:.6f},{clusters[mj]}:{bl_j:.6f})"

        # Compute distances to new cluster
        for k in active:
            if k != mi and k != mj:
                d_ik = D.get((min(mi, k), max(mi, k)), float('inf'))
                d_jk = D.get((min(mj, k), max(mj, k)), float('inf'))
                d_new = (sizes[mi] * d_ik + sizes[mj] * d_jk) / (sizes[mi] + sizes[mj])
                D[(min(next_idx, k), max(next_idx, k))] = d_new

        # Remove old pairs
        for key in list(D.keys()):
            if mi in key or mj in key:
                del D[key]

        active.discard(mi)
        active.discard(mj)
        active.add(next_idx)

        clusters[next_idx] = new_name
        heights[next_idx] = new_height
        sizes[next_idx] = sizes[mi] + sizes[mj]
        next_idx += 1

    root_idx = active.pop()
    return clusters[root_idx] + ";"


def build_nj_tree(names: List[str], dist_matrix: np.ndarray) -> str:
    """Build a Neighbor-Joining tree from a distance matrix. Returns Newick string."""
    n = len(names)
    clusters = {i: names[i] for i in range(n)}
    active = list(range(n))
    next_idx = n

    # Use dict-of-dicts for sparse distance storage
    D = {}
    for i in range(n):
        for j in range(i + 1, n):
            D[(i, j)] = dist_matrix[i, j]

    def get_dist(a, b):
        if a == b: return 0.0
        return D.get((min(a, b), max(a, b)), float('inf'))

    while len(active) > 2:
        m = len(active)
        # Compute row sums
        row_sums = {}
        for i in active:
            row_sums[i] = sum(get_dist(i, j) for j in active if j != i)

        # Find minimum Q
        min_q = float('inf')
        mi, mj = active[0], active[1]

        for idx_a in range(m):
            for idx_b in range(idx_a + 1, m):
                i, j = active[idx_a], active[idx_b]
                q = (m - 2) * get_dist(i, j) - row_sums[i] - row_sums[j]
                if q < min_q:
                    min_q = q
                    mi, mj = i, j

        d_ij = get_dist(mi, mj)
        denom = 2.0 * (m - 2) if m > 2 else 2.0
        delta = (row_sums[mi] - row_sums[mj]) / max(denom, 1e-10)
        bl_i = max(0.0, 0.5 * d_ij + 0.5 * delta)
        bl_j = max(0.0, d_ij - bl_i)

        new_name = f"({clusters[mi]}:{bl_i:.6f},{clusters[mj]}:{bl_j:.6f})"

        # Compute distances to new node
        for k in active:
            if k != mi and k != mj:
                d_new = 0.5 * (get_dist(mi, k) + get_dist(mj, k) - d_ij)
                D[(min(next_idx, k), max(next_idx, k))] = max(0.0, d_new)

        # Remove old pairs
        for key in list(D.keys()):
            if mi in key or mj in key:
                del D[key]

        active.remove(mi)
        active.remove(mj)
        active.append(next_idx)

        clusters[next_idx] = new_name
        next_idx += 1

    # Last two
    if len(active) == 2:
        i, j = active
        bl = get_dist(i, j)
        return f"({clusters[i]}:{bl/2:.6f},{clusters[j]}:{bl/2:.6f});"
    else:
        return clusters[active[0]] + ";"


def assess_structure_quality(pdb_paths: Dict[str, str]) -> Dict[str, Dict]:
    """Assess structure prediction quality using pLDDT scores."""
    quality = {}
    for name, path in pdb_paths.items():
        plddts = parse_pdb_plddt(path)
        if plddts:
            mean_plddt = np.mean(plddts)
            quality[name] = {
                "mean_plddt": float(mean_plddt),
                "min_plddt": float(np.min(plddts)),
                "max_plddt": float(np.max(plddts)),
                "n_residues": len(plddts),
                "confident_fraction": float(np.mean(np.array(plddts) > 70)),
                "quality_tier": (
                    "high" if mean_plddt > 90 else
                    "good" if mean_plddt > 70 else
                    "low" if mean_plddt > 50 else
                    "very_low"
                )
            }
        else:
            quality[name] = {"error": "Could not parse pLDDT"}
    return quality


def read_fasta(fasta_path: str) -> Dict[str, str]:
    """Read FASTA file into {header: sequence} dict."""
    sequences = {}
    current_header = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]
                current_seq = []
            elif line:
                current_seq.append(line)

    if current_header:
        sequences[current_header] = "".join(current_seq)

    return sequences


def cmd_predict(args):
    """Predict structures using ESMFold API."""
    sequences = read_fasta(args.fasta)
    print(f"Loaded {len(sequences)} sequences from {args.fasta}")

    client = ESMFoldClient(args.output_dir, max_length=args.max_length)
    results = client.predict_batch(sequences)

    if args.assess_quality:
        quality = assess_structure_quality(results)
        quality_path = os.path.join(args.output_dir, "quality_report.json")
        with open(quality_path, 'w') as f:
            json.dump(quality, f, indent=2)
        print(f"\nQuality report saved to {quality_path}")
        for name, q in quality.items():
            if "mean_plddt" in q:
                print(f"  {name}: pLDDT={q['mean_plddt']:.1f} ({q['quality_tier']})")

    # Save manifest
    manifest_path = os.path.join(args.output_dir, "manifest.json")
    with open(manifest_path, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"Manifest saved to {manifest_path}")


def cmd_tree(args):
    """Build structure-based phylogenetic tree."""
    # Load structures
    if args.manifest:
        with open(args.manifest) as f:
            pdb_paths = json.load(f)
    else:
        pdb_paths = {}
        for f in Path(args.structures_dir).glob("*.pdb"):
            pdb_paths[f.stem] = str(f)

    if len(pdb_paths) < 3:
        print(f"Error: Need at least 3 structures, found {len(pdb_paths)}")
        sys.exit(1)

    print(f"Building tree from {len(pdb_paths)} structures using criterion: {args.criterion}")

    names, dist_matrix = compute_all_pairwise(pdb_paths, criterion=args.criterion)

    # Replace inf distances with max finite distance * 2
    finite_mask = np.isfinite(dist_matrix)
    if finite_mask.any():
        max_finite = dist_matrix[finite_mask].max()
        dist_matrix[~finite_mask] = max_finite * 2

    if args.method == "upgma":
        newick = build_upgma_tree(names, dist_matrix)
    elif args.method == "nj":
        newick = build_nj_tree(names, dist_matrix)
    else:
        print(f"Unknown method: {args.method}")
        sys.exit(1)

    os.makedirs(os.path.dirname(args.output) or '.', exist_ok=True)
    with open(args.output, 'w') as f:
        f.write(newick + "\n")

    print(f"Tree saved to {args.output}")

    # Save distance matrix
    if args.save_matrix:
        matrix_path = args.output.replace('.newick', '_distances.npz').replace('.nwk', '_distances.npz')
        np.savez_compressed(matrix_path, names=np.array(names), distances=dist_matrix)
        print(f"Distance matrix saved to {matrix_path}")


def cmd_compare(args):
    """Compare sequence-based and structure-based trees."""
    seq_tree = list(Phylo.parse(args.sequence_tree, "newick"))[0]
    struct_tree = list(Phylo.parse(args.structure_tree, "newick"))[0]
    seq_leaves = set(l.name for l in seq_tree.get_terminals())
    struct_leaves = set(l.name for l in struct_tree.get_terminals())
    common = seq_leaves & struct_leaves

    print(f"Sequence tree leaves: {len(seq_leaves)}")
    print(f"Structure tree leaves: {len(struct_leaves)}")
    print(f"Common leaves: {len(common)}")


def main():
    parser = argparse.ArgumentParser(
        description="Protein structure-aware phylogenetics via ESMFold"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    # Predict subcommand
    p_predict = subparsers.add_parser("predict", help="Predict structures using ESMFold API")
    p_predict.add_argument("--fasta", required=True, help="Input FASTA file")
    p_predict.add_argument("--output_dir", required=True, help="Output directory for PDB files")
    p_predict.add_argument("--max_length", type=int, default=400, help="Max sequence length for ESMFold")
    p_predict.add_argument("--assess_quality", action="store_true", help="Assess prediction quality (pLDDT)")

    # Tree subcommand
    p_tree = subparsers.add_parser("tree", help="Build structure-based phylogenetic tree")
    p_tree.add_argument("--structures_dir", help="Directory with PDB files")
    p_tree.add_argument("--manifest", help="JSON manifest from predict command")
    p_tree.add_argument("--output", required=True, help="Output Newick tree file")
    p_tree.add_argument("--criterion", choices=["tm_score", "rmsd", "gdt_ts", "combined"],
                        default="tm_score", help="Structural distance criterion")
    p_tree.add_argument("--method", choices=["upgma", "nj"], default="nj",
                        help="Tree construction method")
    p_tree.add_argument("--save_matrix", action="store_true", help="Save distance matrix")

    # Compare subcommand
    p_compare = subparsers.add_parser("compare", help="Compare sequence-based and structure-based trees")
    p_compare.add_argument("--sequence_tree", required=True, help="Sequence-based tree (Newick)")
    p_compare.add_argument("--structure_tree", required=True, help="Structure-based tree (Newick)")

    args = parser.parse_args()

    if args.command == "predict":
        cmd_predict(args)
    elif args.command == "tree":
        cmd_tree(args)
    elif args.command == "compare":
        cmd_compare(args)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
