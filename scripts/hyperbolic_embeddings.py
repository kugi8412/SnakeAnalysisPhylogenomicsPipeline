#!/usr/bin/env python3
# -*- coding: utf-8 -*- 

"""
Hyperbolic Phylogenetic Tree Embeddings

Embeds phylogenetic trees into the Poincaré ball model of hyperbolic space.
Hyperbolic geometry naturally captures hierarchical/tree-like structures because
the volume of a hyperbolic ball grows exponentially with radius (matching tree
branching), unlike Euclidean space where it grows polynomially.

The Poincaré ball model B^n = {X in R^n : ||x|| < 1} uses the metric:
    ds^2 = (2 / (1 - ||x||^2))^2 * ||dx||^2

Key properties:
- Points near the boundary represent leaves (specialized taxa)
- Points near the origin represent ancestral/root nodes
- Geodesic distances preserve tree distances

Usage:
    python hyperbolic_embeddings.py \
        --tree results/trees/supertree/sco/raw/for_consensus/ML/consensus_tree.treefile \
        --output results/embeddings/phylo_embeddings.npz \
        --dim 64 --epochs 500
"""

import os
import argparse

import numpy as np

from typing import Dict, List, Tuple, Optional


class PoincareBall:
    """Operations in the Poincaré ball model of hyperbolic space."""

    EPS = 1e-7
    MAX_NORM = 1.0 - 1e-5

    @staticmethod
    def mobius_add(x: np.ndarray, y: np.ndarray) -> np.ndarray:
        """Möbius addition in the Poincaré ball: x ⊕ y"""
        x_sq = np.sum(x ** 2)
        y_sq = np.sum(y ** 2)
        xy = np.dot(x, y)

        num = (1.0 + 2.0 * xy + y_sq) * x + (1.0 - x_sq) * y
        den = 1.0 + 2.0 * xy + x_sq * y_sq
        return num / (den + PoincareBall.EPS)

    @staticmethod
    def distance(x: np.ndarray, y: np.ndarray) -> float:
        """Hyperbolic distance in the Poincaré ball."""
        diff = x - y
        diff_sq = np.sum(diff ** 2)
        x_sq = np.sum(x ** 2)
        y_sq = np.sum(y ** 2)

        arg = 1.0 + 2.0 * diff_sq / ((1.0 - x_sq) * (1.0 - y_sq) + PoincareBall.EPS)
        return np.arccosh(np.clip(arg, 1.0 + PoincareBall.EPS, None))

    @staticmethod
    def project(x: np.ndarray) -> np.ndarray:
        """Project point back into the Poincaré ball."""
        norm = np.linalg.norm(x)
        if norm >= PoincareBall.MAX_NORM:
            return x * (PoincareBall.MAX_NORM / (norm + PoincareBall.EPS))
        return x

    @staticmethod
    def exp_map_origin(v: np.ndarray) -> np.ndarray:
        """Exponential map from origin: maps tangent vector to Poincaré ball."""
        v_norm = np.linalg.norm(v)
        if v_norm < PoincareBall.EPS:
            return v
        return np.tanh(v_norm) * (v / v_norm)

    @staticmethod
    def log_map_origin(x: np.ndarray) -> np.ndarray:
        """Logarithmic map to origin: maps point from Poincaré ball to tangent space."""
        x_norm = np.linalg.norm(x)
        if x_norm < PoincareBall.EPS:
            return x
        return np.arctanh(np.clip(x_norm, 0, PoincareBall.MAX_NORM)) * (x / x_norm)

    @staticmethod
    def riemannian_gradient(euclidean_grad: np.ndarray, x: np.ndarray) -> np.ndarray:
        """Convert Euclidean gradient to Riemannian gradient for Poincaré ball."""
        x_sq = np.sum(x ** 2)
        scale = ((1.0 - x_sq) ** 2) / 4.0
        return scale * euclidean_grad


def parse_newick_tree(tree_path: str) -> Tuple[Dict[str, List], Dict[str, float], List[str], str]:
    """
    Parse a Newick tree file and extract:
    - adjacency: node -> list of children
    - distances: (parent, child) edge weights
    - leaves: list of leaf node names
    - root: root node name
    """
    try:
        from Bio import Phylo
        from io import StringIO

        trees = list(Phylo.parse(tree_path, "newick"))
        if not trees:
            raise ValueError(f"[WARNING]: No trees found in {tree_path}")

        tree = trees[0]

        adjacency = {}
        distances = {}
        leaves = []
        node_counter = [0]

        def get_name(clade):
            if clade.name:
                return str(clade.name).strip()
            node_counter[0] += 1
            return f"_internal_{node_counter[0]}"

        def traverse(clade, parent_name=None):
            name = get_name(clade)

            if name not in adjacency:
                adjacency[name] = []

            if parent_name is not None:
                adjacency[parent_name].append(name)
                bl = clade.branch_length if clade.branch_length is not None else 0.1
                distances[(parent_name, name)] = float(bl)

            if clade.is_terminal():
                leaves.append(name)
            else:
                for child in clade.clades:
                    traverse(child, name)

        root_name = get_name(tree.root)
        traverse(tree.root)

        return adjacency, distances, leaves, root_name

    except ImportError:
        return _parse_newick_manual(tree_path)


def _parse_newick_manual(tree_path: str) -> Tuple[Dict, Dict, List, str]:
    """Minimal Newick parser fallback."""
    with open(tree_path) as f:
        nwk = f.read().strip().rstrip(';').strip()

    adjacency = {}
    distances = {}
    leaves = []
    node_id = [0]
    pos = [0]

    def new_node():
        node_id[0] += 1
        return f"_node_{node_id[0]}"

    def parse_subtree():
        if pos[0] >= len(nwk):
            return new_node(), 0.0

        if nwk[pos[0]] == '(':
            pos[0] += 1  # skip (
            name = new_node()
            adjacency[name] = []

            while True:
                child_name, child_dist = parse_subtree()
                adjacency[name].append(child_name)
                distances[(name, child_name)] = child_dist
                if child_name not in adjacency:
                    leaves.append(child_name)

                if pos[0] < len(nwk) and nwk[pos[0]] == ',':
                    pos[0] += 1
                elif pos[0] < len(nwk) and nwk[pos[0]] == ')':
                    pos[0] += 1
                    break
                else:
                    break

            # Parse optional label and distance
            label, dist = _parse_label_dist(nwk, pos)
            if label:
                old_name = name
                name = label
                adjacency[name] = adjacency.pop(old_name)
                for key in list(distances.keys()):
                    if key[0] == old_name:
                        distances[(name, key[1])] = distances.pop(key)

            return name, dist

        else:
            label, dist = _parse_label_dist(nwk, pos)
            node_name = label if label else new_node()
            adjacency[node_name] = []
            return node_name, dist

    def _parse_label_dist(s, p):
        label = ""
        while p[0] < len(s) and s[p[0]] not in ',:;()':
            label += s[p[0]]
            p[0] += 1

        dist = 0.1
        if ':' in label:
            parts = label.split(':')
            label = parts[0]
            try:
                dist = float(parts[1])
            except ValueError:
                dist = 0.1

        return label.strip(), dist

    root_name, _ = parse_subtree()
    leaves = [l for l in leaves if l in adjacency or any(l in v for v in adjacency.values())]

    return adjacency, distances, leaves, root_name


def compute_tree_distances(adjacency: Dict,
                           distances: Dict,
                           nodes: List[str]
                           ) -> Dict[Tuple[str, str], float]:
    """Compute all-pairs shortest path distances on the tree using BFS."""
    from collections import deque

    # Build undirected adjacency
    undirected = {}
    for parent, children in adjacency.items():
        if parent not in undirected:
            undirected[parent] = []
        for child in children:
            if child not in undirected:
                undirected[child] = []

            d = distances.get((parent, child), distances.get((child, parent), 0.1))
            undirected[parent].append((child, d))
            undirected[child].append((parent, d))

    pair_dists = {}
    for source in nodes:
        # BFS on tree (tree has unique paths -> NO Dijkstra)
        visited = {source: 0.0}
        queue = deque([(source, 0.0)])
        while queue:
            current, curr_dist = queue.popleft()
            for neighbor, edge_dist in undirected.get(current, []):
                if neighbor not in visited:
                    new_dist = curr_dist + edge_dist
                    visited[neighbor] = new_dist
                    queue.append((neighbor, new_dist))

        for target in nodes:
            if target != source and target in visited:
                pair_dists[(source, target)] = visited[target]

    return pair_dists


class HyperbolicEmbedder:
    """Embeds tree nodes into the Poincaré ball using Riemannian SGD."""
    def __init__(self, nodes: List[str], dim: int = 64, lr: float = 0.01, seed: int = 42):
        self.nodes = nodes
        self.dim = dim
        self.lr = lr
        self.n = len(nodes)
        self.node_to_idx = {n: i for i, n in enumerate(nodes)}

        rng = np.random.RandomState(seed)
        # Initialize near origin with small norm
        self.embeddings = rng.randn(self.n, dim) * 0.01
        for i in range(self.n):
            self.embeddings[i] = PoincareBall.project(self.embeddings[i])

    def _loss_and_grad(self, pairs: List[Tuple[int, int, float]]) -> Tuple[float, np.ndarray]:
        """Compute stress loss and analytical Riemannian gradients."""
        grad = np.zeros_like(self.embeddings)
        total_loss = 0.0

        for i, j, target_dist in pairs:
            ei = self.embeddings[i]
            ej = self.embeddings[j]

            diff = ei - ej
            diff_sq = np.sum(diff ** 2)
            ei_sq = np.sum(ei ** 2)
            ej_sq = np.sum(ej ** 2)

            denom = (1.0 - ei_sq) * (1.0 - ej_sq) + PoincareBall.EPS
            alpha = 1.0 + 2.0 * diff_sq / denom
            alpha_clamped = max(alpha, 1.0 + PoincareBall.EPS)
            hyp_dist = np.arccosh(alpha_clamped)

            residual = hyp_dist - target_dist
            total_loss += residual ** 2

            # d(arccosh(a))/da = 1/sqrt(a^2 - 1)
            darccosh = 1.0 / (np.sqrt(alpha_clamped ** 2 - 1.0) + PoincareBall.EPS)

            # d(alpha)/d(ei) = 2/(denom) * (2*diff + 2*ei * diff_sq / ((1-||ei||^2) * (1-||ej||^2)))
            factor_a = 2.0 / denom
            dalpha_dei = factor_a * (2.0 * diff + 2.0 * ei * diff_sq / ((1.0 - ei_sq) + PoincareBall.EPS))

            grad_i = 2.0 * residual * darccosh * dalpha_dei
            grad[i] += grad_i
            grad[j] -= grad_i

        return total_loss / len(pairs), grad

    def fit(self, tree_distances: Dict[Tuple[str, str], float],
            epochs: int = 500, batch_size: int = 256,
            verbose: bool = True) -> np.ndarray:
        """Train embeddings using Riemannian SGD."""
        # Build training pairs
        pairs = []
        for (n1, n2), dist in tree_distances.items():
            i = self.node_to_idx.get(n1)
            j = self.node_to_idx.get(n2)
            if i is not None and j is not None:
                pairs.append((i, j, dist))

        if not pairs:
            raise ValueError("No valid training pairs found")

        print(f"[INFO]: Training hyperbolic embeddings: {self.n} nodes, {len(pairs)} pairs, dim={self.dim}")

        best_loss = float('inf')
        best_embeddings = self.embeddings.copy()

        for epoch in range(epochs):
            np.random.shuffle(pairs)

            epoch_loss = 0.0
            n_batches = 0

            for start in range(0, len(pairs), batch_size):
                batch = pairs[start:start + batch_size]
                loss, grad = self._loss_and_grad(batch)
                epoch_loss += loss
                n_batches += 1

                # Riemannian gradient descent
                for i in range(self.n):
                    rgrad = PoincareBall.riemannian_gradient(grad[i], self.embeddings[i])
                    self.embeddings[i] -= self.lr * rgrad
                    self.embeddings[i] = PoincareBall.project(self.embeddings[i])

            avg_loss = epoch_loss / max(n_batches, 1)

            if avg_loss < best_loss:
                best_loss = avg_loss
                best_embeddings = self.embeddings.copy()

            if verbose and (epoch % 50 == 0 or epoch == epochs - 1):
                norms = np.linalg.norm(self.embeddings, axis=1)
                print(f"  Epoch {epoch:4d} | Loss: {avg_loss:.6f} | "
                      f"Norm: mean={norms.mean():.4f} max={norms.max():.4f}")

        self.embeddings = best_embeddings
        return self.embeddings

    def get_embedding(self, node_name: str) -> Optional[np.ndarray]:
        idx = self.node_to_idx.get(node_name)
        if idx is not None:
            return self.embeddings[idx].copy()
        return None

def perturb_tree_distances(tree_distances: Dict[Tuple[str, str], float],
                           noise_scale: float = 0.05,
                           seed: int = None) -> Dict[Tuple[str, str], float]:
    """
    Add multiplicative Gaussian noise to tree distances.
    This creates perturbed embeddings for evolutionary robustness training.

    d_perturbed = d * (1 + N(0, noise_scale))
    """
    rng = np.random.RandomState(seed)
    perturbed = {}
    for key, dist in tree_distances.items():
        noise = rng.normal(0, noise_scale)
        perturbed[key] = max(0.01, dist * (1.0 + noise))
    return perturbed


def generate_perturbed_embeddings(tree_path: str, dim: int = 64,
                                  n_perturbations: int = 5,
                                  noise_scale: float = 0.05,
                                  epochs: int = 300,
                                  seed: int = 42) -> Dict:
    """
    Generate multiple perturbed embeddings for robustness.
    Returns dict with 'base' and 'perturbed' embeddings.
    """
    adjacency, distances, leaves, root = parse_newick_tree(tree_path)
    all_nodes = list(adjacency.keys())
    tree_dists = compute_tree_distances(adjacency, distances, all_nodes)

    # Base embedding
    embedder = HyperbolicEmbedder(all_nodes, dim=dim, seed=seed)
    base_emb = embedder.fit(tree_dists, epochs=epochs, verbose=True)

    result = {
        'nodes': all_nodes,
        'leaves': leaves,
        'dim': dim,
        'base_embeddings': base_emb,
        'perturbed_embeddings': []
    }

    # Perturbed embeddings
    for i in range(n_perturbations):
        perturbed_dists = perturb_tree_distances(tree_dists, noise_scale=noise_scale, seed=seed + i + 1)
        p_embedder = HyperbolicEmbedder(all_nodes, dim=dim, seed=seed + i + 100)
        p_emb = p_embedder.fit(perturbed_dists, epochs=epochs // 2, verbose=False)
        result['perturbed_embeddings'].append(p_emb)
        print(f"  Perturbation {i+1}/{n_perturbations} done")

    return result


def save_embeddings(result: Dict, output_path: str) -> None:
    """Save embeddings to .npz (NumPy) format."""
    os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)

    save_dict = {
        'nodes': np.array(result['nodes'], dtype=object),
        'leaves': np.array(result['leaves'], dtype=object),
        'dim': np.array([result['dim']]),
        'base_embeddings': result['base_embeddings'],
    }

    for i, p_emb in enumerate(result['perturbed_embeddings']):
        save_dict[f'perturbed_{i}'] = p_emb

    np.savez_compressed(output_path, **save_dict)
    print(f"Saved embeddings to {output_path}")


def load_embeddings(path: str) -> Dict:
    """Load embeddings from .npz file."""
    data = np.load(path, allow_pickle=True)
    result = {
        'nodes': list(data['nodes']),
        'leaves': list(data['leaves']),
        'dim': int(data['dim'][0]),
        'base_embeddings': data['base_embeddings'],
        'perturbed_embeddings': []
    }

    i = 0
    while f'perturbed_{i}' in data:
        result['perturbed_embeddings'].append(data[f'perturbed_{i}'])
        i += 1

    return result


def get_species_embedding(result: Dict, species_id: str) -> Optional[np.ndarray]:
    """Get the hyperbolic embedding vector for a species (leaf node)."""
    try:
        idx = result['nodes'].index(species_id)
        return result['base_embeddings'][idx]
    except ValueError:
        # Try matching by prefix
        for i, node in enumerate(result['nodes']):
            if node.startswith(species_id) or species_id.startswith(node):
                return result['base_embeddings'][i]
    return None


def embeddings_to_conditioning_vector(result: Dict, target_species: List[str],
                                      aggregation: str = "mean") -> np.ndarray:
    """
    Create a phylogenetic conditioning vector Z_phylo for the generative model.

    Args:
        result: Embedding result dict
        target_species: List of species IDs to condition on
        aggregation: "mean", "weighted_mean", or "concat"

    Returns:
        Z_phylo vector for cross-attention fusion in the decoder
    """
    vectors = []
    for sp in target_species:
        emb = get_species_embedding(result, sp)
        if emb is not None:
            vectors.append(emb)

    if not vectors:
        return np.zeros(result['dim'])

    vectors = np.array(vectors)

    if aggregation == "mean":
        # Average in tangent space (log map -> mean -> exp map)
        tangent_vectors = np.array([PoincareBall.log_map_origin(v) for v in vectors])
        mean_tangent = tangent_vectors.mean(axis=0)
        return PoincareBall.exp_map_origin(mean_tangent)

    elif aggregation == "weighted_mean":
        # Weight by distance from origin (more evolved = higher weight)
        norms = np.linalg.norm(vectors, axis=1, keepdims=True)
        weights = norms / (norms.sum() + PoincareBall.EPS)
        tangent_vectors = np.array([PoincareBall.log_map_origin(v) for v in vectors])
        weighted_mean = (tangent_vectors * weights).sum(axis=0)
        return PoincareBall.exp_map_origin(weighted_mean)

    elif aggregation == "concat":
        return vectors.flatten()

    return vectors.mean(axis=0)


def compute_distortion(embeddings: np.ndarray, nodes: List[str],
                       tree_distances: Dict[Tuple[str, str], float]) -> Dict[str, float]:
    """Compute embedding quality metrics."""
    node_to_idx = {n: i for i, n in enumerate(nodes)}

    hyp_dists = []
    tree_dists = []

    for (n1, n2), td in tree_distances.items():
        i, j = node_to_idx.get(n1), node_to_idx.get(n2)
        if i is not None and j is not None:
            hd = PoincareBall.distance(embeddings[i], embeddings[j])
            hyp_dists.append(hd)
            tree_dists.append(td)

    hyp_dists = np.array(hyp_dists)
    tree_dists = np.array(tree_dists)

    # Average distortion
    ratios = hyp_dists / (tree_dists + 1e-10)
    avg_distortion = np.mean(np.abs(np.log(ratios + 1e-10)))

    # Stress (Kruskal)
    stress = np.sqrt(np.sum((hyp_dists - tree_dists) ** 2) / np.sum(tree_dists ** 2))

    # Spearman correlation (rank preservation)
    from scipy.stats import spearmanr
    corr, _ = spearmanr(hyp_dists, tree_dists)

    return {
        "avg_distortion": float(avg_distortion),
        "kruskal_stress": float(stress),
        "spearman_correlation": float(corr),
        "n_pairs": len(hyp_dists)
    }


def main():
    parser = argparse.ArgumentParser(
        description="Embed phylogenetic trees into hyperbolic (Poincare ball) space"
    )
    parser.add_argument("--tree", required=True, help="Input Newick tree file")
    parser.add_argument("--output", required=True, help="Output .npz file")
    parser.add_argument("--dim", type=int, default=64, help="Embedding dimension")
    parser.add_argument("--epochs", type=int, default=500, help="Training epochs")
    parser.add_argument("--lr", type=float, default=0.01, help="Learning rate")
    parser.add_argument("--n_perturbations", type=int, default=5, help="Number of perturbed embeddings")
    parser.add_argument("--noise_scale", type=float, default=0.05, help="Perturbation noise scale")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    parser.add_argument("--metrics", action="store_true", help="Compute and display quality metrics")

    args = parser.parse_args()

    print(f"Parsing tree: {args.tree}")
    result = generate_perturbed_embeddings(
        tree_path=args.tree,
        dim=args.dim,
        n_perturbations=args.n_perturbations,
        noise_scale=args.noise_scale,
        epochs=args.epochs,
        seed=args.seed
    )

    save_embeddings(result, args.output)

    print(f"\nEmbedding summary:")
    print(f"  Nodes: {len(result['nodes'])} ({len(result['leaves'])} leaves)")
    print(f"  Dimension: {result['dim']}")
    print(f"  Perturbations: {len(result['perturbed_embeddings'])}")

    if args.metrics:
        adjacency, distances, leaves, root = parse_newick_tree(args.tree)
        all_nodes = list(adjacency.keys())
        tree_dists = compute_tree_distances(adjacency, distances, all_nodes)
        metrics = compute_distortion(result['base_embeddings'], result['nodes'], tree_dists)
        print(f"\nQuality metrics:")
        for k, v in metrics.items():
            print(f"  {k}: {v:.4f}")


if __name__ == "__main__":
    main()
