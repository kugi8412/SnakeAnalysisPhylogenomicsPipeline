#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Comprehensive pytest suite for the SAAP pipeline.

Tests cover:
- select_orthologs.py: parse_header, process_cluster_job (SCO + paralogs)
- hyperbolic_embeddings.py: PoincareBall math, tree parsing, embedding
- protein_structure_phylogeny.py: Kabsch, TM-score, RMSD, GDT-TS, tree builders
- filter_trees.py / merge_trees.py: bootstrap logic
- concat_matrix.py: supermatrix assembly

Run:
    cd SnakeAnalysisPhylogenomicsPipeline-main
    pip install pytest biopython numpy scipy
    pytest tests/ -v
"""

import os
import sys
import pytest
import tempfile
import textwrap

import numpy as np

from unittest.mock import patch, MagicMock

# Add scripts/ to path so we can import pipeline modules
SCRIPTS_DIR = os.path.join(os.path.dirname(__file__), "..", "scripts")
sys.path.insert(0, SCRIPTS_DIR)

# WSL on windows
def _write_temp(suffix: str, content: str) -> str:
    """Write content to a temp file and return path (closed, safe to re-open on Windows)."""
    fd, path = tempfile.mkstemp(suffix=suffix)
    with os.fdopen(fd, 'w') as f:
        f.write(content)
    return path


class TestParseHeader:
    """Tests for select_orthologs.parse_header"""

    def _import(self):
        from select_orthologs import parse_header
        return parse_header

    def test_standard_header(self):
        parse_header = self._import()
        sp, full = parse_header("G01_WP_12345")
        assert sp == "G01"
        assert full == "G01_WP_12345"

    def test_single_underscore(self):
        parse_header = self._import()
        sp, full = parse_header("G25_protein123")
        assert sp == "G25"
        assert full == "G25_protein123"

    def test_no_underscore(self):
        parse_header = self._import()
        sp, full = parse_header("G01")
        assert sp == "G01"
        assert full == "G01"

    def test_empty_string(self):
        parse_header = self._import()
        sp, full = parse_header("")
        assert sp == "UNKNOWN"
        assert full == ""

    def test_multiple_underscores(self):
        parse_header = self._import()
        sp, full = parse_header("G02_WP_001_extra_info")
        assert sp == "G02"
        assert full == "G02_WP_001_extra_info"

    def test_leading_underscore(self):
        parse_header = self._import()
        sp, full = parse_header("_orphan")
        assert sp == ""
        assert full == "_orphan"


class TestProcessClusterJob:
    """Tests for select_orthologs.process_cluster_job"""

    def setup_method(self):
        """Set up a mock SEQ_DB for all tests."""
        import select_orthologs
        select_orthologs.SEQ_DB = {
            "G01_prot1": "MKFLILLFNILCALGSTDYA",
            "G02_prot1": "MKFLILLFNILCALGSTDYB",
            "G02_prot2": "MKFLLLLFNILCALGSTDYB",
            "G03_prot1": "MKFLILLFNILCALGSTDYC",
            "G04_prot1": "MKFLILLFNILCALGSTDYD",
            "G04_prot2": "MKFLILLFNILCALGSTDYE",
            "G04_prot3": "MKFLILLFNILCALGSTDYF",
            "G30_prot1": "MKFLILLFNILCALGSTDYZ",
        }

    def test_singleton_cluster(self):
        from select_orthologs import process_cluster_job
        members = {"G01": ["G01_prot1"]}
        status, data = process_cluster_job(("rep1", members, 2, {"G30"}, "sco"))
        assert status == "SINGLETON"
        assert data is None

    def test_low_species_coverage(self):
        from select_orthologs import process_cluster_job
        members = {
            "G01": ["G01_prot1"],
            "G02": ["G02_prot1"],
        }
        # min_species=3 but only G01 and G02 (1 ingroup after excluding G30)
        # Actually G01 is not in outgroups here, so both are ingroup = 2 < 3
        status, data = process_cluster_job(("rep1", members, 3, {"G30"}, "sco"))
        assert status == "LOW_SPECIES"

    def test_low_species_outgroup_not_counted(self):
        from select_orthologs import process_cluster_job
        members = {
            "G01": ["G01_prot1"],  # outgroup
            "G02": ["G02_prot1"],
            "G30": ["G30_prot1"],  # outgroup
        }
        # Only G02 is ingroup -> 1 < min_species=2
        status, data = process_cluster_job(("rep1", members, 2, {"G01", "G30"}, "sco"))
        assert status == "LOW_SPECIES"

    def test_sco_all_single_copy(self):
        from select_orthologs import process_cluster_job
        members = {
            "G01": ["G01_prot1"],
            "G02": ["G02_prot1"],
            "G03": ["G03_prot1"],
        }
        status, data = process_cluster_job(("rep1", members, 2, {"G30"}, "sco"))
        assert status == "OK"
        rep, seqs = data
        assert rep == "rep1"

        species_in_result = {h for h, _ in seqs}
        assert species_in_result == {"G01", "G02", "G03"}

        seq_dict = {h: pid for h, pid in seqs}
        assert seq_dict["G01"] == "G01_prot1"
        assert seq_dict["G02"] == "G02_prot1"
        assert seq_dict["G03"] == "G03_prot1"

    def test_sco_no_anchor_all_multicopy(self):
        from select_orthologs import process_cluster_job
        members = {
            "G02": ["G02_prot1", "G02_prot2"],
            "G04": ["G04_prot1", "G04_prot2"],
        }
        status, data = process_cluster_job(("rep1", members, 2, set(), "sco"))
        assert status == "NO_ANCHOR"

    def test_sco_with_multicopy_species(self):
        """SCO mode: species with >1 copy should get one representative."""
        from select_orthologs import process_cluster_job
        members = {
            "G01": ["G01_prot1"],       # single copy -> anchor
            "G03": ["G03_prot1"],       # single copy -> anchor
            "G04": ["G04_prot1", "G04_prot2", "G04_prot3"]
        }

        mock_result = MagicMock()
        mock_result.stdout = "G04_prot1\t200\n"
        with patch("select_orthologs.subprocess.run", return_value=mock_result):
            status, data = process_cluster_job(("rep1", members, 2, set(), "sco"))
        assert status == "OK"
        rep, seqs = data
        species_in_result = {h for h, _ in seqs}
        assert "G01" in species_in_result
        assert "G03" in species_in_result
        assert "G04" in species_in_result
        # G04 should have exactly 1 representative
        g04_seqs = [(h, p) for h, p in seqs if h == "G04"]
        assert len(g04_seqs) == 1
        assert g04_seqs[0][1] in ["G04_prot1", "G04_prot2", "G04_prot3"]

    def test_paralogs_mode(self):
        from select_orthologs import process_cluster_job
        members = {
            "G01": ["G01_prot1"],
            "G02": ["G02_prot1", "G02_prot2"],
        }
        status, data = process_cluster_job(("rep1", members, 1, set(), "paralogs"))
        assert status == "OK"
        _, seqs = data
        headers = [h for h, _ in seqs]
        # G01 has 1 copy
        assert "G01" in headers
        # G02 has 2 copies -> _p1 and _p2
        assert "G02_p1" in headers
        assert "G02_p2" in headers
        assert len(seqs) == 3

    def test_paralogs_mode_single_copy_no_suffix(self):
        from select_orthologs import process_cluster_job
        members = {
            "G01": ["G01_prot1"],
            "G03": ["G03_prot1"],
        }
        status, data = process_cluster_job(("rep1", members, 1, set(), "paralogs"))
        assert status == "OK"
        headers = [h for h, _ in data[1]]
        assert "G01" in headers
        assert "G03" in headers
        # No _p suffixes for single-copy
        assert all("_p" not in h for h in headers)

    def test_missing_seq_db_entry_filtered(self):
        """Proteins not in SEQ_DB should be filtered out."""
        import select_orthologs
        from select_orthologs import process_cluster_job
        # Temporarily remove an entry
        saved = select_orthologs.SEQ_DB.pop("G03_prot1", None)
        try:
            members = {
                "G01": ["G01_prot1"],
                "G03": ["G03_prot1"],  # not in SEQ_DB
            }
            status, data = process_cluster_job(("rep1", members, 2, set(), "sco"))
            # Only 1 valid species → LOW_SPECIES
            assert status in ("LOW_SPECIES", "SINGLETON")
        finally:
            if saved:
                select_orthologs.SEQ_DB["G03_prot1"] = saved

    def test_duplicate_members_in_cluster(self):
        """Same protein appearing twice should not double-count."""
        from select_orthologs import process_cluster_job
        members = {
            "G01": ["G01_prot1", "G01_prot1"],  # duplicate
            "G02": ["G02_prot1"],
        }
        # G01 has 2 entries (both same protein) -> treated as multi-copy, needs BLAST
        mock_result = MagicMock()
        mock_result.stdout = "G01_prot1\t200\n"
        with patch("select_orthologs.subprocess.run", return_value=mock_result):
            status, data = process_cluster_job(("rep1", members, 2, set(), "sco"))
        assert status == "OK"
        _, seqs = data
        # G01 should have exactly 1 representative
        g01_entries = [(h, p) for h, p in seqs if h == "G01"]
        assert len(g01_entries) == 1

    def test_outgroup_species_included_in_output(self):
        """Outgroup species should still appear in final sequences."""
        from select_orthologs import process_cluster_job
        members = {
            "G01": ["G01_prot1"],  # outgroup
            "G02": ["G02_prot1"],
            "G03": ["G03_prot1"],
            "G30": ["G30_prot1"],  # outgroup
        }
        status, data = process_cluster_job(("rep1", members, 2, {"G01", "G30"}, "sco"))
        assert status == "OK"
        _, seqs = data
        species = {h for h, _ in seqs}
        assert "G01" in species  # outgroup included
        assert "G30" in species  # outgroup included
        assert "G02" in species
        assert "G03" in species


class TestPoincareBall:
    """Tests for PoincareBall operations."""

    def _import(self):
        from hyperbolic_embeddings import PoincareBall
        return PoincareBall

    def test_distance_same_point(self):
        PB = self._import()
        x = np.array([0.1, 0.2, 0.0])
        d = PB.distance(x, x)
        assert abs(d) < 1e-3  # arccosh(1+eps) = sqrt(2*eps)

    def test_distance_symmetry(self):
        PB = self._import()
        x = np.array([0.1, 0.2])
        y = np.array([0.3, -0.1])
        assert abs(PB.distance(x, y) - PB.distance(y, x)) < 1e-10

    def test_distance_positive(self):
        PB = self._import()
        x = np.array([0.1, 0.2])
        y = np.array([0.3, -0.1])
        assert PB.distance(x, y) > 0

    def test_distance_increases_near_boundary(self):
        """Points near boundary should have larger distances."""
        PB = self._import()
        origin = np.array([0.0, 0.0])
        near = np.array([0.1, 0.0])
        far = np.array([0.9, 0.0])
        d_near = PB.distance(origin, near)
        d_far = PB.distance(origin, far)
        assert d_far > d_near

    def test_project_inside_ball(self):
        PB = self._import()
        x = np.array([0.3, 0.4])
        result = PB.project(x)
        assert np.linalg.norm(result) < 1.0

    def test_project_outside_ball(self):
        PB = self._import()
        x = np.array([5.0, 5.0])
        result = PB.project(x)
        assert np.linalg.norm(result) < 1.0

    def test_exp_map_log_map_roundtrip(self):
        PB = self._import()
        v = np.array([0.3, -0.2, 0.1])
        x = PB.exp_map_origin(v)
        assert np.linalg.norm(x) < 1.0
        v_recovered = PB.log_map_origin(x)
        np.testing.assert_allclose(v, v_recovered, atol=1e-6)

    def test_exp_map_zero(self):
        PB = self._import()
        v = np.array([0.0, 0.0])
        x = PB.exp_map_origin(v)
        np.testing.assert_allclose(x, np.array([0.0, 0.0]), atol=1e-7)

    def test_mobius_add_identity(self):
        PB = self._import()
        x = np.array([0.3, 0.2])
        zero = np.array([0.0, 0.0])
        result = PB.mobius_add(x, zero)
        np.testing.assert_allclose(result, x, atol=1e-6)

    def test_riemannian_gradient_at_origin(self):
        PB = self._import()
        grad = np.array([1.0, 2.0])
        x = np.array([0.0, 0.0])
        rgrad = PB.riemannian_gradient(grad, x)
        # At origin: scale = (1-0)^2 / 4 = 0.25
        np.testing.assert_allclose(rgrad, 0.25 * grad, atol=1e-10)


class TestTreeParsing:
    """Tests for Newick tree parsing."""

    def test_parse_simple_tree(self):
        from hyperbolic_embeddings import parse_newick_tree
        path = _write_temp('.nwk', "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);\n")
        try:
            adj, dists, leaves, root = parse_newick_tree(path)
            assert set(leaves) == {"A", "B", "C", "D"}
            assert len(adj) > 0
        finally:
            os.unlink(path)

    def test_compute_tree_distances(self):
        from hyperbolic_embeddings import parse_newick_tree, compute_tree_distances
        path = _write_temp('.nwk', "(A:1.0,B:2.0);\n")
        try:
            adj, dists, leaves, root = parse_newick_tree(path)
            all_nodes = list(adj.keys())
            tree_dists = compute_tree_distances(adj, dists, all_nodes)
            # Distance A<=>B should be 1.0 + 2.0 = 3.0
            d_ab = tree_dists.get(("A", "B"), tree_dists.get(("B", "A"), None))
            assert d_ab is not None
            assert abs(d_ab - 3.0) < 1e-6
        finally:
            os.unlink(path)


class TestHyperbolicEmbedder:
    """Tests for the embedding optimization."""

    def test_fit_small_tree(self):
        from hyperbolic_embeddings import HyperbolicEmbedder
        nodes = ["A", "B", "C"]
        tree_dists = {
            ("A", "B"): 1.0, ("B", "A"): 1.0,
            ("A", "C"): 2.0, ("C", "A"): 2.0,
            ("B", "C"): 1.5, ("C", "B"): 1.5,
        }
        embedder = HyperbolicEmbedder(nodes, dim=4, lr=0.05, seed=42)
        emb = embedder.fit(tree_dists, epochs=100, verbose=False)
        assert emb.shape == (3, 4)
        # All inside Poincare ball
        norms = np.linalg.norm(emb, axis=1)
        assert norms.max() < 1.0

    def test_get_embedding(self):
        from hyperbolic_embeddings import HyperbolicEmbedder
        nodes = ["X", "Y"]
        embedder = HyperbolicEmbedder(nodes, dim=2, seed=1)
        emb = embedder.get_embedding("X")
        assert emb is not None
        assert len(emb) == 2
        assert embedder.get_embedding("Z") is None


class TestPerturbation:
    def test_perturb_distances(self):
        from hyperbolic_embeddings import perturb_tree_distances
        original = {("A", "B"): 1.0, ("B", "C"): 2.0}
        perturbed = perturb_tree_distances(original, noise_scale=0.1, seed=42)
        assert set(perturbed.keys()) == set(original.keys())
        # Values should differ but stay positive
        for key in original:
            assert perturbed[key] > 0
            assert perturbed[key] != original[key]


class TestKabschSuperposition:
    def test_identical_coords(self):
        from protein_structure_phylogeny import _kabsch_superpose
        coords = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        aligned, ref = _kabsch_superpose(coords, coords)
        rmsd = np.sqrt(np.mean(np.sum((aligned - ref) ** 2, axis=1)))
        assert rmsd < 1e-10

    def test_translated_coords(self):
        from protein_structure_phylogeny import _kabsch_superpose
        coords1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
        coords2 = coords1 + 10.0  # translated
        aligned, ref = _kabsch_superpose(coords1, coords2)
        rmsd = np.sqrt(np.mean(np.sum((aligned - ref) ** 2, axis=1)))
        assert rmsd < 1e-6


class TestStructuralMetrics:
    @pytest.fixture
    def identical_coords(self):
        return np.array([
            [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0], [4.0, 0.0, 0.0], [5.0, 0.0, 0.0],
            [6.0, 0.0, 0.0], [7.0, 0.0, 0.0], [8.0, 0.0, 0.0],
            [9.0, 0.0, 0.0],
        ])

    def test_rmsd_identical(self, identical_coords):
        from protein_structure_phylogeny import compute_rmsd
        rmsd = compute_rmsd(identical_coords, identical_coords)
        assert rmsd < 1e-10

    def test_rmsd_positive(self, identical_coords):
        from protein_structure_phylogeny import compute_rmsd
        perturbed = identical_coords.copy()
        # Non-uniform perturbation
        perturbed[:, 1] += np.arange(len(perturbed)) * 0.5
        rmsd = compute_rmsd(identical_coords, perturbed)
        assert rmsd > 0

    def test_tm_score_identical(self, identical_coords):
        from protein_structure_phylogeny import compute_tm_score
        tm = compute_tm_score(identical_coords, identical_coords)
        assert tm > 0.99

    def test_tm_score_range(self, identical_coords):
        from protein_structure_phylogeny import compute_tm_score
        shifted = identical_coords.copy()
        np.random.seed(42)
        shifted[:, 0] += np.random.randn(10) * 5
        tm = compute_tm_score(identical_coords, shifted)
        assert 0.0 <= tm <= 1.0

    def test_tm_score_short_sequence(self):
        from protein_structure_phylogeny import compute_tm_score
        coords = np.array([[0, 0, 0], [1, 0, 0]], dtype=float)
        tm = compute_tm_score(coords, coords)
        assert tm == 0.0  # too short

    def test_gdt_ts_identical(self, identical_coords):
        from protein_structure_phylogeny import compute_gdt_ts
        gdt = compute_gdt_ts(identical_coords, identical_coords)
        assert gdt > 0.99

    def test_gdt_ts_range(self, identical_coords):
        from protein_structure_phylogeny import compute_gdt_ts
        shifted = identical_coords.copy()
        shifted[:, 0] += np.random.randn(10) * 10
        gdt = compute_gdt_ts(identical_coords, shifted)
        assert 0.0 <= gdt <= 1.0


class TestUPGMATree:
    def test_two_taxa(self):
        from protein_structure_phylogeny import build_upgma_tree
        names = ["A", "B"]
        dist = np.array([[0.0, 2.0], [2.0, 0.0]])
        newick = build_upgma_tree(names, dist)
        assert newick.endswith(";")
        assert "A" in newick
        assert "B" in newick

    def test_three_taxa(self):
        from protein_structure_phylogeny import build_upgma_tree
        names = ["A", "B", "C"]
        dist = np.array([
            [0.0, 1.0, 3.0],
            [1.0, 0.0, 3.0],
            [3.0, 3.0, 0.0],
        ])
        newick = build_upgma_tree(names, dist)
        assert newick.endswith(";")
        # A and B are closest -> should cluster first
        assert "A" in newick and "B" in newick and "C" in newick

    def test_result_parseable(self):
        from protein_structure_phylogeny import build_upgma_tree
        from Bio import Phylo
        from io import StringIO
        names = ["X", "Y", "Z"]
        dist = np.array([
            [0.0, 2.0, 4.0],
            [2.0, 0.0, 3.0],
            [4.0, 3.0, 0.0],
        ])
        newick = build_upgma_tree(names, dist)
        tree = Phylo.read(StringIO(newick), "newick")
        leaves = sorted(l.name for l in tree.get_terminals())
        assert leaves == ["X", "Y", "Z"]


class TestNJTree:
    def test_three_taxa(self):
        from protein_structure_phylogeny import build_nj_tree
        names = ["A", "B", "C"]
        dist = np.array([
            [0.0, 1.0, 3.0],
            [1.0, 0.0, 2.5],
            [3.0, 2.5, 0.0],
        ])
        newick = build_nj_tree(names, dist)
        assert newick.endswith(";")
        assert "A" in newick and "B" in newick and "C" in newick

    def test_four_taxa(self):
        from protein_structure_phylogeny import build_nj_tree
        from Bio import Phylo
        from io import StringIO
        names = ["A", "B", "C", "D"]
        dist = np.array([
            [0, 5, 9, 9],
            [5, 0, 10, 10],
            [9, 10, 0, 8],
            [9, 10, 8, 0],
        ], dtype=float)
        newick = build_nj_tree(names, dist)
        tree = Phylo.read(StringIO(newick), "newick")
        leaves = sorted(l.name for l in tree.get_terminals())
        assert leaves == ["A", "B", "C", "D"]

    def test_two_taxa(self):
        from protein_structure_phylogeny import build_nj_tree
        names = ["A", "B"]
        dist = np.array([[0.0, 4.0], [4.0, 0.0]])
        newick = build_nj_tree(names, dist)
        assert "A" in newick and "B" in newick
        assert newick.endswith(";")


class TestPDBParser:
    def test_parse_ca_coords(self):
        from protein_structure_phylogeny import parse_pdb_ca_coords
        pdb_content = textwrap.dedent("""\
            ATOM      1  N   ALA A   1       1.000   2.000   3.000  1.00 80.00           N
            ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 85.00           C
            ATOM      3  C   ALA A   1       3.000   4.000   5.000  1.00 82.00           C
            ATOM      4  N   GLY A   2       4.000   5.000   6.000  1.00 78.00           N
            ATOM      5  CA  GLY A   2       5.000   6.000   7.000  1.00 90.00           C
            END
        """)
        path = _write_temp('.pdb', pdb_content)
        try:
            coords = parse_pdb_ca_coords(path)
            assert coords.shape == (2, 3)
            np.testing.assert_allclose(coords[0], [2.0, 3.0, 4.0])
            np.testing.assert_allclose(coords[1], [5.0, 6.0, 7.0])
        finally:
            os.unlink(path)

    def test_parse_plddt(self):
        from protein_structure_phylogeny import parse_pdb_plddt
        pdb_content = textwrap.dedent("""\
            ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 85.50           C
            ATOM      5  CA  GLY A   2       5.000   6.000   7.000  1.00 92.30           C
            END
        """)
        path = _write_temp('.pdb', pdb_content)
        try:
            plddts = parse_pdb_plddt(path)
            assert len(plddts) == 2
            assert abs(plddts[0] - 85.50) < 0.01
            assert abs(plddts[1] - 92.30) < 0.01
        finally:
            os.unlink(path)


class TestQualityAssessment:
    def test_assess_quality_tiers(self):
        from protein_structure_phylogeny import assess_structure_quality
        pdb_high = textwrap.dedent("""\
            ATOM      2  CA  ALA A   1       2.000   3.000   4.000  1.00 95.00           C
            ATOM      5  CA  GLY A   2       5.000   6.000   7.000  1.00 92.00           C
            END
        """)
        path = _write_temp('.pdb', pdb_high)
        try:
            quality = assess_structure_quality({"test": path})
            assert quality["test"]["quality_tier"] == "high"
            assert quality["test"]["mean_plddt"] > 90
            assert quality["test"]["n_residues"] == 2
        finally:
            os.unlink(path)


class TestFastaIO:
    def test_read_fasta(self):
        from protein_structure_phylogeny import read_fasta
        content = ">prot1\nMKFLI\nLLFNI\n>prot2\nACDEF\n"
        path = _write_temp('.fasta', content)
        try:
            seqs = read_fasta(path)
            assert len(seqs) == 2
            assert seqs["prot1"] == "MKFLILLFNI"
            assert seqs["prot2"] == "ACDEF"
        finally:
            os.unlink(path)


class TestDistortionMetrics:
    def test_compute_distortion(self):
        from hyperbolic_embeddings import compute_distortion, PoincareBall
        nodes = ["A", "B", "C"]
        # Place points manually
        emb = np.array([
            [0.1, 0.0],
            [0.5, 0.0],
            [-0.3, 0.2],
        ])
        tree_dists = {
            ("A", "B"): PoincareBall.distance(emb[0], emb[1]),
            ("B", "C"): PoincareBall.distance(emb[1], emb[2]),
            ("A", "C"): PoincareBall.distance(emb[0], emb[2]),
        }
        metrics = compute_distortion(emb, nodes, tree_dists)
        # Perfect embedding → distortion should be ~0
        assert metrics["kruskal_stress"] < 0.01
        assert metrics["spearman_correlation"] > 0.99
        assert metrics["n_pairs"] == 3


class TestConditioningVector:
    def test_mean_aggregation(self):
        from hyperbolic_embeddings import embeddings_to_conditioning_vector, PoincareBall
        result = {
            'nodes': ['A', 'B', 'C'],
            'dim': 4,
            'base_embeddings': np.array([
                [0.1, 0.0, 0.0, 0.0],
                [0.0, 0.1, 0.0, 0.0],
                [0.0, 0.0, 0.1, 0.0],
            ])
        }
        vec = embeddings_to_conditioning_vector(result, ['A', 'B'], aggregation='mean')
        assert len(vec) == 4
        assert np.linalg.norm(vec) < 1.0  # inside ball

    def test_missing_species(self):
        from hyperbolic_embeddings import embeddings_to_conditioning_vector
        result = {
            'nodes': ['A', 'B'],
            'dim': 2,
            'base_embeddings': np.array([[0.1, 0.0], [0.0, 0.1]])
        }
        vec = embeddings_to_conditioning_vector(result, ['MISSING'], aggregation='mean')
        np.testing.assert_allclose(vec, [0.0, 0.0], atol=1e-10)


class TestEndToEndEmbedding:
    def test_embed_small_newick(self):
        from hyperbolic_embeddings import (
            parse_newick_tree, compute_tree_distances,
            HyperbolicEmbedder, PoincareBall
        )
        path = _write_temp('.nwk', "((A:0.5,B:0.3):0.2,(C:0.6,D:0.4):0.1);\n")
        try:
            adj, dists, leaves, root = parse_newick_tree(path)
            all_nodes = list(adj.keys())
            tree_dists = compute_tree_distances(adj, dists, all_nodes)

            # Only embed leaves
            leaf_dists = {k: v for k, v in tree_dists.items()
                          if k[0] in leaves and k[1] in leaves}
            embedder = HyperbolicEmbedder(leaves, dim=8, lr=0.05, seed=0)
            emb = embedder.fit(leaf_dists, epochs=200, verbose=False)

            assert emb.shape == (4, 8)
            norms = np.linalg.norm(emb, axis=1)
            assert norms.max() < 1.0

            # A and B should be closer than A and D
            d_ab = PoincareBall.distance(emb[0], emb[1])
            d_ad = PoincareBall.distance(emb[0], emb[3])
            assert d_ab < d_ad
        finally:
            os.unlink(path)

# EDGE CASE TEST

class TestEdgeCases:
    def test_single_protein_cluster_sco(self):
        """Cluster with single protein total should be SINGLETON."""
        import select_orthologs
        select_orthologs.SEQ_DB = {"G01_p1": "MKFL"}
        from select_orthologs import process_cluster_job
        status, _ = process_cluster_job(("rep", {"G01": ["G01_p1"]}, 1, set(), "sco"))
        assert status == "SINGLETON"

    def test_single_protein_cluster_paralogs(self):
        import select_orthologs
        select_orthologs.SEQ_DB = {"G01_p1": "MKFL"}
        from select_orthologs import process_cluster_job
        status, _ = process_cluster_job(("rep", {"G01": ["G01_p1"]}, 1, set(), "paralogs"))
        assert status == "SINGLETON"

    def test_all_outgroup(self):
        """Cluster with only outgroup species should be LOW_SPECIES."""
        import select_orthologs
        select_orthologs.SEQ_DB = {"G01_p1": "MKFL", "G30_p1": "MKFL"}
        from select_orthologs import process_cluster_job
        members = {"G01": ["G01_p1"], "G30": ["G30_p1"]}
        status, _ = process_cluster_job(("rep", members, 1, {"G01", "G30"}, "sco"))
        assert status == "LOW_SPECIES"

    def test_empty_newick_file(self):
        from hyperbolic_embeddings import parse_newick_tree
        path = _write_temp('.nwk', "")
        try:
            with pytest.raises(Exception):
                parse_newick_tree(path)
        finally:
            os.unlink(path)

    def test_upgma_single_taxon(self):
        """UPGMA with single taxon should return just the name."""
        from protein_structure_phylogeny import build_upgma_tree
        names = ["A"]
        dist = np.array([[0.0]])
        newick = build_upgma_tree(names, dist)
        assert "A" in newick

    def test_rmsd_requires_same_shape(self):
        from protein_structure_phylogeny import compute_rmsd
        c1 = np.array([[0, 0, 0]], dtype=float)
        c2 = np.array([[0, 0, 0], [1, 1, 1]], dtype=float)
        with pytest.raises(AssertionError):
            compute_rmsd(c1, c2)
