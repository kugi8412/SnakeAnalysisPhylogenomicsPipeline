#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Automated validation for the small experiment.
Checks for the presence and basic integrity of key output files at each stage of the pipeline.
"""

import os
import sys
import glob


def check(condition, msg):
    status = "PASS" if condition else "FAIL"
    print(f"  [{status}] {msg}")
    return condition


def main():
    passed = 0
    total = 0

    print("=== SAAP Small Experiment Validation ===\n")

    # Stage 1: Proteomes
    print("Stage 1: Data Acquisition")
    for gid in ["G01", "G02", "G20", "G25", "G30"]:
        f = f"results/proteomes/{gid}.faa"
        total += 1
        if check(os.path.exists(f) and os.path.getsize(f) > 100, f"{gid}.faa exists and non-empty"):
            passed += 1

    # Stage 2: Clustering
    print("\nStage 2: Clustering")
    total += 1
    f = "results/clustering/clusters_cluster.tsv"
    if check(os.path.exists(f) and os.path.getsize(f) > 0, "Cluster TSV created"):
        passed += 1

    # Stage 3: Orthologs
    print("\nStage 3: Ortholog Selection")
    families = glob.glob("results/sco/families/*.fasta")
    total += 1
    if check(len(families) >= 5, f"At least 5 gene families (found {len(families)})"):
        passed += 1

    total += 1
    report = "results/sco/orthologs_report.txt"
    if check(os.path.exists(report), "Orthologs report generated"):
        passed += 1

    # Stage 4: Alignment
    print("\nStage 4: Alignment & Trimming")
    alns = glob.glob("results/sco/alignments/*.aln")
    trimmed = glob.glob("results/sco/trimmed/*.trim.aln")
    total += 2
    if check(len(alns) > 0, f"Alignments created ({len(alns)} files)"):
        passed += 1
    if check(len(trimmed) > 0, f"Trimmed alignments created ({len(trimmed)} files)"):
        passed += 1

    # Stage 5: Gene Trees
    print("\nStage 5: Gene Trees")
    trees = glob.glob("results/trees/genetrees/sco/ML/*.treefile")
    total += 1
    if check(len(trees) > 0, f"Gene trees built ({len(trees)} files)"):
        passed += 1

    non_empty = [t for t in trees if os.path.getsize(t) > 0]
    total += 1
    if check(len(non_empty) > 0, f"Non-empty gene trees ({len(non_empty)})"):
        passed += 1

    # Stage 6: Species Tree
    print("\nStage 6: Species Tree")
    species_trees = glob.glob("results/trees/supertree/**/*.newick", recursive=True)
    species_trees += glob.glob("results/trees/supertree/**/*.treefile", recursive=True)
    total += 1
    if check(len(species_trees) > 0, f"Species trees created ({len(species_trees)})"):
        passed += 1

    # Stage 7: Benchmarks
    print("\nStage 7: Benchmarks")
    total += 1
    bm = "benchmarks/MASTER_BENCHMARK.csv"
    if check(os.path.exists(bm) and os.path.getsize(bm) > 0, "Master benchmark CSV created"):
        passed += 1

    # Summary
    print(f"\n{'='*40}")
    print(f"Results: {passed}/{total} checks passed")
    if passed == total:
        print("All checks PASSED!")
    else:
        print(f"WARNING: {total - passed} check(s) FAILED")
        sys.exit(1)


if __name__ == "__main__":
    main()
