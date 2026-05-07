# Small Experiment: Pipeline Validation on 5 Genomes

## Purpose
Quick end-to-end validation of the SAAP pipeline using a minimal dataset (5 bacterial genomes) to verify all stages work correctly before running the full 30-genome analysis.

## Genome Selection

| ID  | Species | Phylum | Rationale |
|-----|---------|--------|-----------|
| G01 | *Methanobrevibacter smithii* | Euryarchaeota | **Outgroup** (archaeal) |
| G02 | *Bacteroides thetaiotaomicron* | Bacteroidetes | Representative Bacteroidetes |
| G20 | *Faecalibacterium prausnitzii* | Firmicutes | Representative Firmicutes |
| G25 | *Escherichia coli* | Proteobacteria | Representative Proteobacteria (well-annotated) |
| G30 | *Thermotoga maritima* | Thermotogae | **Outgroup** (deep-branching) |

This selection covers 3 distinct phyla + 2 outgroups, providing a phylogenetically diverse test set that should produce a meaningful tree topology.

## Setup

### 1. Create test sample file

```bash
cd SnakeAnalysisPhylogenomicsPipeline-main
```

Create `config/samples_small.csv`:
```csv
id,species,strain,accession,bioproject,type,source
G01,Methanobrevibacter_smithii,OF4,GCF_000016525.1,PRJNA224116,Outgroup,Control
G02,Bacteroides_thetaiotaomicron,VPI_5482,GCF_000011065.1,PRJNA224116,Bacteroidetes,SLE
G20,Faecalibacterium_prausnitzii,M21/2,GCF_000154385.1,PRJNA18203,Firmicutes,Commensal
G25,Escherichia_coli,K-12_substr.MG1655,GCF_000005845.2,PRJNA225,Proteobacteria,RA
G30,Thermotoga_maritima,JGI,GCF_000230655.2,PRJNA63051,Outgroup,Control
```

### 2. Create test config

Create `config/config_small.yaml`:
```yaml
samples_file: "config/samples_small.csv"
workflow_mode: "supertree"
supertree_method: "fasturec"

experiments:
  include_paralogs: False
  filter_bootstraps: False
  min_tree_support: 0.3

clustering:
  min_seq_id: 0.25
  coverage: 0.6
  min_species: 3       # only 3 ingroup species
  outgroups: ["G01", "G30"]

alignment:
  mafft_args: "--auto --quiet"
  trimal_method: "-automated1"

consensus:
  minsup: 0.5
  threads: 2

gene_tree_params:
  ML:
    flags: "-m JTT+G4 -bb 500"
  NJ:
    flags: "-fastest -quiet"

tree:
  bootstrap: 50
  gene_tree_method: "ML"
  model: "LG+G"

tree_methods: ["ML"]

resources:
  download: 1
  clustering: 2
  orthologs: 2
  mafft: 1
  trimal: 1
  genetree: 1
  supertree: 2
```

### 3. Run the experiment

```bash
conda env create -f envs/sapp.yaml
conda activate sapp
snakemake --configfile config/config_small.yaml --cores 4 -p
```

## Validation Checklist

### Stage 1: Data Acquisition
- [ ] All 5 `.faa` files created in `results/proteomes/`
- [ ] Each file has size > 100 bytes
- [ ] Headers follow format `>GXX_ProteinID`

```bash
ls -la results/proteomes/G{01,02,20,25,30}.faa
grep -c "^>" results/proteomes/*.faa
```

### Stage 2: Clustering
- [ ] `results/clustering/all_proteins.fasta` created and non-empty
- [ ] `results/clustering/clusters_cluster.tsv` created
- [ ] Clusters contain proteins from multiple species

```bash
wc -l results/clustering/clusters_cluster.tsv
cut -f2 results/clustering/clusters_cluster.tsv | cut -d'_' -f1 | sort -u
```

### Stage 3: Ortholog Selection
- [ ] `results/sco/families/` directory populated
- [ ] `results/sco/orthologs_report.txt` shows reasonable stats
- [ ] At least 10+ gene families extracted

```bash
ls results/sco/families/ | wc -l
cat results/sco/orthologs_report.txt
```

### Stage 4: Alignment & Trimming
- [ ] `.aln` files in `results/sco/alignments/`
- [ ] `.trim.aln` files in `results/sco/trimmed/`
- [ ] Aligned files have consistent sequence lengths

```bash
ls results/sco/alignments/ | head
ls results/sco/trimmed/ | head
# Verify alignment length consistency
head -2 results/sco/trimmed/*.trim.aln | grep -c "^>"
```

### Stage 5: Gene Trees
- [ ] `.treefile` outputs in `results/trees/genetrees/sco/ML/`
- [ ] Tree files are non-empty and parseable

```bash
ls results/trees/genetrees/sco/ML/ | wc -l
# Quick parse test
python3 -c "
from Bio import Phylo
import glob
for f in glob.glob('results/trees/genetrees/sco/ML/*.treefile')[:3]:
    t = list(Phylo.parse(f, 'newick'))
    print(f'{f}: {len(t)} trees, {len(t[0].get_terminals())} leaves')
"
```

### Stage 6: Species Tree
- [ ] Consensus tree created (if applicable)
- [ ] Supertree/Fasturec tree created
- [ ] Final tree contains all 5 species

```bash
find results/trees/supertree/ -name "*.treefile" -o -name "*.newick" | head
python3 -c "
from Bio import Phylo
import glob
for f in glob.glob('results/trees/supertree/**/*.newick', recursive=True):
    t = list(Phylo.parse(f, 'newick'))
    if t:
        leaves = [l.name for l in t[0].get_terminals()]
        print(f'{f}: leaves={leaves}')
"
```

### Stage 7: Benchmarks
- [ ] `benchmarks/MASTER_BENCHMARK.csv` created
- [ ] Contains timing data for each pipeline stage

```bash
cat benchmarks/MASTER_BENCHMARK.csv
```

## Expected Tree Topology

Based on known bacterial phylogeny, the expected tree should be:
```
         ┌── Thermotoga maritima (outgroup)
    ┌────┤
    │    └── Methanobrevibacter smithii (outgroup)
────┤
    │    ┌── Bacteroides thetaiotaomicron
    └────┤
         │    ┌── Escherichia coli
         └────┤
              └── Faecalibacterium prausnitzii
```

Key phylogenetic relationships to verify:
1. Both outgroups (G01, G30) should separate from ingroup species
2. Bacteroidetes (G02) should branch before Firmicutes+Proteobacteria split
3. Firmicutes (G20) and Proteobacteria (G25) should be more closely related to each other than to Bacteroidetes

## Automated Validation Script

Save as `scripts/validate_small_experiment.py`:

```python
#!/usr/bin/env python3
"""Automated validation for the small experiment."""

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
```

## Troubleshooting

| Issue | Likely Cause | Fix |
|-------|-------------|-----|
| Download fails for a genome | NCBI datasets CLI not installed or rate-limited | Install via conda: `conda install -c conda-forge ncbi-datasets-cli` |
| 0 gene families extracted | `min_species` too high for test set | Ensure `min_species: 3` in config |
| Empty gene trees | Too few sequences in family after cleaning | Check `results/sco/trimmed/` for empty files |
| Fasturec fails | Not installed or not on PATH | See README installation instructions |
| IQ-TREE model error | Model not compatible with alignment | Try `-m MFP` for automatic model selection |

## Cleanup

```bash
# Remove all results from the small experiment
rm -rf results/ benchmarks/
```

---

## Stage 8: Protein Structure Phylogeny (ESMFold)

### Background

Traditional phylogenetic trees are built from **sequence alignments** — they detect evolutionary divergence in primary amino acid order. However, protein **3D structure** is more conserved than sequence: two proteins can diverge beyond detectable sequence similarity yet retain identical folds. Structure-based phylogenetics therefore reveals deeper evolutionary relationships invisible to sequence methods.

This pipeline uses the **ESMFold API** (Meta's protein language model) to predict 3D structures from sequences alone, then computes structural distances to build phylogenetic trees.

### Structural Distance Criteria

| Criterion | Formula | Range | Interpretation |
|-----------|---------|-------|----------------|
| **TM-score** | $\text{TM} = \frac{1}{L} \sum_{i} \frac{1}{1 + (d_i / d_0)^2}$ | 0–1 | >0.5 = same fold, >0.17 = same superfold |
| **RMSD** | $\text{RMSD} = \sqrt{\frac{1}{N} \sum_{i} \|\mathbf{r}_i - \mathbf{r}_i'\|^2}$ | 0–∞ Å | Lower = more similar; <2Å = near-identical |
| **GDT-TS** | $\frac{1}{4}(P_1 + P_2 + P_4 + P_8)$ | 0–1 | Fraction of Cα within 1/2/4/8 Å |
| **Combined** | $0.5 \cdot \text{TM} + 0.3 \cdot \text{GDT} + 0.2 \cdot \frac{1}{1 + \text{RMSD}}$ | 0–1 | Multi-metric consensus |

Where:
- $d_i$ = distance between aligned Cα atoms after Kabsch superposition
- $d_0 = 1.24 \sqrt[3]{L-15} - 1.8$ (length normalization)
- $P_t$ = fraction of residues within threshold $t$ Å

### Quality Assessment (pLDDT)

ESMFold outputs per-residue confidence (pLDDT) in the B-factor column:

| Tier | pLDDT Range | Meaning |
|------|-------------|---------|
| High | >90 | Very reliable structure |
| Good | 70–90 | Backbone reliable, side-chains approximate |
| Low | 50–70 | Topology uncertain |
| Very low | <50 | Likely disordered / unreliable |

**Only structures with mean pLDDT > 70 should be used for tree construction.**

### Running the Structure Pipeline (Small Experiment)

#### Step 1: Pick a conserved gene family

After the main pipeline runs, select a well-represented family:
```bash
# Find families with all 5 species (3 ingroup + 2 outgroup)
for f in results/sco/families/*.fasta; do
    n=$(grep -c "^>" "$f")
    if [ "$n" -eq 5 ]; then
        echo "$f ($n seqs)"
    fi
done | head -5
```

#### Step 2: Predict structures via ESMFold
```bash
FAMILY="G01_WP_011011408"  # example — replace with actual family name

python3 scripts/protein_structure_phylogeny.py predict \
    --fasta results/sco/families/${FAMILY}.fasta \
    --output_dir results/structures/sco/${FAMILY}/ \
    --max_length 400 \
    --assess_quality
```

Expected output:
```
[1/5] Processing G01_WP_011011408
  [G01_WP_011011408] Predicting structure (287 residues)...
  [G01_WP_011011408] Structure saved to results/structures/sco/.../G01_WP_011011408.pdb
...
Quality report saved to results/structures/sco/.../quality_report.json
  G01_...: pLDDT=82.3 (good)
  G02_...: pLDDT=91.1 (high)
```

#### Step 3: Build structure-based tree
```bash
python3 scripts/protein_structure_phylogeny.py tree \
    --manifest results/structures/sco/${FAMILY}/manifest.json \
    --output results/trees/structure/sco/${FAMILY}/tm_score.newick \
    --criterion tm_score \
    --method nj \
    --save_matrix
```

#### Step 4: Compare with sequence-based tree
```bash
python3 scripts/protein_structure_phylogeny.py compare \
    --sequence_tree results/trees/genetrees/sco/ML/${FAMILY}.treefile \
    --structure_tree results/trees/structure/sco/${FAMILY}/tm_score.newick
```

### Structure Validation Checklist

- [ ] PDB files created in `results/structures/sco/<family>/`
- [ ] `quality_report.json` shows mean pLDDT > 70 for most structures
- [ ] `manifest.json` maps all 5 proteins to PDB files
- [ ] Structure tree (`.newick`) is parseable and has 5 leaves
- [ ] Distance matrix (`.npz`) is saved alongside tree

```bash
# Verify structure outputs
ls results/structures/sco/${FAMILY}/*.pdb | wc -l
python3 -c "
import json
with open('results/structures/sco/${FAMILY}/quality_report.json') as f:
    q = json.load(f)
for name, info in q.items():
    if 'mean_plddt' in info:
        print(f'{name}: pLDDT={info[\"mean_plddt\"]:.1f} ({info[\"quality_tier\"]})')
"
```

### Expected Structural Relationships

For conserved housekeeping proteins (e.g., ribosomal, metabolic enzymes):
- **Ingroup** structures should cluster together (TM-score > 0.5)
- **Outgroup** structures may be more divergent structurally
- Structure tree topology should largely agree with sequence tree for well-folded proteins
- Disagreements between sequence and structure trees may indicate **convergent evolution** or **remote homology**

### Via Snakemake

Structure prediction and tree building are also available as Snakemake rules:
```bash
# Predict structure for a specific family
snakemake results/structures/sco/FAMILY_NAME/manifest.json \
    --configfile config/config_small.yaml --cores 1

# Build structure tree
snakemake results/trees/structure/sco/FAMILY_NAME/tm_score.newick \
    --configfile config/config_small.yaml --cores 1
```

---

## Stage 9: Hyperbolic Phylogenetic Embeddings

After the species tree is built, embed it into Poincaré ball space for conditioning the generative model.

```bash
python3 scripts/hyperbolic_embeddings.py \
    --tree results/trees/supertree/sco/raw/for_consensus/ML/consensus_tree.treefile \
    --output results/embeddings/sco/raw/ML/phylo_embeddings.npz \
    --dim 32 --epochs 300 --n_perturbations 3 --metrics
```

### Validation
- [ ] `phylo_embeddings.npz` created and non-empty
- [ ] Quality metrics: Spearman correlation > 0.8, Kruskal stress < 0.3
- [ ] All leaf nodes have embeddings with norm < 1.0 (inside Poincaré ball)

```bash
python3 -c "
import numpy as np
data = np.load('results/embeddings/sco/raw/ML/phylo_embeddings.npz', allow_pickle=True)
nodes = list(data['nodes'])
emb = data['base_embeddings']
norms = np.linalg.norm(emb, axis=1)
print(f'Nodes: {len(nodes)}, Leaves: {len(data[\"leaves\"])}')
print(f'Norms: min={norms.min():.4f}, mean={norms.mean():.4f}, max={norms.max():.4f}')
assert norms.max() < 1.0, 'Embeddings outside Poincaré ball!'
print('All embeddings inside Poincaré ball: OK')
n_pert = sum(1 for k in data.files if k.startswith('perturbed_'))
print(f'Perturbations: {n_pert}')
"
```

---

## Alternative Mode Testing

The main experiment uses SCO + supertree + ML. These variants test the remaining code paths.

### Variant A: Paralogs Mode
```bash
# Create config with paralogs enabled
cp config/config_small.yaml config/config_small_paralogs.yaml
```
Edit `config/config_small_paralogs.yaml`:
```yaml
experiments:
  include_paralogs: True
  filter_bootstraps: False
```
```bash
snakemake results/paralogs/orthologs_report.txt \
    --configfile config/config_small_paralogs.yaml --cores 4 -p
```
- [ ] `results/paralogs/families/` populated
- [ ] Families may contain `_p1`, `_p2` suffixed headers (multi-copy genes)
- [ ] `orthologs_report.txt` shows results in paralogs mode

### Variant B: NJ Gene Trees
```bash
cp config/config_small.yaml config/config_small_nj.yaml
```
Edit: `gene_tree_method: "NJ"` in config.
```bash
snakemake results/trees/genetrees/sco/NJ/ \
    --configfile config/config_small_nj.yaml --cores 4 -p
```
- [ ] `.treefile` outputs in `results/trees/genetrees/sco/NJ/`
- [ ] FastTree runs without errors

### Variant C: Bootstrap Filtering
```bash
cp config/config_small.yaml config/config_small_filter.yaml
```
Edit:
```yaml
experiments:
  filter_bootstraps: True
  min_tree_support: 0.3
```
```bash
snakemake results/trees/processed/sco/ML/filtered_list.txt \
    --configfile config/config_small_filter.yaml --cores 4 -p
```
- [ ] `filtered_list.txt` created
- [ ] Contains paths to trees passing the bootstrap threshold
- [ ] Filtered consensus/supertree trees are built

### Variant D: Supermatrix Mode
```bash
cp config/config_small.yaml config/config_small_supermatrix.yaml
```
Edit: `workflow_mode: "supermatrix"` and add `tree_methods: ["ML"]`.
```bash
snakemake results/trees/supermatrix/ML/final_tree.treefile \
    --configfile config/config_small_supermatrix.yaml --cores 4 -p
```
- [ ] `results/supermatrix/supermatrix.phy` created (PHYLIP format)
- [ ] `results/supermatrix/partitions.txt` created
- [ ] `results/trees/supermatrix/ML/final_tree.treefile` built
