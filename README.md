<div align="center">
  <img src="https://github.com/kugi8412/SnakeAnalysisPhylogenomicsPipeline/blob/main/SAPP_logo.png?raw=true" alt="SAAP Logo" width="700"/>
</div>

## Summary

![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)
![Biopython](https://img.shields.io/badge/Biopython-50C878?style=for-the-badge&logo=python&logoColor=white)
![Bash](https://img.shields.io/badge/Bash-3c4549?style=for-the-badge&logo=gnubash&logoColor=white)
![Conda](https://img.shields.io/badge/Conda-44A833?style=for-the-badge&logo=anaconda&logoColor=white)
![IQ-TREE](https://img.shields.io/badge/IQ--TREE-ML-blue?style=for-the-badge)


**SAAP** is a modular phylogenomic pipeline designed to automate the reconstruction of species trees from proteomic data. Unlike standard pipelines that rely on pre-defined marker genes, SAAP constructs trees from whole-proteome clustering, ensuring a comprehensive evolutionary signal.

It integrates state-of-the-art tools to process raw protein sequences into high-quality phylogenetic trees using **Consensus** or **Supertree** approaches. The pipeline is optimized for deep phylogenetic relationships (e.g., Bacteria/Archaea) by utilizing Maximum Likelihood methods with protein-specific substitution models via **IQ-TREE**.

<p align="center">
  <em>Example output: Reference taxonomy vs. Supertree inference.</em>
</p>

<p align="center">
  <img src="https://github.com/kugi8412/SnakeAnalysisPhylogenomicsPipeline/blob/main/Example.png?raw=true"
       alt="Example"
       width="700">
</p>


## Workflow Overview

1.  **Data Acquisition:** Fetching proteomes from NCBI Datasets.
2.  **Clustering:** Ortholog detection using **MMseqs2**.
3.  **Alignment:** Accurate MSA using **MAFFT**.
4.  **Trimming:** Removing spurious regions with **trimAI**.
5.  **Gene Trees:** Constructing individual gene trees using **IQ-TREE** (Maximum Likelihood) with automatic model selection (ModelFinder).
6.  **Species Tree:**
    * **Consensus:** Using **IQ-TREE** (Majority Rule / Extended Consensus).
    * **Supertree:** Using **Fasturec** or **Supermatrix** approaches.

## Installation

### 1. Clone repository

```bash
git clone https://github.com/kugi8412/SnakeAnalysisPhylogenomicsPipeline.git
cd SAPP
```

### 2. Install FASTUREC

```bash
sudo apt update
sudo apt install build-essential wget
git clone https://bitbucket.org/pgor17/fasturec.git
cd fasturec
make
export PATH="/opt/fasturec/bin:${PATH}"
```

### 3. Create conda environment
```bash
conda env create -f envs/environment.yml
conda activate sapp
```

### 4. Define config and samples
You can change snakefile to choose your own files.
Config contains parameters and methods used in phylogenetics, while samples contain species to be processed.

### 5. Run & Enjoy
```bash
snakemake --cores 8 --use-conda
```
