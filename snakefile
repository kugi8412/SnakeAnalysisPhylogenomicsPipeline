# SAAP
configfile: "config/config.yaml"


import pandas as pd


samples_df = pd.read_csv(config["samples_file"])
SAMPLES = samples_df['id'].tolist()
MODE = config["workflow_mode"]
METHODS = config["tree_methods"]
FILTERS = ["raw"]
DATASETS = ["sco"]

if config["experiments"]["include_paralogs"]:
    DATASETS.append("paralogs")

if config["experiments"]["filter_bootstraps"]:
    FILTERS.append("filtered")

tree_targets = []
if MODE == "supermatrix":
    tree_targets = expand("results/trees/supermatrix/{method}/final_tree.treefile", method=METHODS)
elif MODE == "supertree":
    combinations = expand(
        "{dataset}/{filter}/{method}",
        dataset=DATASETS, filter=FILTERS, method=config["tree"]["gene_tree_method"]
    )
    tree_targets.extend([f"results/trees/supertree/{c}/fasturec_tree.newick" for c in combinations])
    tree_targets.extend([f"results/trees/supertree/{c}/consensus_tree.treefile" for c in combinations])

rule all:
    input: "benchmarks/MASTER_BENCHMARK.csv"

rule download_genomes:
    output: flag=touch("results/proteomes/.download_complete"), proteomes=expand("results/proteomes/{sample}.faa", sample=SAMPLES)
    params: samples=config["samples_file"], out_dir="results/proteomes"
    threads: config["resources"]["download"]
    conda: "envs/sapp.yaml"
    benchmark: "benchmarks/download_genomes.tsv"
    script: "scripts/download_genomes.py"

rule concat_proteomes:
    input: "results/proteomes/.download_complete"
    output: "results/clustering/all_proteins.fasta"
    params: src_dir="results/proteomes", script="scripts/concat_files.sh"
    shell: "bash {params.script} '{params.src_dir}/*.faa' {output}"

rule cluster_proteins:
    input: "results/clustering/all_proteins.fasta"
    output: tsv="results/clustering/clusters_cluster.tsv"
    params: out_prefix="results/clustering/clusters", tmp_dir="results/clustering/tmp", min_id=config["clustering"]["min_seq_id"], cov=config["clustering"]["coverage"], script="scripts/run_clustering.sh"
    threads: config["resources"]["clustering"]
    conda: "envs/sapp.yaml"
    benchmark: "benchmarks/clustering.tsv"
    shell: "bash {params.script} {input} {params.out_prefix} {params.tmp_dir} {params.min_id} {params.cov} {threads}"

checkpoint extract_families:
    input: clusters="results/clustering/clusters_cluster.tsv", fasta="results/clustering/all_proteins.fasta"
    output: out_dir=directory("results/{dataset}/families"), stats="results/{dataset}/orthologs_report.txt"
    params: min_species=config["clustering"]["min_species"], outgroups=config["clustering"]["outgroups"], mode="{dataset}"
    threads: config["resources"]["orthologs"]
    conda: "envs/sapp.yaml"
    benchmark: "benchmarks/extract_{dataset}.tsv"
    script: "scripts/select_orthologs.py"

rule run_mafft:
    input: "results/{dataset}/families/{family}.fasta"
    output: "results/{dataset}/alignments/{family}.aln"
    params: args=config["alignment"]["mafft_args"], script="scripts/run_mafft.sh"
    threads: config["resources"]["mafft"]
    conda: "envs/sapp.yaml"
    benchmark: "benchmarks/mafft_{dataset}_{family}.tsv"
    shell: "bash {params.script} {input} {output} {threads} '{params.args}'"

rule run_trimal:
    input: "results/{dataset}/alignments/{family}.aln"
    output: "results/{dataset}/trimmed/{family}.trim.aln"
    params: method=config["alignment"]["trimal_method"], script="scripts/run_trimal.sh"
    threads: config["resources"]["trimal"]
    conda: "envs/sapp.yaml"
    benchmark: "benchmarks/trimal_{dataset}_{family}.tsv"
    shell: "bash {params.script} {input} {output} '{params.method}'"

# Gene trees
rule build_gene_trees:
    input: aln="results/{dataset}/trimmed/{family}.trim.aln"
    output: tree="results/trees/genetrees/{dataset}/{method}/{family}.treefile"
    params: method="{method}", flags=lambda w: config["gene_tree_params"][w.method]["flags"]
    threads: config["resources"]["genetree"]
    conda: "envs/sapp.yaml"
    benchmark: "benchmarks/genetree_{dataset}_{method}_{family}.tsv"
    script: "scripts/run_genetree.py"

# Helpery functions
def get_families(wildcards):
    checkpoint_output = checkpoints.extract_families.get(dataset=wildcards.dataset).output.out_dir
    import glob, os
    return [os.path.basename(f).replace(".fasta", "") for f in glob.glob(os.path.join(checkpoint_output, "*.fasta"))]

def get_raw_trees(wildcards):
    families = get_families(wildcards)
    return expand("results/trees/genetrees/{dataset}/{method}/{family}.treefile", 
                  dataset=wildcards.dataset, method=wildcards.method, family=families)


# Filtering trees
rule filter_gene_trees:
    input: 
        trees = get_raw_trees
    output: 
        list_file = "results/trees/processed/{dataset}/{method}/filtered_list.txt"
    params: 
        threshold = config["experiments"]["min_tree_support"],
        samples_csv = config["samples_file"],
        checker_script = "scripts/filter_trees.py"
    conda: "envs/sapp.yaml"
    shell:
        """
        > {output.list_file}
        for tree in {input.trees}; do
            python3 {params.checker_script} "$tree" {params.threshold} {params.samples_csv} >> {output.list_file}
        done
        """

tree_targets = []
METHOD = config.get("tree", {}).get("gene_tree_method", "NJ")

if MODE == "supertree":
    # SCO -> Consensus (Raw + Filtered)
    tree_targets.append(f"results/trees/supertree/sco/raw/for_consensus/{METHOD}/consensus_tree.treefile")
    tree_targets.append(f"results/trees/supertree/sco/filtered/for_consensus/{METHOD}/consensus_tree.treefile")
    # SCO -> Supertree (Raw + Filtered)
    tree_targets.append(f"results/trees/supertree/sco/raw/for_supertree/{METHOD}/fasturec_tree.newick")
    tree_targets.append(f"results/trees/supertree/sco/filtered/for_supertree/{METHOD}/fasturec_tree.newick")
    # Paraloges -> Supertree Fasturec (Raw + Filtered)
    tree_targets.append(f"results/trees/supertree/paralogs/raw/for_supertree/{METHOD}/fasturec_tree.newick")
    tree_targets.append(f"results/trees/supertree/paralogs/filtered/for_supertree/{METHOD}/fasturec_tree.newick")


rule prepare_supertree_input:
    input: 
        trees = get_raw_trees
    output: 
        merged = "results/trees/supertree/{dataset}/{filter}/{purpose}/{method}/all_genetrees.newick"
    params:
        purpose = "{purpose}",
        filter_type = "{filter}",
        min_species = config["clustering"]["min_species"],
        outgroups = ",".join(config["clustering"]["outgroups"]),
        samples_csv = config["samples_file"],
        bs_thresh = config["experiments"]["min_tree_support"],
        script = "scripts/merge_trees.py"
    conda: "envs/sapp.yaml"
    shell:
        """
        export MIN_SPECIES="{params.min_species}"
        export SAMPLES_CSV="{params.samples_csv}"
        export OUTGROUPS="{params.outgroups}"
        python3 {params.script} {output.merged} {params.purpose} {params.filter_type} {params.bs_thresh} {input.trees}
        """

# FASTUREC
rule build_fasturec:
    input:
        trees = "results/trees/supertree/{dataset}/{filter}/for_supertree/{method}/all_genetrees.newick"
    output:
        tree = "results/trees/supertree/{dataset}/{filter}/for_supertree/{method}/fasturec_tree.newick"
    params:
        script = "scripts/run_fasturec_pipeline.sh",
        post_process = "scripts/process_fasturec_output.py"
    shell:
        """
        bash {params.script} "{input.trees}" "{output.tree}" "{params.post_process}"
        """

# Consensus IQ-TREE
rule build_consensus:
    input:
        trees = "results/trees/supertree/{dataset}/{filter}/for_consensus/{method}/all_genetrees.newick"
    output:
        tree = "results/trees/supertree/{dataset}/{filter}/for_consensus/{method}/consensus_tree.treefile"
    params:
        minsup = config["consensus"]["minsup"],
        threads = config["resources"]["supertree"],
        script = "scripts/run_consensus.sh"
    conda: "envs/sapp.yaml"
    shell:
        """
        bash {params.script} "{input.trees}" "{output.tree}" "{params.minsup}" "{params.threads}"
        """

# SUPERMATRIX
def get_trimmed_msas(wildcards):
    families = get_families(wildcards) 
    return expand("results/orthologs/trimmed/{family}.trim.aln", family=families)


rule concat_matrix:
    input: get_trimmed_msas
    output: matrix="results/supermatrix/supermatrix.phy", partitions="results/supermatrix/partitions.txt"
    conda: "envs/sapp.yaml"
    script: "scripts/concat_matrix.py"

rule build_supermatrix_tree:
    input: matrix="results/supermatrix/supermatrix.phy"
    output: "results/trees/supermatrix/{method}/final_tree.treefile"
    params: 
        method="{method}", 
        model=config["tree"]["model"], 
        bootstrap=config["tree"]["bootstrap"], 
        script="scripts/run_supertree.sh"
    threads: config["resources"]["supertree"]
    conda: "envs/sapp.yaml"
    shell: "bash {params.script} {input.matrix} {output} {params.method} {threads} '{params.model}' {params.bootstrap}"

rule gather_benchmarks:
    input: tree_targets
    output: "benchmarks/MASTER_BENCHMARK.csv"
    script: "scripts/gather_benchmarks.py"
