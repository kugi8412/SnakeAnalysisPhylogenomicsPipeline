#!/usr/bin/env python3
# download_genomes.py
# -*- coding: utf-8 -*- 


import os
import sys
import shutil
import subprocess
import pandas as pd


if "snakemake" in locals() and snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

def sanitize_header(fasta_path: str,
                    sample_id: str,
                    output_path: str
                    ) -> bool:
    try:
        with open(fasta_path, "r") as f_in, open(output_path, "w") as f_out:
            for line in f_in:
                if line.startswith(">"):
                    clean_id = line.strip()[1:].split()[0]
                    f_out.write(f">{sample_id}_{clean_id}\n")
                else:
                    f_out.write(line)
        return True
    except Exception as e:
        print(f"[{sample_id}] Error writing file: {e}", file=sys.stderr)
        return False


def run_datasets_cli(identifier: str, 
                     temp_zip: str, 
                     temp_dir: str, 
                     sample_id: str, 
                     output_file: str) -> bool:
    cmd = [
        "datasets", "download", "genome", "accession", identifier,
        "--include", "protein",
        "--filename", temp_zip
    ]
    
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        subprocess.run(["unzip", "-o", "-q", temp_zip, "-d", temp_dir], check=True)
        found = False
        for root, _, files in os.walk(temp_dir):
            for file in files:
                if file.endswith(".faa") or file.endswith(".GP_AA"):
                    source_path = os.path.join(root, file)
                    if sanitize_header(source_path, sample_id, output_file):
                        found = True
                    break
            if found: break
            
        return found

    except subprocess.CalledProcessError:
        return False
    except Exception as e:
        print(f"[{sample_id}] Unexpected error processing {identifier}: {e}", file=sys.stderr)
        return False
    finally:
        if os.path.exists(temp_zip): os.remove(temp_zip)
        if os.path.exists(temp_dir): shutil.rmtree(temp_dir)


def download_genomes(samples_csv: str, output_dir: str) -> None:
    df = pd.read_csv(samples_csv).fillna("")
    os.makedirs(output_dir, exist_ok=True)

    for _, row in df.iterrows():
        sample_id = str(row['id'])
        accession = str(row.get('accession', '')).strip()
        bioproject = str(row.get('bioproject', '')).strip()
        species = str(row.get('species', 'Unknown Species')).strip()
        output_file = os.path.join(output_dir, f"{sample_id}.faa")

        if os.path.exists(output_file) and os.path.getsize(output_file) > 100:
            print(f"[{sample_id}] File exists. Skipping.", file=sys.stderr)
            continue

        temp_zip = f"{sample_id}.zip"
        temp_dir = f"temp_{sample_id}"
        success = False

        # Accession
        if accession:
            print(f"[{sample_id}] Attempt 1: Downloading by Accession {accession}", file=sys.stderr)
            success = run_datasets_cli(accession, temp_zip, temp_dir, sample_id, output_file)
        
        # BioProject
        if not success:
            if bioproject:
                print(f"[{sample_id}] [WARNING]: Accession failed or missing.", file=sys.stderr)
                print(f"[{sample_id}] Fallback: Trying BioProject '{bioproject}' for species '{species}'", file=sys.stderr)
                success = run_datasets_cli(bioproject, temp_zip, temp_dir, sample_id, output_file)
            else:
                print(f"[{sample_id}] Skipped Fallback: No BioProject provided.", file=sys.stderr)

        if not success:
            print(f"[{sample_id}] [ERROR]: Could not retrieve genome (Accession: {accession}, BioProject: {bioproject})", file=sys.stderr)
            open(output_file, 'w').close()
        else:
            print(f"[{sample_id}] SUCCESS.", file=sys.stderr)


if __name__ == "__main__":
    download_genomes(snakemake.params.samples, snakemake.params.out_dir)
