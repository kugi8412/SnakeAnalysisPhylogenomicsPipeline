#!/usr/bin/env python3
# download_genomes.py
# -*- coding: utf-8 -*- 


import os
import sys
import shutil
import subprocess

import pandas as pd


if snakemake.log:
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

def download_genomes(samples_csv: str,
                     output_dir: str
                    ) -> None:
    df = pd.read_csv(samples_csv).fillna("")
    os.makedirs(output_dir, exist_ok=True)
    for _, row in df.iterrows():
        sample_id = row['id']
        accession = str(row['accession']).strip()
        output_file = os.path.join(output_dir, f"{sample_id}.faa")

        if os.path.exists(output_file) and os.path.getsize(output_file) > 100:
            print(f"[{sample_id}] File exists. Skipping.", file=sys.stderr)
            continue

        print(f"[{sample_id}] Downloading {accession}", file=sys.stderr)
        temp_zip = f"{sample_id}.zip"
        temp_dir = f"temp_{sample_id}"
        try:
            cmd = [
                "datasets", "download", "genome", "accession", accession,
                "--include", "protein",
                "--filename", temp_zip
            ]
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            subprocess.run(["unzip", "-o", "-q", temp_zip, "-d", temp_dir], check=True)
            found = False
            for root, _, files in os.walk(temp_dir):
                for file in files:
                    if file.endswith(".faa") or file.endswith(".GP_AA"):
                        source_path = os.path.join(root, file)
                        sanitize_header(source_path, sample_id, output_file)
                        found = True
                        break
            
            if not found:
                print(f"[{sample_id}] WARNING: No protein file found inside zip!", file=sys.stderr)
                open(output_file, 'w').close()

        except subprocess.CalledProcessError as e:
            print(f"[{sample_id}] DOWNLOAD ERROR. Cmd: {' '.join(cmd)}", file=sys.stderr)
            open(output_file, 'w').close()   
        except Exception as e:
            print(f"[{sample_id}] ERROR: {e}", file=sys.stderr)
            open(output_file, 'w').close()
        finally:
            if os.path.exists(temp_zip): os.remove(temp_zip)
            if os.path.exists(temp_dir): shutil.rmtree(temp_dir)

if __name__ == "__main__":
    download_genomes(snakemake.params.samples, snakemake.params.out_dir)
