#!/usr/bin/env python3
# gather_benchmarks.py
# -*- coding: utf-8 -*-


import os
import pandas as pd
import glob


benchmark_dir = "benchmarks"
output_file = snakemake.output[0]

all_data = []

for file in glob.glob(os.path.join(benchmark_dir, "*.tsv")):
    if not os.path.isfile(file):
        continue
    
    process_name = os.path.basename(file).replace(".tsv", "")
    
    try:
        df = pd.read_csv(file, sep="\t")
        if not df.empty:
            row = df.iloc[0].to_dict()
            row["Process"] = process_name
            all_data.append(row)
    except Exception as e:
        print(f"Skipping empty or bad file: {file}")

if all_data:
    final_df = pd.DataFrame(all_data)
    cols = ["Process", "s", "h:m:s", "max_rss", "max_vms", "io_in", "io_out"]
    cols = [c for c in cols if c in final_df.columns]
    final_df = final_df[cols]
    final_df.rename(columns={"s": "Time_Seconds", "max_rss": "Peak_Memory_MB"}, inplace=True)
    final_df.to_csv(output_file, index=False)
    print(f"Benchmark summary saved to {output_file}")
else:
    print("[WARNING]: No benchmark data found.")
    with open(output_file, 'w') as f:
        f.write("Process,Time_Seconds,Peak_Memory_MB\n")
