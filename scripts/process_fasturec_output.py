#!/usr/bin/env python3
# process_fasturec_output.py
# -*- coding: utf-8 -*-


import sys
import argparse


def process_fasturec_output(input_file: str,
                            output_file: str
                            ) -> None:
    """ Processes a FASTUREC output file.
    """
    try:
        with open(input_file, 'r') as infile:
            lines = [line.strip() for line in infile if line.strip()]

        if not lines:
            print(f"Warning: Input file {input_file} is empty.")
            return

        best_line = lines[-1]

        if " " in best_line:
            parts = best_line.lstrip().split(' ', 1)
            if len(parts) == 2 and parts[0].replace('.', '', 1).isdigit():
                 processed_line = parts[1].strip()
            else:
                 processed_line = best_line
        else:
            processed_line = best_line

        if not processed_line.endswith(';'):
            processed_line += ';'

        with open(output_file, 'w') as outfile:
            outfile.write(processed_line + '\n')

        print(f"Success: Selected the last tree and saved to: {output_file}")

    except Exception as e:
        print(f"An error occurred: {e}")
        sys.exit(1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Processes Fasturec output, selecting only the last tree.')
    parser.add_argument('--fasturec_tree', required=True, help='Path to a file containing Fasturec output.')
    parser.add_argument('--output_tree', required=True, help='Desired output file path.')
    args = parser.parse_args()

    process_fasturec_output(args.fasturec_tree, args.output_tree)
