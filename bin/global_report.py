#!/bin/python

'''USAGE: python global_report.py -kl /path/to/kraken_log -st /path/to/strainge_file -lr /path/to/linreport_file -sm /path/to/sourmash_file -sc /path/to/strainscan_file -s sample_id -o /path/to/final_summary_rep
ort.csv'''

import pandas as pd
import os
import argparse
from pathlib import Path
import re
import numpy as np

# Function to extract relevant data from strainge results
def extract_strainge(file):
    try:
        strainge_df = pd.read_csv(file, sep='\t')
        strainge_data = strainge_df.loc[:, ['strain', 'gcov']]
        return strainge_data
    except Exception as e:
        print(f"Error reading strainge file {file}: {e}")
        return pd.DataFrame()

# Function to extract relevant data from LIN report
def extract_linreport(file):
    try:
        linreport_df = pd.read_csv(file, sep=',')
        linreport_filtered = linreport_df[(linreport_df['Percentage_assigned_reads'] > 0) & linreport_df['Percentage_assigned_reads'].notna()]
               
        linreport_data = linreport_filtered.loc[:, ['LINgroup_Name', 'Assigned_reads']]
        return linreport_data
    except Exception as e:
        print(f"Error reading LIN report file {file}: {e}")
        return pd.DataFrame()

# Function to extract relevant data from sourmash output
def extract_sourmash(file):
    try:
        sourmash_df = pd.read_csv(file, sep=',')
        sourmash_data = sourmash_df.loc[:, ['name', 'f_orig_query']]
        return sourmash_data
    except Exception as e:
        print(f"Error reading sourmash file {file}: {e}")
        return pd.DataFrame()

# Function to extract relevant data from strainscan output
def extract_strainscan(file):
    try:
        strainscan_df = pd.read_csv(file, sep='\t')
        strainscan_data = strainscan_df.loc[:, ['Strain_Name', 'Predicted_Depth (Enet)']]
        return strainscan_data
    except Exception as e:
        print(f"Error reading strainscan file {file}: {e}")
        return pd.DataFrame()

# Function to extract classified reads from Kraken log
def extract_kraken_classified(log_file):
    with open(log_file, 'r') as file:
        log_content = file.read()
    
    # Use regular expressions to find the number and percentage of classified reads
    classified_pattern = re.compile(r'^\s*(\d+)\s+sequences classified \(([\d\.]+)%\)', re.MULTILINE)

    classified_match = classified_pattern.search(log_content)

    if classified_match:
        classified_reads = classified_match.group(1)
        classified_percentage = classified_match.group(2)

        return {
            "classified_reads": classified_reads,
            "classified_percentage": classified_percentage
        }
    else:
        raise ValueError("Unable to find classified and unclassified read counts in the log file.")

def main(kraken_file, strainge_file, linreport_file, sourmash_file, strainscan_dir, sample_id, output_file):
    # Initialize variables to store extracted data
    strainge_data = extract_strainge(strainge_file)
    linreport_data = extract_linreport(linreport_file)
    sourmash_data = extract_sourmash(sourmash_file)
    strainscan_data = extract_strainscan(Path(strainscan_dir) / 'final_report.txt')

    # Extract classified reads from kraken report file
    kraken_data = pd.DataFrame()
    if kraken_file:
        try:
            kraken_results = extract_kraken_classified(kraken_file)
            kraken_data = pd.DataFrame({
                'Classified_Reads': [kraken_results['classified_reads']],
                'Classified_Percentage': [kraken_results['classified_percentage']]
            }, dtype=str)
        except ValueError as e:
            print(e)

    # Create a dataframe to store Sample ID 
    sample_summary = pd.DataFrame({'Sample': [sample_id]})

    # Combine data from each tool
    combined_data = pd.concat([sample_summary, kraken_data, strainge_data, linreport_data, sourmash_data, strainscan_data], axis=1)

    # Save the final summary report to a CSV file
    combined_data.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize results from multiple tools into a single report.")
    parser.add_argument('-kl', '--kraken_log', required=False, help="Path to the kraken log file.")
    parser.add_argument('-st', '--strainge', required=False, help="Path to the strainge results file.")
    parser.add_argument('-lr', '--linreport', required=False, help="Path to the LIN report file.")
    parser.add_argument('-sm', '--sourmash', required=False, help="Path to the sourmash results file.")
    parser.add_argument('-sc', '--strainscan', required=False, help="Path to the strainscan output directory")
    parser.add_argument('-s', '--sample', required=True, help="Sample_ID: unique identifier for the sample to summarize.")
    parser.add_argument('-o', '--output', required=True, help="Path and name of the final summary report file.")
    
    args = parser.parse_args()
    
    main(args.kraken_log, args.strainge, args.linreport, args.sourmash, args.strainscan, args.sample, args.output)

