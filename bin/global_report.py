#!/bin/python

'''USAGE:python global_report.py -i /path/to/samples/ -o /path/to/final_summary_report.csv'''


import pandas as pd
import glob
import os
import argparse

# Function to extract relevant data from strainge results
def extract_strainge(file):
    strainge_df = pd.read_csv(file, sep='\t')
    strainge_data = strainge_df.loc[:, ['strain', 'gcov']]
    return strainge_data

# Function to extract relevant data from LIN report
def extract_linreport(file):
    linreport_df = pd.read_csv(file, sep=',')
    linreport_filtered = linreport_df[linreport_df['Percentage_assigned_reads'] != 0]
    linreport_data = linreport_filtered.loc[:, ['LINgroup_Name', 'Assigned_reads']]
    return linreport_data

# Function to extract relevant data from sourmash output
def extract_sourmash(file):
    sourmash_df = pd.read_csv(file, sep=',')
    sourmash_data = sourmash_df.loc[:, ['name', 'f_orig_query']]
    return sourmash_data

# Function to extract relevant data from strainscan output
def extract_strainscan(file):
    strainscan_df = pd.read_csv(file, sep='\t')
    strainscan_data = strainscan_df.loc[:, ['Strain_Name', 'Predicted_Depth (Enet)']]
    return strainscan_data

def main(base_dir, sample_id, output_file):

    # Initialize variables to store extracted data
    strainge_data = None
    linreport_data = None
    sourmash_data = None
    strainscan_data = None
    
    # Check if 'strainge_results' file exists for the specified sample_id
    strainge_file = os.path.join(base_dir, f"{sample_id}.strains.tsv")
    if os.path.exists(strainge_file):
        strainge_data = extract_strainge(strainge_file)

    # Check if 'LINreport' file exists for the specified sample_id
    linreport_file = os.path.join(base_dir, f"{sample_id}.LINreport.txt")
    if os.path.exists(linreport_file):
        linreport_data = extract_linreport(linreport_file)

    # Check if 'sourmash_output' file exists for the specified sample_id
    sourmash_file = os.path.join(base_dir, f"{sample_id}.gather")
    if os.path.exists(sourmash_file):
        sourmash_data = extract_sourmash(sourmash_file)

    # Check if 'strainscan_output' file exists for the specified sample_id
    strainscan_file = os.path.join(base_dir, 'strainscan_output', 'final_report.txt')
    if os.path.exists(strainscan_file):
        strainscan_data = extract_strainscan(strainscan_file)
        
    # Create a dictionary to store data for this sample
    sample_summary = {'Sample': sample_id}
        
    # Add data from each tool to the sample summary dictionary
    if strainge_data is not None and not strainge_data.empty:
        sample_summary.update(strainge_data.iloc[0].to_dict())
    if linreport_data is not None and not linreport_data.empty:
        sample_summary.update(linreport_data.iloc[0].to_dict())
    if sourmash_data is not None and not sourmash_data.empty:
        sample_summary.update(sourmash_data.iloc[0].to_dict())
    if strainscan_data is not None and not strainscan_data.empty:
        sample_summary.update(strainscan_data.iloc[0].to_dict())

    # Convert dictionary to DataFrame
    summary_df = pd.DataFrame([sample_summary])  # Wrap sample_summary in a list

    # Save the final summary report to a CSV file
    summary_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize results from multiple tools into a single report.")
    parser.add_argument('-i', '--input', required=True, help="Path to the directory containing sample results.")
    parser.add_argument('-s', '--sample', required=True, help="Sample_ID: unique identifier for the sample to summarize.")
    parser.add_argument('-o', '--output', required=True, help="Path and name of the final summary report file.")
    
    args = parser.parse_args()
    
    main(args.input, args.sample, args.output)