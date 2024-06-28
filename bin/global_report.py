#!/bin/python

'''USAGE: python global_report.py -st /path/to/strainge_file -lr /path/to/linreport_file -sm /path/to/sourmash_file -sc /path/to/strainscan_file -s sample_id -o /path/to/final_summary_report.csv'''

import pandas as pd
import os
import argparse

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
        linreport_filtered = linreport_df[linreport_df['Percentage_assigned_reads'] != 0]
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

def main(strainge_file, linreport_file, sourmash_file, strainscan_file, sample_id, output_file):

    # Initialize variables to store extracted data
    strainge_data = pd.DataFrame()
    linreport_data = pd.DataFrame()
    sourmash_data = pd.DataFrame()
    strainscan_data = pd.DataFrame()
    
    # Extract data from strainge file
    if strainge_file:
        strainge_data = extract_strainge(strainge_file)

    # Extract data from LINreport file
    if linreport_file:
        linreport_data = extract_linreport(linreport_file)

    # Extract data from sourmash file
    if sourmash_file:
        sourmash_data = extract_sourmash(sourmash_file)

    # Extract data from strainscan file
    if strainscan_file:
        strainscan_data = extract_strainscan(strainscan_file)
        
    # Create a dictionary to store data for this sample
    sample_summary = {'Sample': sample_id}
        
    # Add data from each tool to the sample summary dictionary
    if not strainge_data.empty:
        sample_summary.update(strainge_data.iloc[0].to_dict())
    if not linreport_data.empty:
        sample_summary.update(linreport_data.iloc[0].to_dict())
    if not sourmash_data.empty:
        sample_summary.update(sourmash_data.iloc[0].to_dict())
    if not strainscan_data.empty:
        sample_summary.update(strainscan_data.iloc[0].to_dict())

    # Convert dictionary to DataFrame
    summary_df = pd.DataFrame([sample_summary])  # Wrap sample_summary in a list

    # Save the final summary report to a CSV file
    summary_df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Summarize results from multiple tools into a single report.")
    parser.add_argument('-st', '--strainge', required=False, help="Path to the strainge results file.")
    parser.add_argument('-lr', '--linreport', required=False, help="Path to the LIN report file.")
    parser.add_argument('-sm', '--sourmash', required=False, help="Path to the sourmash results file.")
    parser.add_argument('-sc', '--strainscan', required=False, help="Path to the strainscan output directory")
    parser.add_argument('-s', '--sample', required=True, help="Sample_ID: unique identifier for the sample to summarize.")
    parser.add_argument('-o', '--output', required=True, help="Path and name of the final summary report file.")
    
    args = parser.parse_args()
    
    main(args.strainge, args.linreport, args.sourmash, args.strainscan, args.sample, args.output)
