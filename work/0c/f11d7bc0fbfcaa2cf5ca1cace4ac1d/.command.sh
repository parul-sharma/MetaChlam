#!/bin/bash -ue
python -c "
import pandas as pd
import glob

summary_files = glob.glob('*_summary.csv')
combined_df = pd.concat([pd.read_csv(f) for f in summary_files], ignore_index=True)
combined_df.to_csv('Combined_summary.csv', index=False)
        "
