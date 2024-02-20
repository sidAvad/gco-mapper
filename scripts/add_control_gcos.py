import argparse
import pandas as pd
import numpy as np

def add_control_gcos(window_df,control_gcos_df,output_file):
    assert len(window_df) > 1 #This only makes sense to run once windows have been created
    window_df['gco_control']=[[] for _ in range(len(window_df))]
    i,j=0,0
    while i < len(window_df) and j < len(control_gcos_df):
        print(j)
        print(len(control_gcos_df))
        gco_row=control_gcos_df.iloc[j]
        gco_midpoint = (gco_row['first'] +  gco_row['last'])/2
        if gco_midpoint < window_df.iloc[i]['window_end_marker']:
            j+=1
        elif gco_midpoint > window_df.iloc[i]['window_end_marker']:
            i+=1
        else:
            window_df.iloc[i]['gco_control'].append({'pop':gco_row['pop'],\
                    'first_marker':gco_row['first'],\
                    'last_marker':gco_row['last'],\
                    'diffs':gco_row['diffs'],\
                    'bg':gco_row['bg']\
                    })
            j+=1 #advance gco because there may be more than one in the window.
            
    window_df.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Process scanning file.")
    parser.add_argument('-w', '--window_file', required=True, help="Input file path.")
    parser.add_argument('-c', '--control_gcos_file', required=True, help="Input file path.")
    parser.add_argument('-o', '--output_file', required=True, help="Output file path.")
    
    args = parser.parse_args()

    window_df = pd.read_table(args.window_file)
    control_gcos_df = pd.read_csv(args.control_gcos_file, sep=r'\s+')
    add_control_gcos(window_df, control_gcos_df,args.output_file)
    #Save the processed data to the output file

if __name__ == "__main__":
    main()
