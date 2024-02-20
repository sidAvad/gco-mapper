import pandas as pd
import random
import sys

def randomize_gcotables(input_file):
    # Read the input text file into a DataFrame
    df = pd.read_csv(input_file,sep=r'\s+',index_col=0,engine='python')
    
    # Compute the range of the 'last' column
    last_range = df['last'].iloc[-1] - df['first'].iloc[0]
    
    # Compute the average length
    avg_length = last_range / len(df)
    
    # Create a new DataFrame with selected values
    new_data = []
    for index, row in df.iterrows():
        random_first = random.randint(df['first'].iloc[0], df['last'].iloc[-1])
        new_row = {
            'hap': row['hap'],
            'chrom': row['chrom'],
            'first': random_first,
            'last': random_first + int(avg_length)
        }
        new_data.append(new_row)
    
    new_df = pd.DataFrame(new_data)
    new_df = new_df[['chrom','hap','first','last']]
    
    # Write the new DataFrame to the specified output file
    new_df.to_csv(output_file, sep='\t', index=False,)

if __name__ == '__main__':
    if len(sys.argv) != 3:
    	print("Usage: python randomize_gcotables.py <input_file> <output_file>")
    	sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    randomize_gcotables(input_file)
