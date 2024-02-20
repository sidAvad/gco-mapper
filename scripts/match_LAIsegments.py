import numpy as np
import argparse
from tqdm import tqdm

#def read_phgeno2numpy(infile):
#    #Define a generator function to yield integers from the file
#    def generator():
#        with open(infile, "r") as file:
#            for line in file:
#                # Process each line character by character
#                for char in line.strip():
#                    yield bool(char)
#    
#    # Use np.fromiter() to create a NumPy array from the generator
#    array = np.fromiter(generator(), dtype=bool)
#    return(array)


def read_input_file(file_path):
    # Read the input file
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Remove newline characters from each line and convert each character to a boolean
    data = [[bool(int(char)) for char in line.strip()] for line in lines]

    # Convert the nested list into a NumPy array
    array = np.array(data, dtype=bool)

    return array


def match_query(query, search_strings):

    matches = []
    for search_string in tqdm(search_strings, desc="Matching"):
        differences = query ^ search_string
        diff_positions_relative=np.where(differences)
        match_ratio = sum(differences)
        print('{}'.format(match_ratio) + '\t' + '\t'.join(map(str,diff_positions_relative[0])))
        matches.append((search_string, match_ratio))
        
    matches.sort(key=lambda x: x[1])

    return matches

def main(holdout_file,  admixedfile, output_file):
    #Read YRIholdout, CEUholdout, and admixedfile into numpy arrays
    #holdout = read_input_file(holdout_file)
    #admixed = read_input_file(admixedfile)
    search_strings=np.genfromtxt(holdout_file,dtype=int, delimiter=1, comments=None).astype(bool)
    query_string=np.genfromtxt(admixedfile,dtype=int, delimiter=1, comments=None).astype(bool)
    #search_strings = holdout.T
    #query_string=admixed.T[0]
    matches = match_query(query_string, search_strings)
    best_match,distance = matches[0]
    with open(output_file, "w") as output:
        output.write("".join([str(int(boolean)) for boolean in best_match]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the closest match in YRIholdout or CEUholdout for each range query.")
    parser.add_argument("-p","--panel", help="Path to the holdout file")
    parser.add_argument("-a","--admixedfile", help="Path to the admixedfile")
    parser.add_argument("-o","--output", help="Path to the output file")

    args = parser.parse_args()

    main(args.panel, args.admixedfile, args.output)
