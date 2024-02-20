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


#def match_query(query, search_strings):
#    query_int = int(query, 2)
#
#    matches = []
#    for search_string in search_strings:
#        search_int = int(search_string, 2)
#        difference = query_int ^ search_int
#        match_ratio = bin(difference).count('1')
#        matches.append((search_string, match_ratio))
#
#    matches.sort(key=lambda x: x[1])
#
#    return matches

def match_query(query, search_strings):

    matches = []
    for i,search_string in enumerate(tqdm(search_strings, desc="Matching")):
        difference = query ^ search_string
        match_ratio = sum(difference)
        matches.append((i, match_ratio))

    matches.sort(key=lambda x: x[1])

    return matches

def main(holdout_file,  admixedfile, output_file):
    # Read YRIholdout, CEUholdout, and admixedfile into numpy arrays
    #holdout = read_input_file(holdout_file)
    #admixed = read_input_file(admixedfile)
    holdout=np.genfromtxt(holdout_file,dtype=int, delimiter=1, comments=None).astype(bool)
    holdout_full=np.genfromtxt(holdout_full_file,dtype=int, delimiter=1, comments=None).astype(bool)
    admixed=np.genfromtxt(admixedfile,dtype=int, delimiter=1, comments=None).astype(bool)
    search_strings = holdout.T
    query_string=admixed.T[0]
    print(query_string)
    matches = match_query(query_string, search_strings)
    best_match_index = matches[0][0]
    best_match = holdout_full[best_match_index]
    with open(output_file, "w") as output:
        output.write("\n".join([str(int(boolean)) for boolean in best_match]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the closest match in YRIholdout or CEUholdout for each range query.")
    parser.add_argument("-p","--panel", help="Path to the trimmed holdout file")
    parser.add_argument("-pf","--panel_full", help="Path to the full holdout file")
    parser.add_argument("-a","--admixedfile", help="Path to the admixedfile")
    parser.add_argument("-o","--output", help="Path to the output file")

    args = parser.parse_args()

    main(args.panel, args.admixedfile, args.output)
