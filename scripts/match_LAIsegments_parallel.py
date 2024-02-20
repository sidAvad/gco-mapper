import numpy as np
import argparse
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

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
        difference = query ^ search_string
        match_ratio = sum(difference)
        matches.append((search_string, match_ratio))

    matches.sort(key=lambda x: x[1])

    return matches

def process_file(yri_panel_file, ceu_panel_file, admixed_file, output_file):
    if panel_select==0:
        holdout = read_input_file(yri_panel_file)
    else:
        holdout = read_input_file(ceu_panel_file)

    admixed = read_input_file(admixed_file)
    search_strings = holdout.T
    query_string = admixed.T[0]
    matches = match_query(query_string, search_strings)
    best_match = matches[0][0]

    with open(output_file, "w") as output:
        output.write("\n".join([str(int(boolean)) for boolean in best_match]))

def main(yri_panel_file, ceu_panel_file, admixed_files, output_files):
    with ProcessPoolExecutor() as executor:
        futures = []
        for admixed_file, output_file in zip(admixed_files, output_files):
            futures.append(executor.submit(process_file, panel_file, admixed_file, output_file))
        
        # Wait for all futures to complete
        for future in tqdm(futures, desc="Processing", total=len(futures)):
            future.result()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the closest match in YRIholdout or CEUholdout for each range query.")
    parser.add_argument("-y", "--yri_panel", help="Path to the holdout file")
    parser.add_argument("-c", "--ceu_panel", help="Path to the holdout file")
    parser.add_argument("-a", "--admixedfiles", nargs='+', help="Paths to the admixed files")
    parser.add_argument("-o", "--outputfiles", nargs='+', help="Paths to the output files")

    args = parser.parse_args()

    main(args.panel, args.admixedfiles, args.outputfiles)

