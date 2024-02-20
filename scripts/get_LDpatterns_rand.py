import re
import numpy as np
import argparse
from scipy.stats import pearsonr
import pandas as pd 

class ancestryFrame():

    def __init__(self,ancestry_frame_path):
        self.df=pd.read_table(ancestry_frame_path,header=None,sep='\s+',engine='python',names=['start','stop','anc'])

    def find_anc(self,pos):
        for i,row in self.df.iterrows():
            if (row.start < pos) & (row.stop > pos):
                return(row.anc)
            else:
                continue
        return(None) 


def parse_column_phgeno(filename):
    with open(filename, 'r') as f:
        column = f.read().strip().replace('\n', '')
    return column

def parse_table_line(line):
    return list(map(int, re.split(r'\s+', line.strip())))

def parse_panel_line(line):
    return [int(char) for char in line]

def load_panel(file):
    with open(file) as f:
        return [line.strip() for line in f]

def find_correlation_index(panel, diff_row, start, end):
    max_rsq = -1
    max_index = -1
    for i in range(start, end):
        panel_row = parse_panel_line(panel[i])
        correlation, _ = pearsonr(panel_row, diff_row)
        rsq=correlation*correlation
        if rsq > max_rsq:
            max_rsq = rsq
            max_index = i
    return max_index, max_rsq

def count_subpanel_columns(panel, indexes,possible_triplets):
    subpanel = [parse_panel_line(panel[i]) for i in indexes]
    col_counts = {}
    num_columns = len(subpanel[0])

    for j in range(num_columns):
        col_str = ''.join(str(subpanel[i][j]) for i in range(len(indexes)))
        if col_str not in col_counts:
            col_counts[col_str] = 0
        col_counts[col_str] += 1


    # Create a list of counts corresponding to the entries in possible_triplets
    counts = [col_counts.get(triplet, 0) for triplet in possible_triplets]

    return counts

def main(gco_table_file, ancestry_file, panel1_file, panel0_file, recomb_phgeno_file, recomb_gco_phgeno_file, output_file):
    with open(gco_table_file) as f:
        gco_table = [line for line in f]

    panel1 = load_panel(panel1_file)
    panel0 = load_panel(panel0_file)
    recomb_phgeno=parse_column_phgeno(recomb_phgeno_file)
    recomb_gco_phgeno=parse_column_phgeno(recomb_gco_phgeno_file)
    
    print(len(panel1),len(panel0),len(recomb_phgeno),len(recomb_gco_phgeno))

    updated_gco_table = []
    #Create the list of all 8 possible triplet strings
    possible_triplets = ['{0:03b}'.format(x) for x in range(8)]
    ancestry_data=ancestryFrame(ancestry_file)
    
    for line in gco_table[1:]:
        chrom,hap,start,stop = parse_table_line(line)
        bg=ancestry_data.find_anc(start)
        panel = panel1 if bg == 1 else panel0
        try:
            recomb_range = [recomb_phgeno[i] for i in range(start, stop+1)]
            recomb_gco_range = [recomb_gco_phgeno[i] for i in range(start, stop+1)]
            differences = [i for i in range(start, stop) if recomb_phgeno[i] != recomb_gco_phgeno[i]]
        except:
            print('Error, geno files out of range - ending scan before end of gcolist')
            break

        result_row = []
        for diff_pos in differences:
            diff_row = parse_panel_line(panel[diff_pos])
            before_index,before_max_rsq = find_correlation_index(panel, diff_row, max(0, diff_pos-1050), diff_pos-50)
            after_index,after_max_rsq = find_correlation_index(panel, diff_row, diff_pos+50, min(len(panel), diff_pos+1050))
            focal_hap = ''.join([recomb_gco_phgeno[before_index], recomb_gco_phgeno[diff_pos], recomb_gco_phgeno[after_index]])
            zero_frequency = sum(1 for x in panel[diff_pos] if x == '0') / sum(1 for x in panel[diff_pos])
            counts = count_subpanel_columns(panel,[before_index,diff_pos,after_index],possible_triplets)
            result_row.append((focal_hap, zero_frequency, before_max_rsq, after_max_rsq, counts))

        updated_gco_table.append((line, result_row))

    with open(output_file, 'w') as outfile:
        outfile.write('index' + '\t' + '\t'.join(gco_table[0].strip('\n').split()) + '\tfocal_hap\tZero_Frequency\tbefore_max_rsq\tafter_max_rsq\t' + '\t'.join(possible_triplets) + '\n')
        for line, result_row in updated_gco_table:
            for focal_hap, zero_frequency, before_max_rsq, after_max_rsq, counts in result_row:
                outfile.write('\t'.join(line.strip('\n').split()) + '\t' + focal_hap + '\t' + str(zero_frequency) + '\t' + str(before_max_rsq) + '\t' + str(after_max_rsq) + '\t' + '\t'.join(map(str, counts)) + '\n')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process specified files and update gco_table")
    parser.add_argument("gco_table_file", help="Path to gco_table file")
    parser.add_argument("ancestry_file", help="Path to parsed rfmix output with suffix ; .msp.tab.column_{hap}.split.txt")
    parser.add_argument("panel0_file", help="Path to panel0 file")
    parser.add_argument("panel1_file", help="Path to panel1 file")
    parser.add_argument("recomb_phgeno_file", help="Path to recomb_phgeno file")
    parser.add_argument("recomb_gco_phgeno_file", help="Path to recomb_gco_phgeno file")
    parser.add_argument("outfile", help="Path to recomb_gco_phgeno file")

    args = parser.parse_args()

    result = main(args.gco_table_file, args.ancestry_file, args.panel1_file, args.panel0_file, args.recomb_phgeno_file, args.recomb_gco_phgeno_file,args.outfile)

