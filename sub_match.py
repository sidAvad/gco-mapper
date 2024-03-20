import pandas as pd
from pandas import DataFrame
import numpy as np
from src.sample_admixed_genotypes import *
from src.match_in_windows2 import *
import warnings
from typing import *
import ast
from tqdm import tqdm


def filter_and_centerWindows(windowDF:DataFrame,window_size:int) -> DataFrame:
    windowDF['exact_matches']=windowDF['exact_matches'].apply(lambda x: ast.literal_eval(x))
    windowDF['diff_markers']=windowDF['diff_markers'].apply(lambda x: ast.literal_eval(x))
    sorted_filtered_results=windowDF[((windowDF.edit_distance == 1) | (windowDF.edit_distance == 2)) & (windowDF.exact_matches.map(len)==0) & (windowDF.diff_markers.map(len)>0)].sort_values(by=['start_index'])
    sorted_filtered_results['subwindow_start_index']=sorted_filtered_results['diff_markers'].apply(lambda x: x[0]-(window_size//4)-1)
    sorted_filtered_results['subwindow_end_index']=sorted_filtered_results['diff_markers'].apply(lambda x: x[0]+(window_size//4)-1)
    return(sorted_filtered_results)

class subScanner(Scanner):
    
    def __init__(self, sorted_filtered_windows:DataFrame, test=False):
        if test==True:
            self.windows = sorted_filtered_windows.iloc[:50]
        else:
            self.windows = sorted_filtered_windows
    
    def scan(self):
        for start_index,end_index in tqdm(zip(self.windows['subwindow_start_index'].unique(),self.windows['subwindow_end_index'].unique())): #loop over all windows that begin in a particular position
            haplist=self.windows[self.windows.subwindow_start_index==start_index].hap.tolist()
            YRIlines,CEUlines,ADMIXlines=[],[],[]

            #find the appropriate line to read from based on constructed window
            i=0
            while i < start_index:
                self.Y.readline()
                self.C.readline()
                self.A.readline()
                i+=1

            while (i >= start_index) and (i <= end_index):
                print(i,start_index,end_index)
                YRIlines.append([int(x) for x in self.Y.readline().strip()])
                CEUlines.append([int(x) for x in self.C.readline().strip()])
                ADMIXlines.append([int(x) for x in self.A.readline().strip()])
                i+=1 
            
            if any(len(x)==0 for x in ADMIXlines) or any(len(x)==0 for x in YRIlines) or any(len(x)==0 for x in CEUlines):
                warnings.warn("One or more lists contain an empty list; prematurely reached end of file")
            
            #convert to genotypes            
            YRIgenotypes=[list(i) for i in zip(*YRIlines)] 
            CEUgenotypes=[list(i) for i in zip(*CEUlines)] 
            ADMIXgenotypes=[list(i) for i in zip(*ADMIXlines)]

            root_node = BKTreeNode.make_empty()
            # Insert vectors into the BK tree
            for i,genotype in enumerate(YRIgenotypes + CEUgenotypes):
                bk_tree_insert(root_node, ['hap{}'.format(i)], [int(x) for x in genotype], hamming_distance)
        
            # Look up nearest neighbors for each ADMIXgenotype
            self.matches=[]
            for query_vector in [ADMIXgenotypes[i-1] for i in haplist]: #match only for the haplotypes in the filtered windows (at edit distance 1/2)
                results, dist_best, exact_matches=bk_tree_lookup2(root_node, query_vector, hamming_distance)
                self.matches.append({\
                        'start_index':start_index,\
                        'end_index':end_index,\
                        'admixed_genotype':query_vector,\
                        'matches':[x.vector for x in results],\
                        'edit_distance':dist_best, \
                        'exact_matches':[x.vector for x in exact_matches]\
                        })


def main():
    parser = argparse.ArgumentParser(description="match admixed haplotypes with reference panel in windows using BK Trees")
    parser.add_argument('-W', '--windowDF_path', required=True, help="dataframe containing window matches")
    parser.add_argument('-w', '--window_size', type=int, required=True, help="window size to run matching with")
    parser.add_argument('-a', '--admixed_file', required=True, help="file containing admixed strings, one per line")
    parser.add_argument('-Y', '--YRIfile', required=True, help="phgeno file corresponding to reference YRI")
    parser.add_argument('-C', '--CEUfile', required=True, help="phgeno file corresponding to reference CEU")
    parser.add_argument('-s', '--snpfile', required=True, help="snp file corresponding to chromosome")
    parser.add_argument('-o', '--output_prefix', required=True, help="Prefix for output files")
    args = parser.parse_args()
    
    windowDF=pd.read_table(args.windowDF_path)
    sorted_filtered_windows=filter_and_centerWindows(windowDF,args.window_size)
    scanner=subScanner(args.window_size,sorted_filtered_windows)
    scanner.open_data({'YRIfile':args.YRIfile, 'CEUfile':args.CEUfile, 'admixedfile':args.admixed_file, 'snpfile':args.snpfile})
    scanner.scan()
    scanner.close_data()
    scanner.write(args.output_prefix)

if __name__ == "__main__":
    main()
