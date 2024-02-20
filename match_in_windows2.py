import pandas as pd
import numpy as np
from src.sample_admixed_genotypes import *
from src.BKTree import *
import warnings

def hamming_distance(vector1, vector2):
    assert len(vector1) == len(vector2)
    return sum(el1 != el2 for el1, el2 in zip(vector1, vector2))

class windows():
    """
    iterator that takes in open file objects and steps through windows
    """

    def __init__(self,Y,C,A,snpDF,window_size,test):
        self.Y=Y
        self.C=C
        self.A=A
        self.snpDF=snpDF
        self.window_size=window_size
        self.snp_df_positions=list(self.snpDF['phys_pos'])
        self.test=test

    def __iter__(self):
    
        self.last_pos=self.snp_df_positions[-1]
        
        self.start=self.snp_df_positions[0]
        self.end=self.start + self.window_size
        
        #self.start_marker=np.searchsorted(self.snp_df_positions,self.start)
        #self.end_marker=np.searchsorted(self.snp_df_positions,self.end)
        
        return self

    def __next__(self):
        """
        return iterated state - updates state to markers in current range, but next start and end positions
        """

        self.start_index=np.searchsorted(self.snp_df_positions,self.start)
        self.end_index=np.searchsorted(self.snp_df_positions,self.end)
        
        if (self.end >= self.last_pos) or (self.end_index <= self.start_index):
            raise StopIteration

        if (self.test==True) and (self.start_index > 30000):
            raise StopIteration

        YRIlines,CEUlines,ADMIXlines=[],[],[]
        for i in range(self.start_index,self.end_index): #it is important not to read the end marker itself to avoid overlaps
            YRIlines.append([int(x) for x in self.Y.readline().strip()])
            CEUlines.append([int(x) for x in self.C.readline().strip()])
            ADMIXlines.append([int(x) for x in self.A.readline().strip()])
        
        if any(len(x)==0 for x in ADMIXlines) or any(len(x)==0 for x in YRIlines) or any(len(x)==0 for x in CEUlines):
            warnings.warn("One or more lists contain an empty list; prematurely reached end of file")
        
        self.YRIgenotypes=[list(i) for i in zip(*YRIlines)] 
        self.CEUgenotypes=[list(i) for i in zip(*CEUlines)] 
        self.ADMIXgenotypes=[list(i) for i in zip(*ADMIXlines)]
        self.start=self.end
        self.end=self.start + self.window_size
        
        return self

    def match(self):
        """
        Creates matches along the window, one for each query genotype in the admixed set 
        """
        # Initialize an empty BKTreeNode
        root_node = BKTreeNode.make_empty()

        # Insert vectors into the BK tree
        for i,genotype in enumerate(self.YRIgenotypes + self.CEUgenotypes):
            bk_tree_insert(root_node, ['hap{}'.format(i)], [int(x) for x in genotype], hamming_distance)
        
        # Look up nearest neighbors for each ADMIXgenotype
        self.matches=[]
        for query_vector in self.ADMIXgenotypes:
            results, dist_best, exact_matches=bk_tree_lookup2(root_node, query_vector, hamming_distance)
            self.matches.append({\
                    'start_index':self.start_index,\
                    'end_index':self.end_index,\
                    'admixed_genotype':query_vector,\
                    'matches':[x.vector for x in results],\
                    'edit_distance':dist_best, \
                    'exact_matches':[x.vector for x in exact_matches]\
                    })
class Scanner():
    """
    class that holds data for doing string matching and inference
    """
    def __init__(self,window_size,test=False):
        self.window_size=window_size
        self.test=test
        pass

    def open_data(self,files):
        self.Y=open(files['YRIfile'],'r')
        self.C=open(files['CEUfile'],'r')
        self.A=open(files['admixedfile'],'r')
        self.snpDF=pd.read_table(files['snpfile'],header=None,\
                names=['id','chr','gen_pos','phys_pos','ref','alt'],sep=r'\s+',index_col=0)
        return()

    def scan(self):
        self.windows=windows(self.Y,self.C,self.A,self.snpDF,self.window_size,self.test)
        results=[]
        for window in iter(self.windows):
            window.match()
            results_window=[]
            for match in window.matches: #Loop over admixed set
                results_window.append({'start_index':match['start_index'],'end_index':match['end_index']\
                ,'admixed_string':''.join([str(x) for x in match['admixed_genotype']])\
                ,'matches':match['matches'],'edit_distance':match['edit_distance'],'exact_matches':match['exact_matches']})
            results.append(results_window)

        self.results=[list(i) for i in zip(*results)] #This tranposed list of lists contains an element for a admixed genotype across all windows

    def close_data(self):
        self.Y.close()
        self.C.close()
        self.A.close()

    def write(self,outprefix):
        for i,result in enumerate(self.results):#loop over each admixed individual
            outfile=outprefix+'_hap{}'.format(i+1)+'.txt'
            result_df=pd.DataFrame(result)
            result_df.to_csv(outfile,sep='\t',index=False)

def main():
    parser = argparse.ArgumentParser(description="match admixed haplotypes with reference panel in windows using BK Trees")
    parser.add_argument('-a', '--admixed_file', required=True, help="file containing admixed strings, one per line")
    parser.add_argument('-w', '--window_size', type=int, required=True, help="window size to run matching with")
    parser.add_argument('-Y', '--YRIfile', required=True, help="phgeno file corresponding to reference YRI")
    parser.add_argument('-C', '--CEUfile', required=True, help="phgeno file corresponding to reference CEU")
    parser.add_argument('-s', '--snpfile', required=True, help="snp file corresponding to chromosome")
    parser.add_argument('-o', '--output_prefix', required=True, help="Prefix for output files")
    args = parser.parse_args()
    
    scanner=Scanner(args.window_size)
    scanner.open_data({'YRIfile':args.YRIfile, 'CEUfile':args.CEUfile, 'admixedfile':args.admixed_file, 'snpfile':args.snpfile})
    scanner.scan()
    scanner.close_data()
    scanner.write(args.output_prefix)
    
if __name__ == "__main__":
    main()
