import pandas as pd
import numpy as np
from src.sample_admixed_genotypes import *
from src.BKTree import *
import warnings


def hamming_distance2(vector1, vector2):
    assert len(vector1) == len(vector2)
    diffs=[el1 != el2 for el1, el2 in zip(vector1, vector2)]
    return sum(diffs),np.where(diffs)

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

        self.start_index = min(range(len(self.snp_df_positions)), key=lambda index: abs(self.snp_df_positions[index] - self.start))
        self.end_index = min(range(len(self.snp_df_positions)), key=lambda index: abs(self.snp_df_positions[index] - self.end))
        
        if (self.end >= self.last_pos) or (self.end_index <= self.start_index):
            raise StopIteration

        if (self.test==True) and (self.start_index > 3000):
            raise StopIteration

        YRIlines,CEUlines,ADMIXlines=[],[],[]
        for i in range(self.start_index,self.end_index): #this includes start markers index and excludes the end marker index
            YRIlines.append([int(x) for x in self.Y.readline().strip()])
            CEUlines.append([int(x) for x in self.C.readline().strip()])
            ADMIXlines.append([int(x) for x in self.A.readline().strip()])
        
        if any(len(x)==0 for x in ADMIXlines) or any(len(x)==0 for x in YRIlines) or any(len(x)==0 for x in CEUlines):
            warnings.warn("One or more lists contain an empty list; prematurely reached end of file")
        
        self.YRIgenotypes=[list(i) for i in zip(*YRIlines)] 
        self.CEUgenotypes=[list(i) for i in zip(*CEUlines)] 
        self.ADMIXgenotypes=[list(i) for i in zip(*ADMIXlines)]
        self.start=self.end+1
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
            bk_tree_insert(root_node, ['hap{}'.format(i)], [int(x) for x in genotype], hamming_distance2)
        
        # Look up nearest neighbors for each ADMIXgenotype
        self.matches=[]
        for query_vector in self.ADMIXgenotypes:
            results, dist_best, exact_matches=bk_tree_lookup3(root_node, query_vector, hamming_distance2)
            self.matches.append({\
                    'start_marker':self.start_marker,\
                    'end_marker':self.end_marker,\
                    'admixed_genotype':query_vector,\
                    'matches':[x.vector for x in results],\
                    'edit_distance':dist_best[0], \
                    'edit_positions':dist_best[1], \
                    'exact_matches':[x.vector for x in exact_matches]\
                    })

    def sub_match(self,window_size,window_index,hap_index):
            """
            This function finds approximate matches for a specific genotype within the admixed population,
            using a designated window and a BK-tree search algorithm to evaluate potential matches based on
            hamming distance. It inputs the desired window size, the window's index, and the haplotype's
            index to find and compare the target genotype against reference genotypes within the same range.

            :param window_size: The size of the window in base pairs to consider for matching.
            :param window_index: Index identifying the current window being processed.
            :param hap_index: Index of the specific haplotype within the admixed population to match.
            """
            query_vector=self.ADMIXgenotypes[hap_index][]



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
        results_allwindows=[]
        results_subwindows=[]
        for w,window in enumerate(iter(self.windows)):
            window.match()
            results_window=[]
            results_subwindow=[]
            for i,match in enumerate(window.matches): #Loop over admixed set
                if match['edit_distance']==1:
                    window.sub_match(50000,i,w)
                    results_subwindow.append({'hap_id':i+1,'window_no':w+1,'start':match['start_marker'],'end':match['end_marker']\
                    ,'admixed_string':''.join([str(x) for x in match['admixed_genotype']])\
                    ,'matches':match['matches'],'edit_distance':match['edit_distance'],'edit_positions':match['edit_positions']\
                    ,'exact_matches':match['exact_matches']})
                else:
                    results_window.append({'hap_id':i+1,'window_no':w+1,'start':match['start_marker'],'end':match['end_marker']\
                    ,'admixed_string':''.join([str(x) for x in match['admixed_genotype']])\
                    ,'matches':match['matches'],'edit_distance':match['edit_distance'],'edit_positions':match['edit_positions']\
                    ,'exact_matches':match['exact_matches']})
            
            results_allwindows.append(results_window)
            results_subwindows.append(results_subwindow)

        self.results_allwindows=[list(i) for i in zip(*results_allwindows)] #This tranposed list of lists contains an element for a admixed genotype across all windows
        self.results_subwindows=[list(i) for i in zip(*results_subwindows)] #This tranposed list of lists contains an element for a admixed genotype across all windows

    def close_data(self):
        self.Y.close()
        self.C.close()
        self.A.close()

    def write(self,outprefix):
        for i,result in enumerate(self.results_allwindows):#loop over each admixed individual
            outfile=outprefix+'_hap{}'.format(i+1)+'.txt'
            result_df=pd.DataFrame(result)
            result_df.to_csv(outfile,sep='\t',index=False)

        for i,result in enumerate(self.results_subwindows):#loop over each admixed individual
            outfile=outprefix+'_id{}'.format(i+1)+'.txt'
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
