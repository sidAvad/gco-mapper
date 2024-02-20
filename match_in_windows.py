import pandas as pd
import numpy as np
from src.sample_admixed_genotypes import *

def match_query(query_string, search_strings):
    
    query=np.array([int(i) for i in query_string]).astype(bool)

    matches = []
    for search_string in search_strings:
        search=np.array([int(i) for i in search_string]).astype(bool)
        differences = query ^ search
        diff_positions_relative=np.where(differences)
        match_ratio = sum(differences)
        #print('{}'.format(match_ratio) + '\t' + '\t'.join(map(str,diff_positions_relative[0])))
        matches.append((search_string, match_ratio))
        
    matches.sort(key=lambda x: x[1])

    return matches

class Scanner():
    """
    class that holds data for doing string matching and inference
    """

    @classmethod
    def read_data(cls,files):
        with open(files['YRIfile']) as f:
            cls._YRI_lines=f.readlines()
        with open(files['CEUfile']) as f:
            cls._CEU_lines=f.readlines()
        cls._snp_table=pd.read_table(files['snpfile'],header=None,\
                names=['id','chr','gen_pos','phys_pos','ref','alt'],sep=r'\s+',index_col=0)
    
    
    @classmethod 
    def scan(cls,window_size):
        """
        scan creates windows along the genome and stores all the founder haplotypes
        """
        snp_df_positions=list(cls._snp_table['phys_pos'])
        end=False
        window_start=snp_df_positions[0]
        last_pos=snp_df_positions[-1]
        cls._windows=[]
        while not end:
            window_end=window_start + window_size
            window_start_marker=np.searchsorted(snp_df_positions,window_start)
            window_end_marker=np.searchsorted(snp_df_positions,window_end)
           

            #Store all haplotypes
            if window_start_marker == 0:
                haps_YRI=[s[window_start_marker:window_end_marker] for s in cls._YRI_lines]
                haps_CEU=[s[window_start_marker:window_end_marker] for s in cls._CEU_lines]
            else:
                haps_YRI=[s[window_start_marker+1:window_end_marker] for s in cls._YRI_lines]
                haps_CEU=[s[window_start_marker+1:window_end_marker] for s in cls._CEU_lines]
            
            cls._windows.append({'start':window_start_marker,'end':window_end_marker,'haps_YRI':haps_YRI,'haps_CEU':haps_CEU})
            
            if window_end >= last_pos:
                end=True
            
            window_start=window_end+1
    
    @classmethod
    def get_windows(cls):
        return cls._windows
    
    def __init__(self,admixedstring):
        self.windows=self.__class__.get_windows()
        self.admixed_line=admixedstring 
    
    def match(self):
        """
        Creates ordered list of matches along the window
        """
        for window in self.windows:
            
            if window['start']==0:
                window['admixed_string']=self.admixed_line[window['start']:window['end']]
            else:
                window['admixed_string']=self.admixed_line[window['start']+1:window['end']]
            
            assert len(window['admixed_string']) == len(window['haps_YRI'][0])
            
            window['matches']\
                    =match_query(window['admixed_string'],window['haps_YRI']+window['haps_CEU'])
    
    def write(self,outfile):
        windows_df=pd.DataFrame(self.windows)
        windows_df.to_csv(outfile,sep='\t',columns=['start','end','admixed_string','matches'])

def main():
    parser = argparse.ArgumentParser(description="sample recomb breakpoints and splice in gcos with GCbias")
    parser.add_argument('-a', '--admixedstring_file', required=True, help="file containing admixed strings, one per line")
    parser.add_argument('-Y', '--YRIfile', required=True, help="snp file corresponding to chromosome")
    parser.add_argument('-C', '--CEUfile', required=True, help="snp file corresponding to chromosome")
    parser.add_argument('-s', '--snpfile', required=True, help="snp file corresponding to chromosome")
    parser.add_argument('-o', '--output_prefix', required=True, help="Prefix for output files")
    args = parser.parse_args()
    
    Scanner.read_data({'YRIfile':args.YRIfile, 'CEUfile':args.CEUfile, 'snpfile':args.snpfile})
    Scanner.scan(50000)
   
    #with open(args.admixedstring_file) as f:
    #    i=1
    #    for admixstring in f:
    #        scanner = Scanner(admixstring)
    #        scanner.match()
    #        outfile=args.output_prefix + '_hap{}.txt'.format(i)
    #        scanner.write(outfile)
    #        i+=1
    
    with open(args.admixedstring_file) as f:
        admixstring=f.readline()

    scanner = Scanner(admixstring)
    scanner.match()
    outfile=args.output_prefix
    scanner.write(outfile)
    
if __name__ == "__main__":
    main()
