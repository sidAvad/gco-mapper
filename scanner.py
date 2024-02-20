import argparse
import pandas as pd
import numpy as np

class Scanner:
    def __init__(self, input_file, control_gcos_file, snp_file):
        
        self.input_df = pd.read_csv(input_file, sep=r'\s+')
        self.control_gcos_df = pd.read_csv(control_gcos_file, sep=r'\s+')
        required_columns = ['hap', 'chrom', 'pop', 'first', 'last', 'diffs', 'bg']
        if not all(col in self.input_df.columns for col in required_columns):
            raise ValueError("Input data does not have all required columns.")
        
        self.snp_df = pd.read_table(snp_file, sep=r'\s+',index_col=0,header=None,\
                names=['chr','gen_pos','phys_pos','ref','alt'])
        self.snp_df_positions=self.snp_df['phys_pos'].tolist()
        del self.snp_df

        self.window_df=pd.DataFrame()

    #Create windows around gene-conversion one-by-one
    def create_windows_around_gcos(self,window_length):
        """
        Create a dataframe with windows around gcos of variable length
        """
        #Get all windows around gcos 
        window_df_list_of_dicts=[]
        i=0
        while i < len(self.input_df)-1:
            row=self.input_df.iloc[i]
            first_physical=self.snp_df_positions[row['first']]
            last_physical=self.snp_df_positions[row['last']]
            gco_mid=(first_physical + last_physical)//2
            window_start,window_end=gco_mid-(window_length//2),gco_mid+(window_length//2)
            window_start_marker=np.searchsorted(self.snp_df_positions,window_start)
            window_end_marker=np.searchsorted(self.snp_df_positions,window_end)
            if last_physical > window_end:
                print('gco is longer than window')
                i=i+1
                continue

            #Create window
            row_dict={'hap':row['hap'],'chrom':row['chrom'],'window_start':window_start,\
                    'window_end':window_end,'window_start_marker':window_start_marker,\
                    'window_end_marker':window_end_marker,
                    'gcos':None}

            
            #print(i,first_physical,window_end) 
            #Loop over gcos that are within the window and add them to a list of dicts 
            gco_dicts=[]
            while (first_physical < window_end):
                #print(i)
                #print('within while loop i = {} , length of input = {}, first physical = {}\
                #        and window end = {}'.format(i,len(self.input_df),first_physical,window_end))
                #Store gene-conversion information in dictionaries
                gco_dicts.append({'pop':row['pop'],'first':self.snp_df_positions[row['first']]\
                        ,'first_marker':row['first'],'last':self.snp_df_positions[row['last']]\
                        ,'last_marker':row['last'],'diffs':row['diffs'],'bg':row['bg']})
                #update within-window loop only if not at the ned
                if i >= len(self.input_df)-1:
                    #print('breaking because i ({}) is less than equal to length of input minus 1  {}'.format(i,len(self.input_df)-1))
                    break
                else:
                    i=i+1
                    #print('i after increment = {}'.format(i))
                    row=self.input_df.iloc[i] 
                    first_physical=self.snp_df_positions[row['first']]
            
            #Store all gco information in row_dict and append it to df_list_of_dicts
            row_dict['gcos']=gco_dicts
            window_df_list_of_dicts.append(row_dict)

        #Create window df
        self.window_df=pd.DataFrame(window_df_list_of_dicts,\
                columns=['chrom','hap','window_start','window_end','window_start_marker','window_end_marker','gcos'])



def main():
    parser = argparse.ArgumentParser(description="Process scanning file.")
    parser.add_argument('-i', '--inputfile', required=True, help="Input file path.")
    parser.add_argument('-c', '--control_gcos_file', required=True, help="Input file path.")
    parser.add_argument('-s', '--snpfile', required=True, help="Map file path.")
    parser.add_argument('-w', '--windowlength', type=int, required=True, help="Map file path.")
    parser.add_argument('-o', '--output_file', required=True, help="Output file path.")
    args = parser.parse_args()

    scanner = Scanner(args.inputfile, args.control_gcos_file, args.snpfile)
    scanner.create_windows_around_gcos(args.windowlength)
    scanner.add_control_gcos()
    #Save the processed data to the output file
    scanner.window_df.to_csv(args.output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()
