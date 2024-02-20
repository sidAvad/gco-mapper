import pandas as pd
import numpy as np
import sys
import ast
import collections


def tabulate_match(match_df,gco_df):
    match_df['gcos'] = [[]]*len(match_df)
    j=0
    current_gco=gco_df.iloc[j]
    for i,row in match_df.iterrows():
        gcolist=[]
        while (current_gco['start_marker'] <= row['end_index']+1) and (j < len(gco_df)-1):
            #print(current_gco['start_marker'],row['end'],j)
            gcolist.append(current_gco.to_dict())
            j+=1
            current_gco=gco_df.iloc[j]

        match_df.at[i, 'gcos'] = gcolist
        
        #if j == len(gco_df)-1: #Break if we reach last gco before last window (generally happens)
        #    print('last gco reached before last window')
        #    break
        
    return(match_df)

def consensus_string(strings):
    consensus = ''
    for i in range(len(strings[0])):
        column = [s[i] for s in strings]
        counter = collections.Counter(column)
        consensus += counter.most_common(1)[0][0]
    return consensus

def find_diffs(row):
    diff_markers=[row['start_index']+i+1 for i in range(len(row['admixed_string'])) if \
                    row['admixed_string'][i] != row['consensus_seq'][i]]
    gco_list=row['gcos']
    diffs_in_gco=0
    if row['true_gco_count']>0:
        for gco in gco_list:
            for diff in diff_markers:
                if (diff >= gco['start_marker']) and (diff <= gco['end_marker']):
                    diffs_in_gco+=1
                #print(gco['start_marker'],gco['end_marker'],diff,diffs_in_gco)
   
    return pd.Series([diff_markers,diffs_in_gco])

def get_consensus_seq(l_str):
    l=ast.literal_eval(l_str)
    l_s=[''.join([str(e) for e in x]) for x in l]
    consensus_ed_seq=consensus_string(l_s)
    return(consensus_ed_seq)
    
def compute_statistics(df):
    df['true_gco_count']=df['gcos'].apply(lambda x: len(x))
    df['nsnps_in_window']=df['end_index']-df['start_index']
    df['consensus_seq']=df['matches'].apply(lambda x: get_consensus_seq(x))
    df[['diff_markers','diffs_in_gco']]=df.apply(lambda x: find_diffs(x),axis=1)

    return()

def main():
    haps=int(sys.argv[1])
    window_size=int(sys.argv[2])
    print("collating results for {} haps and for window size = {}".format(haps,window_size))
    results=[]
    for hap in range(1,haps+1):
        print(hap)
        matchfile="rsync_data/matches/300_individuals/matches_wbackgroundgc_windowsize{}_hap{}.txt".format(window_size,hap)
        gcofile="rsync_data/simulated_data/300_individuals/admixed_w_backgroundgc_gcotable_{}.txt".format(hap)
        match_df=pd.read_table(matchfile)
        gco_df=pd.read_table(gcofile,index_col=0)
        match_df_joined=tabulate_match(match_df,gco_df)
        match_df_joined['exact_matches']=match_df_joined['exact_matches'].apply(lambda x: ast.literal_eval(x))
        compute_statistics(match_df_joined)
        match_df_joined['diffs_from_gcos']=match_df_joined['gcos'].apply(lambda x : [ast.literal_eval(e['diff_marker_nos']) for e in x] if len(x) > 0 else None)
        match_df_joined['gco_start_marker']=match_df_joined['gcos'].apply(lambda x : [e['start_marker'] for e in x] if len(x) > 0 else None)
        match_df_joined['gco_end_marker']=match_df_joined['gcos'].apply(lambda x : [e['end_marker'] for e in x] if len(x) > 0 else None)
        match_df_joined['hap']=hap
        results.append(match_df_joined)

    results_df=pd.concat(results)
    results_df.to_csv("results/300_individuals/results_windowsize{}.txt".format(window_size),na_rep='NULL',sep='\t')

if __name__=='__main__':
    main()

