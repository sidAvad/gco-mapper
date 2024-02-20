import pandas as pd
import numpy as np
from src import admixString_generators as admixStrGen
import re
import os
import argparse
import random

def keep_gco(gco_string,background_string,refs,alts,gcbias,coding_dict={'ref':'1','alt':'0'}):
    gc_gco=[refs[i] if x==coding_dict['ref'] else alts[i] for i,x in enumerate(gco_string)]
    gc_background=[refs[i] if x==coding_dict['ref'] else alts[i] for i,x in enumerate(background_string)]
   
    #Roll 70% dice if gco > background, else roll 30% dice
    if sum([(x=='G') or (x=='C') for x in gc_gco]) >= sum([(x=='G') or (x=='C') for x in gc_background]):
        if random.random() < gcbias:
            print('retaining gco : g/c content of transmitted allele >= g/c content of background')
            return(True)
        else:
            print('dropping gco : g/c content of transmitted allele >= g/c content of background')
            return(False)
    else:
        if random.random() < (1-gcbias):
            print('retaining gco : g/c content of transmitted allele < g/c content of background')
            return(True)
        else:
            print('dropping gco : g/c content of transmitted allele < g/c content of background')
            return(False)

    
class Hap():
    """
    base class that contains bpfiles, and holds the simulated genotype founder data
    contains methods for formatting bpfiles, sampling segments based on recombination bpfile
    and splicing in gene-conversions accounting for GCbias
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
    def get_data(cls):
        return cls._YRI_lines,cls._CEU_lines,cls._snp_table
    
    def __init__(self,bpstr,gcotable):
        self.bpstr=bpstr[:-1]
        self.gcotable=gcotable
        self.YRI_lines,self.CEU_lines,self.snp_table=self.__class__.get_data()
    
    def generate_segments(self,admixstring):
        """
        generate information necessary to sample segments
        Pull out bpline and parse it into list
        of dictionaries ; create YRI and CEU map from admixture parameters
        """
        haplist = re.findall(r'(\d):', self.bpstr)
        markers = [int(x) for x in re.findall(r':(\d+)', self.bpstr)]
        assert len(haplist)==len(markers)
        
        self.YRI_map\
                =[x for y in \
                [(str(2*(i+1)-1),str(2*(i+1))) for i in range(len(admixstring)) if admixstring[i]=='0']\
                for x in y ] #stores list of haplotypes that are YRI
        self.CEU_map\
                =[x for y in \
                [(str(2*(i+1)-1),str(2*(i+1))) for i in range(len(admixstring)) if admixstring[i]=='1']\
                for x in y ] #stores list of haplotypes that are CEU

        self.segments=[{'hap':str(int(haplist[i])+1),'end_m':markers[i]} for i in range(len(haplist))]
    
    def sample(self):
        """
        Sample admixed genotype after generating segments
        """
        if not self.segments:
            print("error: please generate_segments() first")
            exit()
        self.admix_line="" 
        start=0
        for segment in self.segments:
            #sample from appropriate population and individual
            if segment['hap'] in self.YRI_map:
                sample_id=self.YRI_map.index(segment['hap'])
                write_string=self.YRI_lines[sample_id][start:segment['end_m']]
                self.admix_line+=write_string
                start=segment['end_m'] #this is equivalent to the index of the marker ahead
            elif segment['hap'] in self.CEU_map:
                sample_id=self.CEU_map.index(segment['hap'])
                write_string=self.CEU_lines[sample_id][start:segment['end_m']]
                self.admix_line+=write_string
                start=segment['end_m'] #this is equivalent to the index of the marker ahead
            else:
                print(segment['hap'])
                print("hap id is not in YRI or CEU map")
                exit()
        print(len(self.admix_line))
        print(len(self.YRI_lines[0]))
        assert len(self.admix_line) == len(self.YRI_lines[0])
        return()

    def splice(self,gcbias):
        """
        splice in gene-conversions after generating admixed recombinant genotype
        """
        refs=self.snp_table['ref'].tolist()
        alts=self.snp_table['alt'].tolist()
        if not self.admix_line:
            print("please run .sample(admixstring) fist")
            exit()
        self.admix_line_w_gcos=self.admix_line
        indices_to_delete=[]
        self.gcotable['diff_marker_nos'] = None
        for index,gco in self.gcotable.iterrows():
            hap_id=str(gco['hap'])
            if hap_id in self.YRI_map:
                sample_id=self.YRI_map.index(hap_id)
                write_string=self.YRI_lines[sample_id][int(gco['start_marker'])-1:int(gco['end_marker'])] #this way of indexing ensures start and end marker inclusive
                replace_string=self.admix_line_w_gcos[int(gco['start_marker'])-1:int(gco['end_marker'])]
                if keep_gco(write_string,replace_string,refs,alts,gcbias):
                    diff_marker_nos=[int(gco['start_marker'])+i for i in range(len(write_string)) \
                            if write_string[i]!=replace_string[i]]
                    self.gcotable.at[index,'diff_marker_nos']=diff_marker_nos
                    self.admix_line_w_gcos=self.admix_line_w_gcos[:int(gco['start_marker'])-1] \
                            + write_string \
                            + self.admix_line_w_gcos[int(gco['end_marker']):]
                else:
                    indices_to_delete.append(index)
                    continue
            
            elif hap_id in self.CEU_map:
                sample_id=self.CEU_map.index(hap_id)
                write_string=self.CEU_lines[sample_id][int(gco['start_marker'])-1:int(gco['end_marker'])]
                replace_string=self.admix_line_w_gcos[int(gco['start_marker'])-1:int(gco['end_marker'])]
                if keep_gco(write_string,replace_string,refs,alts,gcbias):
                    diff_marker_nos=[int(gco['start_marker'])+i for i in range(len(write_string)) \
                            if write_string[i]!=replace_string[i]]
                    self.gcotable.at[index,'diff_marker_nos']=diff_marker_nos
                    self.admix_line_w_gcos=self.admix_line_w_gcos[:int(gco['start_marker'])-1] \
                            + write_string \
                            + self.admix_line_w_gcos[int(gco['end_marker']):]
                else:
                    indices_to_delete.append(index)
                    continue
            
            else:
                print(hap_id)
                print("hap id is not in YRI or CEU map")
                exit()
        
        self.gcotable=self.gcotable.drop(indices_to_delete)
        return()
    
    def write_final(self,outfile):
        with open(outfile,'a+') as f:
            f.writelines(self.admix_line_w_gcos)
    
    def write_recomb_only(self,outfile2):
        with open(outfile2,'a+') as f:
            f.writelines(self.admix_line)
    
    def write_gcotable(self,gco_outfile):
        self.gcotable.drop(self.gcotable.filter(regex="Unname"),axis=1, inplace=True)
        self.gcotable.to_csv(gco_outfile,sep='\t')

def main():
    parser = argparse.ArgumentParser(description="sample recomb breakpoints and splice in gcos with GCbias")
    parser.add_argument('-b', '--recomb_bpfile', required=True, help="recombination bp file \
            one string for each individual separated by newline") 
    parser.add_argument('-g', '--gco_bp_file_list', required=True, help="gco bp file (tabular)")
    parser.add_argument('-y', '--YRIfile', required=True, help="snp file corresponding to chromosome")
    parser.add_argument('-c', '--CEUfile', required=True, help="snp file corresponding to chromosome")
    parser.add_argument('-s', '--snpfile', required=True, help="snp file corresponding to chromosome")
    parser.add_argument('-o', '--outprefix', required=True, help="output prefix for recomb only, recomb with gco")
    args = parser.parse_args()
    
    Hap.read_data({'YRIfile':args.YRIfile, 'CEUfile':args.CEUfile, 'snpfile':args.snpfile})
    
    with open(args.recomb_bpfile) as f:
        bplines=f.readlines()
    
    gcotables=[] #Store all the gcotables
    with open(args.gco_bp_file_list) as f:
        for filename in f:
            gcotables.append(pd.read_table(filename[:-1]))
   
    assert len(gcotables)==len(bplines)
    
    admixstr=admixStrGen.generateAdmixStrings(6,6,0.8,0.8)[0] #will generate admixstring deterministically unless seed argument is specified 
    #Iterate over bpstrs and gcotables simultaneously ; sample and splice
    for i,bpstr in enumerate(bplines): #TODO:Refactor this to use with construct for writing
        hap = Hap(bpstr,gcotables[i])
        hap.generate_segments(admixstr)
        hap.sample()
        hap.splice(0.7)
        hap.write_recomb_only(args.outprefix+'_recomb.txt')
        hap.write_gcotable(args.outprefix+'_gcotable_{}.txt'.format(i+1))
        hap.write_final(args.outprefix+'_recomb_gco.txt')
        
if __name__ == "__main__":
    main()


