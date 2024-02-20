#eigenstrat (packed or unpacked) to vcf
#Writes to stdout, unlike many of these scripts. 
#Usage: python eigenstrat2vcf.py -r root -s snp_root [options]
#Data files are root.ind and root.geno, and a separate snp file
 
from __future__ import division, print_function
import numpy as np
import argparse

class eigenstrat():
    
    def __init__(self, file_root, snp_file_root):
        self.inds=load_ind_file(file_root)
        self.snps=load_snp_file(snp_file_root)
        self.geno=load_geno_file(file_root)
        assert len(self.snps)==len(self.geno)

def load_snp_file(file_root):

    snpList=[]
    with open(file_root + ".snp") as f:
        lines=f.readlines()

    for l in lines:
        values=l.strip("\n").split()
        keys=["ID","CHR","MAP","POS","REF","ALT"]
        snp=dict(zip(keys,values))
        snpList.append(snp)

    return(snpList)

def load_ind_file(file_root):

    with open(file_root + ".ind") as f:
        inds=[x.split('\t')[0] for x in f.read().split("\n")[:-1]]

    return(inds)

def load_geno_file(file_root):
    
    geno_array=[]
    
    with open(file_root+".geno") as f:
        lines=f.readlines()
    
    for line in lines:
        geno_array.append([int(x) for x in line.strip('\n')])

    return(geno_array)

GT_DICT={2:"0/0", 1:"0/1", 0:"1/1", 9:"./."}

################################################################################

def parse_options():
    """
    argparse
    """
    parser=argparse.ArgumentParser()
    parser.add_argument('-r', '--root', type=str, default="", help=
                        "Root for eigenstrat files - i.e {root.snp, root.geno, root.ind}")
    parser.add_argument('-s', '--snp_root', type=str, default="", help=
                        "Root for snp file - root.snp")
    return parser.parse_args()

################################################################################

def main(options):
    """
    Convert
    """

    data=eigenstrat(options.root,options.snp_root)

    #Write header. 
    print("##fileformat=VCFv4.0")
    print("##source=eig2vcf.py")
    print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">")
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"+"\t".join(data.inds))

    #Now line by line write data
    for i,s in enumerate(data.snps):
        this_snp=s
        line="\t".join([this_snp["CHR"], str(this_snp["POS"]), this_snp["ID"], 
                         this_snp["REF"], this_snp["ALT"], "100", "PASS", ".", "GT" ])
        line=line+"\t"+"\t".join([GT_DICT[x] for x in data.geno[i]])
        print(line)
        
################################################################################

if __name__=="__main__":
    options=parse_options()
    main(options)