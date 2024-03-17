import re
import sys 
import pickle
import math
from tqdm import tqdm

def find_homozygotes(inputfile,outfile):
    with open(inputfile) as f:
        inlines = f.readlines()

    genotype_list=[]
    for line in tqdm(inlines): #loop over snps
        linelist = re.findall('..',line)
        genotype_list.append([{'ind':(j+1),'genotype':str(int(elem[0]) + int(elem[1]))} for j,elem in enumerate(linelist)])
    
    with open(outfile,'wb') as f:
        pickle.dump(genotype_list,f)


if __name__ == "__main__":
    inputfile=sys.argv[1]
    outfile=sys.argv[2]

    find_homozygotes(inputfile,outfile)

