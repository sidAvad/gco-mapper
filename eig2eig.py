import re
import sys 

def find_homozygotes(inputfile):
    with open(inputfile) as f:
        inlines = f.readlines()

    genotype_list=[]
    for i,line in enumerate(inlines): #loop over snps
        linelist = re.findall('..',line)
        genotype_list.append([{'hap':math.ceil(i+1/2),'genotype':str(int(elem[0]) + int(elem[1]))} for elem in linelist])
    
    with open(outfile

def eig2eig(inputfile):
    '''
    Simple function to convert phased eignestrat data (.phgeno file)
    into unphased eigenstrat data (.geno file)

    '''
    with open(inputfile) as f:
        inlines = f.readlines()

    outlines = []
    for line in inlines:
        linelist = re.findall('..',line)
        newlist = [str(int(elem[0]) + int(elem[1])) for elem in linelist]
        newline = ''.join(newlist) 
        outlines.append(newline) 
    
    writefile = inputfile[:-6] + 'geno'
    with open(writefile,'w') as f:
        f.write('\n'.join(outlines))


    return()


if __name__ == "__main__":

    eig2eig(sys.argv[1])

