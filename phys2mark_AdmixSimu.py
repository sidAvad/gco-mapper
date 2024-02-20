import re
import numpy as np
import pandas as pd

def phys2mark_ancestryfile(inlines,snpDF):
    """
    Use with output of rfmix that's been modified with parse_msp.py
    """
    snpDF_pos = np.array(snpDF["pos"])
    outlines=[]
    for inline in inlines:
        physlist = [int(x) for x in re.findall(r':(\d+)', inline)]
        markers = [x + 1 for x in np.searchsorted(snpDF_pos, physlist)]
        outline = " ".join(["{}:{}".format(x[0],x[1]) for x in list(zip(haplist,markers))]) 
        outlines.append(outline)

    return(outlines)

def phys2mark_lines(inlines,snpDF):

    snpDF_pos = np.array(snpDF["pos"])
    outlines=[]
    for inline in inlines:
        haplist = re.findall(r'(\d+):', inline)
        physlist = [int(x) for x in re.findall(r':(\d+)', inline)]
        markers =[x+1 for x in np.searchsorted(snpDF_pos, physlist)] #this tags marker number behind the recombination
        outline = " ".join(["{}:{}".format(x[0],x[1]) for x in list(zip(haplist,markers))]) 
        outlines.append(outline)

    return(outlines)

def find_markers_within_gco(row,phys_pos_array):
    return(np.where((phys_pos_array>=row['start'])*(phys_pos_array<=row['end']))[0] + 1) #these tag markers within the gco

def phys2mark_gco_v2(infile,snpDF):
    '''
    finds flanking markers as well as markers within gco (the former is mostly for debuggin)
    '''
    phys_pos_array = np.array(snpDF["pos"])
    gcoDF=pd.read_table(infile)
    gcoDF['markers_within']=gcoDF[['start','end']].apply(lambda x: find_markers_within_gco(x,phys_pos_array),axis=1)
    gcoDF['start_marker']=gcoDF['markers_within'].apply(lambda x: int(x[0]) if len(x)>0 else None)
    gcoDF['end_marker']=gcoDF['markers_within'].apply(lambda x: int(x[-1]) if len(x)>0 else None)

    return(gcoDF)

def phys2mark_gco_v1(infile,snpDF):
    '''
    finds flanking markers of gcos - not very useful when using filtered markers
    because some gcos don't contain markers in them
    '''
    snpDF_pos = np.array(snpDF["pos"])
    gcoDF=pd.read_table(infile)
    gcoDF['start_marker']=gcoDF['start'].apply(lambda x: np.searchsorted(snpDF_pos, x) + 1)
    gcoDF['end_marker']=gcoDF['end'].apply(lambda x: np.searchsorted(snpDF_pos, x) + 1)

    return(gcoDF)

if __name__ == '__main__':

    import argparse
    import re

    parser = argparse.ArgumentParser(description="takes in a single chromosome asu.bp file (or gco.bp -- that is with haplotype ids instead of pops) and replaces physical position with marker number everywhere")

    # add arguments
    parser.add_argument("--inputfile", "-i", required=True, help="set input file (must be a bp file)")
    parser.add_argument("--snpfile", "-s", required=True, help="set snpfile")
    parser.add_argument("--gcomode", "-g", action='store_true')
    
    # read arguments from the command line
    args = parser.parse_args()
    # read static arguments 

        
    #if inlines[0] == "YRI CEU\n":
    #    header = inlines[0]
    #    inlines=inlines[1:]

    snpDF = pd.read_table(args.snpfile,sep=r"\s+",header=None) #This code can be finicky with the snp files depending on pandas version
    snpDF.columns = ['rsid','chr','sexavg','pos','ref','alt']
    
    writefile = args.inputfile[:-2] + 'mrk.bp'
    if args.gcomode:
        outDF = phys2mark_gco_v2(args.inputfile,snpDF)
        outDF.to_csv(writefile,sep='\t')
    else:
        with open(args.inputfile,'r') as f:
            inlines = f.readlines()
        
        outlines = phys2mark_lines(inlines, snpDF)
        with open(writefile,'w+') as f:
            #f.write(header)
            f.writelines("\n".join(outlines))
    
    
