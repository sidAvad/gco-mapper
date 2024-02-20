import re

def splitbp_lines(bplines, outfile_prefix):
    '''
    Args:
    bpfile : Input admixsimu format bpfile  ( or gen file)

    Returns:
    NULL. Writes out admixSimu format bp/gen files ( one for each chromosome ) 
    '''

#{{{     
    #Loop over chromosomes
    for chrom in range(1,23):
        regex_chr = '\s'  + str(chrom) + '\|\d+\s(\d+:\d+\s)+'
        writelines=[]
        for bpline in bplines:
            match = str(re.search(regex_chr,bpline).group(0))[1:-1]
            writelines.append(match) #Get capturing group, remove trailing whitespace 
        outfile = outfile_prefix[:-3] + '_chr{}.bp'.format(chrom)
        with open(outfile,"w") as writefile:
            writefile.write("\n".join(map(str,writelines)))
#}}}
    return(writelines)


def splitbp_lines_gco(bplines, outfile_prefix):
    '''
    Args:
    bpfile : Input admixsimu format bpfile  ( or gen file)

    Returns:
    NULL. Writes out admixSimu format bp/gen files ( one for each chromosome ) 
    '''

#{{{     
    #Loop over chromosomes
    for chrom in range(1,23):
        regex_chr = '\s'  + str(chrom) + '\|\d+\s(\d+\/\d+:\d+\s)+'
        writelines=[]
        for bpline in bplines:
            match_obj=re.search(regex_chr,bpline)
            if match_obj==None:
                regex_chr_single = '(' + str(chrom) + '\|\d+)\s\s'
                match = str(re.search(regex_chr_single,bpline).group(1))
            else:
                match = str(re.search(regex_chr,bpline).group(0))[1:-1]
            writelines.append(match) #Get capturing group, remove trailing whitespace 
        outfile = outfile_prefix[:3] + '_chr{}.bp'.format(chrom)
        with open(outfile,"w") as writefile:
            writefile.write("\n".join(map(str,writelines)))
#}}}
    return(writelines)


if __name__ == '__main__':

    import argparse
    import re

    parser = argparse.ArgumentParser(description="takes in a asu.gen(bp) file and splits it into chromosomes, giving us the final admixsimu formatted files")

    # add arguments
    parser.add_argument("--inputfile", "-i", required=True, help="set input file (must be a bp file)")
    parser.add_argument("--gcomode", "-g", action='store_true')

    # read arguments from the command line
    args = parser.parse_args()
    # read static arguments 

    with open(args.inputfile,'r') as f:
        inlines = f.readlines()
    
    if not args.gcomode:
        print('regularmode')   
        splitbp_lines(inlines, args.inputfile) 
    else:
        print('gcomode')
        splitbp_lines_gco(inlines,args.inputfile)
