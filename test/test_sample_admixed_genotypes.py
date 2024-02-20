import unittest
import pandas as pd
from src.sample_admixed_genotypes import *
from src import admixString_generators as strGen
import importlib

class TestHap(unittest.TestCase):
    def setUp(self):
        self.files = {
            'gcofile': "rsync_data/bpfiles/gcos_1_0.gco_chr1.mrk.bp",
            'YRIfile': "rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.phgeno.trnsp",
            'CEUfile': "rsync_data/source_data/backgroundgc/OOA4J17_CEU_chr1.sim.filt.phgeno.trnsp",
            'snpfile': "rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp"
        }
        
        Hap.read_data(self.files)
        
        with open("rsync_data/bpfiles/gcos_1_0.recomb_chr1.mrk.bp") as f:
            bpstr=f.readline()
        gcotable=pd.read_table(self.files['gcofile'])
        self.hap = Hap(bpstr,gcotable)      
        self.admixstr=strGen.generateAdmixStrings(6,6,0.8,0.8)[0]
        self.write_string_test='0'*7
        self.replace_string_test='1'*7
        self.refs_test=['G']*7
        self.alts_test=['A']*7

    def test_keep_gco(self):
        """
        Test expected proportions when write and read strings are all 1's and 0's
        """
        self.result_keep_gco=[]
        for run in range(1000):
            self.result_keep_gco.append(keep_gco(self.write_string_test,self.replace_string_test,self.refs_test,self.alts_test,gcbias=0.7))
        print(sum(self.result_keep_gco)*100/len(self.result_keep_gco))

    def test_read_data(self):
        #TODO:Implement useful tests - check column names for snpdf, check nlines for YRI and CEU are equal and print number
        assert isinstance(self.hap.bpstr, str)
        assert isinstance(self.hap.gcotable, pd.DataFrame)
        assert isinstance(self.hap.YRI_lines, list)
        assert isinstance(self.hap.CEU_lines, list)
        assert isinstance(self.hap.snp_table, pd.DataFrame)
    
    def test_write(self):
        self.hap.generate_segments(self.admixstr)
        self.hap.sample()
        self.hap.splice(0.7)
        self.hap.write_recomb_only('testrecomb.txt')
        self.hap.write_final('testrecomb_gco.txt')
        self.hap.write_gcotable('test_gcotable.txt')

if __name__ == '__main__':
    
    #importlib.reload(samp);
    #samp.hap.read_data(files);
    #hap1=samp.hap(bpstr);
    #admixstr=strGen.generateAdmixStrings(6,6,0.8,0.8)[0];
    #hap1.generate_segments(admixstr);
    #hap1.sample();
    unittest.main() 
