import unittest
import pandas as pd
from src.sub_match import *
from src import admixString_generators as strGen
import importlib

class TestsubScanner(unittest.TestCase):
    def setUp(self):
        self.files = {
            'YRIfile':"rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.hold.filt.phgeno",
            'CEUfile':"rsync_data/source_data/backgroundgc/OOA4J17_CEU_chr1.hold.filt.phgeno",
            'snpfile':"rsync_data/source_data/backgroundgc/OOA4J17_YRI_chr1.sim.filt.snp",
            'admixedfile':"rsync_data/simulated_data/100_individuals_testgc/admixed_w_backgroundgc_recomb_gco.phgeno",
        }
        windowDF=pd.read_table("results/100_individuals_testgc/results_windowsize200000.txt")
        sorted_filtered_windows=filter_and_centerWindows(windowDF,200000)
        self.scanner=subScanner(sorted_filtered_windows,test=True) #Test will only run till marker=3000
        self.scanner.open_data(self.files)
    
    def test_scan(self):
        self.scanner.scan()
        self.scanner.close_data()
        self.scanner.write('test_sub_match.txt')

if __name__ == '__main__':
    
    unittest.main() 