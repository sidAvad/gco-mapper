import unittest
import pandas as pd
from src.match_in_windows2 import *
from src import admixString_generators as strGen
import importlib

class TestWindows(unittest.TestCase):
    def setUp(self):
        self.files = {
            'YRIfile': "rsync_data/source_data/OutOfAfrica_4J17_YRI_chr1.sim.filt.phgeno",
            'CEUfile': "rsync_data/source_data/OutOfAfrica_4J17_CEU_chr1.sim.filt.phgeno",
            'admixedfile': "testadmixed_all.phgeno",
            'snpfile': "rsync_data/source_data/OutOfAfrica_4J17_YRI_chr1.sim.filt.snp",
        }
        self.scanner=Scanner(50000,test=True) #Test will only run till marker=3000
        self.scanner.open_data(self.files)
    
    def test_scan(self):
        self.scanner.scan()
        self.scanner.close_data()
        self.scanner.write('matchout')

if __name__ == '__main__':
    
    unittest.main() 
    
