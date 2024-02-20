import unittest
import pandas as pd
from src.match_in_windows import *
from src import admixString_generators as strGen
#import importlib

class TestHap(unittest.TestCase):
    def setUp(self):
        self.files = {
            'YRIfile': "rsync_data/source_data/OutOfAfrica_4J17_YRI_chr1.sim.filt.phgeno.trnsp",
            'CEUfile': "rsync_data/source_data/OutOfAfrica_4J17_CEU_chr1.sim.filt.phgeno.trnsp",
            'snpfile': "rsync_data/source_data/OutOfAfrica_4J17_YRI_chr1.sim.filt.snp",
        }
        Scanner.read_data(self.files)
        Scanner.scan(50000)
        with open('testadmixed.txt') as f:
            admixed_line=f.readlines()[0]
        self.scanner=Scanner(admixed_line)

    def test_match(self):
        self.scanner.match()
        self.scanner.write('matchout.txt') 
    
if __name__ == '__main__':
    
    unittest.main() 
