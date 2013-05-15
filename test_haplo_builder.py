import unittest
import logging
from haplo_builder import HaploBuilder
from numpy import logical_and

class TestHaploBuilder(unittest.TestCase):
  
  def setUp(self):
    #logging.basicConfig(level=logging.DEBUG)
    self.haplo_diff = HaploBuilder("test_haplos.txt")
    self.all_haplo_diff = HaploBuilder("HLAB.txt")
    
  def test_filter_haplos(self):
    filtered = self.haplo_diff.filter_haps(["B*51240301"], ['0'])
    self.assertEqual(1, len(filtered))
    
    filtered = self.haplo_diff.filter_haps(self.haplo_diff.allele_index, ['0'])
    self.assertEqual(3, len(filtered))
    
  def test_filter_haplos_2(self):
    filtered = self.haplo_diff.filter_haps(self.haplo_diff.allele_index, ['0','2'])
    self.assertEqual(5, len(filtered))
    
  def test_check_1_hap(self):
     important_haplos = ['B*07420101']
     unimportant_haplos = [hap for hap in self.haplo_diff.allele_index if hap not in important_haplos]
     checked = self.haplo_diff.check_haplos(important_haplos, unimportant_haplos)
     self.assertEqual(1, len(checked[important_haplos[0]]))
     self.assertTrue('3:G' in checked[important_haplos[0]])
  
  def test_filter(self):
    haps = ['B*51240301', 'B*07420101', 'B*27090101', 'B*55070101', 'B*35700101']
    snps = ['2']
    checked = self.haplo_diff.filter_haps(haps, snps)
    self.assertEqual(3, len(checked))
    
    haps = ['B*51240301', 'B*07420101', 'B*27090101', 'B*55070101', 'B*35700101']
    snps = ['2', '0']
    checked = self.haplo_diff.filter_haps(haps, snps)
    self.assertEqual(4, len(checked))
    
    haps = ['B*15020101', 'B*35700101', 'B*51560101', 'B*07760101', 'B*50060101', 'B*08340101', 'B*52060201', 'B*08130101']
    snps = ['3', '4', '2']
    checked = self.haplo_diff.filter_haps(haps, snps)
    self.assertEqual(1, len(checked))
  
  def test_check_2_haplos(self):
    important_haplos = ['B*51240301', 'B*55070101']
    unimportant_haplos = [hap for hap in self.haplo_diff.allele_index if hap not in important_haplos]
    checked = self.haplo_diff.check_haplos(important_haplos, unimportant_haplos)
    self.assertEqual(2, len(checked))
    
  def test_check_4_haplos(self):
    important_haplos = ['B*51240301', 'B*55070101', 'B*07420101', 'B*27090101']
    unimportant_haplos = [hap for hap in self.haplo_diff.allele_index if hap not in important_haplos]
    checked = self.haplo_diff.check_haplos(important_haplos, unimportant_haplos)
    self.assertEqual(4, len(checked))    
    
  def test_check_hard_haplos(self):
     important_haplos = ['B*51240301', 'B*55070101', 'B*35700101']
     unimportant_haplos = [hap for hap in self.all_haplo_diff.allele_index if hap not in important_haplos]
     checked = self.all_haplo_diff.check_haplos(important_haplos, unimportant_haplos)
     self.assertEqual(3, len(checked))  
     self.assertEqual(3, len(checked['B*51240301']))  
     self.assertEqual(3, len(checked['B*55070101']))  
     self.assertTrue('Unidentifiable' in checked['B*35700101'])
     
  def test_all_check_1(self):
    important_haplos = ['B*51500101']
    unimportant_haplos = [hap for hap in self.all_haplo_diff.allele_index if hap not in important_haplos]
    checked = self.all_haplo_diff.check_haplos(important_haplos, unimportant_haplos)
    self.assertEqual(1, len(checked))   
    self.assertEqual(1, len(checked['B*51500101']))   
     
  def test_all_check_5(self):
    important_haplos = ['B*51500101','B*52170101','B*07760101','B*56080101','B*50060101']
    unimportant_haplos = [hap for hap in self.all_haplo_diff.allele_index if hap not in important_haplos]
    checked = self.all_haplo_diff.check_haplos(important_haplos, unimportant_haplos)
    called_snps = list(set([snp for snp_list in checked.values() for snp in snp_list]))
    self.assertEqual(5, len(checked))
    self.assertEqual(8, len(called_snps))
  
if __name__ == '__main__':
    unittest.main()