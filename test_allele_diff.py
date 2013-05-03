import unittest
import logging
from allele_diff import AlleleDiff 

class TestAlleleDiff(unittest.TestCase):
  
  def setUp(self):
    #logging.basicConfig(level=logging.DEBUG)
    self.haplo_diff = AlleleDiff("test_haplos.txt")
    self.all_haplo_diff = AlleleDiff("HLAB.txt")
    
  def test_filter_haplos_include(self):
    filtered = self.haplo_diff.filter_haplos(["B*51240301"], [{'snp': '0:A'}])
    self.assertEqual(1, len(filtered))
    
    filtered = self.haplo_diff.filter_haplos(self.haplo_diff.haplos.keys(), [{'snp': '0:A'}])
    self.assertEqual(2, len(filtered))
    
  def test_filter_haplos_2_include(self):
    filtered = self.haplo_diff.filter_haplos(self.haplo_diff.haplos.keys(), [{'snp': '0:A'}, {'snp': '2:C'}])
    self.assertEqual(1, len(filtered))
    
  def test_filter_haplos(self):
    filtered = self.haplo_diff.filter_haplos(["B*51240301"], [{'snp': '0:A'}])
    self.assertEqual(1, len(filtered))
    
    filtered = self.haplo_diff.filter_haplos(self.haplo_diff.haplos.keys(), [{'snp': '0:A'}])
    self.assertEqual(2, len(filtered))
    
    filtered = self.haplo_diff.filter_haplos(self.haplo_diff.haplos.keys(), [{'snp': '11:G'}])
    self.assertEqual(2, len(filtered))
    
  def test_filter_haplos_2_exclude(self):
    filtered = self.haplo_diff.filter_haplos(self.haplo_diff.haplos.keys(), [{'snp': '0:A'}, {'snp': '2:C'}], False)
    self.assertEqual(4, len(filtered))
    
  def test_check_1_allele(self):
     important_haplos = ['B*55070101']
     unimportant_haplos = [allele for allele in self.haplo_diff.haplos.keys() if allele not in important_haplos]
     checked = self.haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
     self.assertEqual(2, len(checked[0]['snps']))
    
  def test_check_2_haplos(self):
    important_haplos = ['B*55070101', 'B*07420101']
    unimportant_haplos = [allele for allele in self.haplo_diff.haplos.keys() if allele not in important_haplos]
    checked = self.haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
    self.assertEqual(2, len(checked))
    
  def test_check_4_haplos(self):
    important_haplos = ['B*51240301', 'B*55070101', 'B*07420101', 'B*27090101']
    unimportant_haplos = [allele for allele in self.haplo_diff.haplos.keys() if allele not in important_haplos]
    checked = self.haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
    self.assertEqual(4, len(checked))    
    
  def test_all_check_1(self):
    important_haplos = ['B*55070101']
    #important_haplos = ['B*51240301']
    unimportant_haplos = [allele for allele in self.all_haplo_diff.haplos.keys() if allele not in important_haplos]
    checked = self.all_haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
    self.assertEqual(1, len(checked))   
        
  def test_all_check_4(self):
    important_haplos = ['B*51240301', 'B*55070101', 'B*07420101', 'B*27090101']
    unimportant_haplos = [allele for allele in self.all_haplo_diff.haplos.keys() if allele not in important_haplos]
    checked = self.all_haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
    self.assertEqual(4, len(checked))
  
  # def test_all_check_10(self):
  #   important_haplos = ['B*51240301', 'B*55070101', 'B*35700101', 'B*07420101', 'B*27090101', 'B*51560101', 'B*07760101', 'B*50060101', 'B*08340101', 'B*52060201']
  #   unimportant_haplos = [allele for allele in self.all_haplo_diff.haplos.keys() if allele not in important_haplos]
  #   checked = self.all_haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
  #   self.assertEqual(10, len(checked))
    
  # def test_all_check(self):
  #   for important_allele in self.all_haplo_diff.haplos.keys():
  #     important_haplos = [important_allele]
  #     unimportant_haplos = [allele for allele in self.all_haplo_diff.haplos.keys() if allele not in important_haplos]
  #     checked = self.all_haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
  #     print checked
  #     if len(checked) != 1:
  #       print "ERROR", important_haplos, checked
  #     else:
  #       print "no error", important_haplos
    
if __name__ == '__main__':
    unittest.main()