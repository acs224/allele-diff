import unittest
import logging
from allele_diff import AlleleDiff 

class TestAlleleDiff(unittest.TestCase):
  
  def setUp(self):
    logging.basicConfig(level=logging.DEBUG)
    self.allele_diff = AlleleDiff("test_alleles.txt")
    self.all_allele_diff = AlleleDiff("HLAB.txt")
    
  def test_filter_alleles_include(self):
    filtered = self.allele_diff.filter_alleles(["B*51240301"], [{'snp': '0:A'}])
    self.assertEqual(1, len(filtered))
    
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '0:A'}])
    self.assertEqual(2, len(filtered))
    
  def test_filter_alleles_2_include(self):
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '0:A'}, {'snp': '2:C'}])
    self.assertEqual(1, len(filtered))
    
  def test_filter_alleles(self):
    filtered = self.allele_diff.filter_alleles(["B*51240301"], [{'snp': '0:A'}])
    self.assertEqual(1, len(filtered))
    
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '0:A'}])
    self.assertEqual(2, len(filtered))
    
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '11:G'}])
    self.assertEqual(2, len(filtered))
    
  def test_filter_alleles_2_exclude(self):
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '0:A'}, {'snp': '2:C'}], False)
    self.assertEqual(4, len(filtered))
    
  def test_check_1_allele(self):
     important_alleles = ['B*55070101']
     unimportant_alleles = [allele for allele in self.allele_diff.alleles.keys() if allele not in important_alleles]
     checked = self.allele_diff.check_alleles(important_alleles, unimportant_alleles)
     self.assertEqual(2, len(checked[0]['snps']))
    
  def test_check_2_alleles(self):
    important_alleles = ['B*55070101', 'B*07420101']
    unimportant_alleles = [allele for allele in self.allele_diff.alleles.keys() if allele not in important_alleles]
    checked = self.allele_diff.check_alleles(important_alleles, unimportant_alleles)
    self.assertEqual(2, len(checked))
    
  def test_check_4_alleles(self):
    important_alleles = ['B*51240301', 'B*55070101', 'B*07420101', 'B*27090101']
    unimportant_alleles = [allele for allele in self.allele_diff.alleles.keys() if allele not in important_alleles]
    checked = self.allele_diff.check_alleles(important_alleles, unimportant_alleles)
    self.assertEqual(4, len(checked))    
    
  # def test_all_check_1(self):
  #   important_alleles = ['B*55070101']
  #   unimportant_alleles = [allele for allele in self.all_allele_diff.alleles.keys() if allele not in important_alleles]
  #   checked = self.all_allele_diff.check_alleles(important_alleles, unimportant_alleles)
  #   print checked
  #   self.assertEqual(2, len(checked))
    
if __name__ == '__main__':
    unittest.main()