import unittest
from allele_diff import AlleleDiff 

class TestAlleleDiff(unittest.TestCase):
  
  def setUp(self):
    self.allele_diff = AlleleDiff("test_alleles.txt")
    
  def test_filter_alleles_include(self):
    filtered = self.allele_diff.filter_alleles(["B*51240301"], [{'snp': '0:A'}])
    self.assertEqual(1, len(filtered))
    
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '0:A'}])
    self.assertEqual(2, len(filtered))
    
  def test_filter_alleles_2_include(self):
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '0:A'}, {'snp': '2:C'}])
    self.assertEqual(1, len(filtered))
    
  def test_filter_alleles_exclude(self):
    filtered = self.allele_diff.filter_alleles(["B*51240301"], [{'snp': '0:A'}], False)
    self.assertEqual(0, len(filtered))
    
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '0:A'}], False)
    self.assertEqual(len(self.allele_diff.alleles.keys()) - 2, len(filtered))
    
  def test_filter_alleles_2_exclude(self):
    filtered = self.allele_diff.filter_alleles(self.allele_diff.alleles.keys(), [{'snp': '0:A'}, {'snp': '2:C'}], False)
    self.assertEqual(len(self.allele_diff.alleles.keys()) - 3, len(filtered))
    
  # def test_check_1_allele(self):
  #   important_alleles = ['B*55070101']
  #   unimportant_alleles = [allele for allele in self.allele_diff.alleles.keys() if allele not in important_alleles]
  #   checked = self.allele_diff.check_alleles2(important_alleles, unimportant_alleles)
  #   print checked
    
  def test_check_2_alleles(self):
    important_alleles = ['B*55070101', 'B*07420101', "B*27090101", 'B*51240301']
    unimportant_alleles = [allele for allele in self.allele_diff.alleles.keys() if allele not in important_alleles]
    checked = self.allele_diff.check_alleles2(important_alleles, unimportant_alleles)
    print checked
    
if __name__ == '__main__':
    unittest.main()