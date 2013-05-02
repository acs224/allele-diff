import string
import operator
from collections import defaultdict
from itertools import combinations

class AlleleDiff:
  def __init__(self, important_alleles, allele_file="HLAB.txt", wildtype="B*15020101"):
    self.alleles = self.parse_alleles(allele_file)
    self.important_alleles = important_alleles
    self.unimportant_alleles = [allele for allele in self.alleles.keys() if allele not in self.important_alleles]
    self.num_important = len(self.important_alleles)
    self.num_unimportant = len(self.unimportant_alleles)
    #wildtype_seq = self.alleles[wildtype].translate(string.maketrans('ACGT.', '_____'))
    self.snp_map = defaultdict(lambda: defaultdict(list))
    self.all_alleles = defaultdict(list)
    for allele, seq in self.alleles.iteritems():
      for item in seq:
        self.snp_map[item]['present'].append(allele)
      #self.snp_map[item]['absent'] = [ab_allele for ab_allele in self.alleles.keys() if ab_allele not in self.snp_map[item]['present']]
      #important_ratio = len(self.snp_map[item]['important_alleles']) / float(self.num_important)
      #self.snp_map[item]['important_ratio'] = important_ratio
      #unimportant_ratio = len(self.snp_map[item]['unimportant_alleles']) / float(self.num_unimportant)
      #self.snp_map[item]['unimportant_ratio'] = unimportant_ratio
      #print allele, 
      #self.all_alleles[allele] = self.snp_locs(wildtype_seq, self.alleles[allele])
    #print self.snp_map
    #sorted_snps = sorted([(snp, allele_map['important_ratio'] -  allele_map['unimportant_ratio']) for snp, allele_map in self.snp_map.iteritems()], key=operator.itemgetter(1))
    #print sorted_snps
    for snp, allele_map in self.snp_map.iteritems():
      important_snp_allele = [allele for allele in allele_map['present'] if allele in important_alleles]
      absent_snp_allele = [allele for allele in self.alleles.keys() if allele not in allele_map['present']]
      unimportant_snp_allele = [allele for allele in absent_snp_allele if allele in self.unimportant_alleles]
      print snp, len(important_snp_allele), (len(unimportant_snp_allele) / float(self.num_unimportant)), len(absent_snp_allele)
    
  def check_alleles(self, important):
    self.important_alleles = {}
    self.unimportant_alleles = {}
    important_allele_mapping = {}
    for allele, seq in self.all_alleles.iteritems():
      if allele in important:
        self.important_alleles[allele] = seq
      else:
        self.unimportant_alleles[allele] = seq
    
    for important_allele, seq in self.important_alleles.iteritems():
      alleles = self.test_allele_combinations(important_allele, seq)
      print important_allele, alleles
      important_allele_mapping[important_allele] = alleles
    
    return important_allele_mapping
  
  def does_not_contain_alleles(self, allele_combo, test_alleles):
    return not all(allele in test_alleles for allele in allele_combo)
    
  def test_allele_combinations(self, important_allele, seq):
    for i in range(len(seq)):
      for allele_combo in combinations(self.important_alleles[important_allele], i):
        if all(self.does_not_contain_alleles(allele_combo, v) for v in self.unimportant_alleles.values()):
          return allele_combo
  
  def parse_alleles(self, allele_file):
    allele_definitions = {}
    with open(allele_file) as f:
      for i, line in enumerate(f):
        if i < 3:
          continue
        allele_def = line.split()
        #allele_definitions[allele_def[0]] = "".join(allele_def[1:])
        allele_definitions[allele_def[0]] = ["%s:%s" % (i, base) for i, base in enumerate(allele_def[1:])]
    return allele_definitions
    
  def snp_locs(self, allele1, allele2):
    return ["%s:%s" % (i, allele2[i]) for i, base in enumerate(allele1) if base != allele2[i] ]
    
  
if __name__ == "__main__":
  
  allele_diff = AlleleDiff(["B*51240301", "B*55070101", "B*35700101", "B*07420101", "B*27090101"])
  allele_mapping = allele_diff.check_alleles(["B*51240301", "B*55070101", "B*35700101", "B*07420101", "B*27090101"])
  #print allele_diff.unimportant_alleles
