import string
import operator
from collections import defaultdict
from itertools import combinations
from operator import itemgetter

class AlleleDiff:
  def __init__(self, allele_file="HLAB.txt", wildtype="B*15020101"):
    self.alleles = self.parse_alleles(allele_file, wildtype)
 
  def filter_alleles(self, alleles, snp_filter, include=True):
    if len(snp_filter) < 1:
      return alleles
    filtered_alleles = []
    if not include:
      filtered_alleles = alleles
    for allele in alleles:
      seq = self.alleles[allele]
      to_filter = [snp['snp'] for snp in snp_filter]
      if include and all([snp in seq for snp in to_filter]):
        filtered_alleles.append(allele)
      elif not include and any([snp in seq for snp in to_filter]):
        filtered_alleles.remove(allele)
        
    return filtered_alleles
    
  # def check_alleles(self, important, unimportant, snp_filter=[]):
  #   
  #   print important, snp_filter
  #   important_alleles = self.filter_alleles(important, snp_filter)
  #   remaining_important = [allele for allele in important if allele not in important_alleles]
  #   if not important or not important_alleles:
  #     return []
  #   filtered_snps = [temp_snp['snp'] for temp_snp in snp_filter]
  #   
  #   print important_alleles, remaining_important
  #   unimportant_alleles = self.filter_alleles(unimportant, snp_filter, False)
  #   if len(snp_filter) > 1 and len(important_alleles) == 1:
  #     return [{'allele': important_alleles[0], 'snps': filtered_snps}] + self.check_alleles(remaining_important, unimportant)
  #     
  #   #print "FILTER:", len(important), snp_filter
  #   
  #   all_alleles = important_alleles + unimportant_alleles
  #   #print "IMPORTANT:", len(important_alleles), len(unimportant_alleles), len(all_alleles)
  #   num_important = len(important_alleles)
  #   num_unimportant = len(unimportant_alleles)
  #   snp_map = defaultdict(lambda: defaultdict(list))
  #   for allele in all_alleles:
  #     seq = self.alleles[allele]
  #     for item in seq:
  #       snp_map[item]['present'].append(allele)
  #       
  #   max_snp = None
  #   for snp, allele_map in snp_map.iteritems():
  #     important_snp_allele = [allele for allele in allele_map['present'] if allele in important_alleles]
  #     unimportant_snp_allele = [allele for allele in  allele_map['present'] if allele in unimportant_alleles]
  #     #print num_important, num_unimportant
  #     important_ratio = len(important_snp_allele) / float(num_important)
  #     unimportant_ratio = len(unimportant_snp_allele) / float(num_unimportant)
  #     raw_score = important_ratio - unimportant_ratio
  #     #print snp, raw_score
  #     #print snp, len(allele_map['present'])
  #     #print len(important_alleles), len(important_snp_allele), len(unimportant_alleles), len(unimportant_snp_allele)
  #     if '_' not in snp and '|' not in snp:
  #       # if snp not in filtered_snps and raw_score == 1:
  #       #   print "Found defining allele: %s %s %s" % (snp, important_ratio, unimportant_ratio)
  #       #   new_important_alleles = [allele for allele in important_alleles if allele not in important_snp_allele]
  #       #   filtered_snps.append(snp)
  #       #   return [{'allele': important_snp_allele[0], 'snps': filtered_snps}
  #       if not max_snp or (snp not in filtered_snps and max_snp['score'] < raw_score):
  #         max_snp = {'snp': snp, 'score': raw_score}
  #         #print snp, len(allele_map['present'])
  #         #print len(important_alleles), len(important_snp_allele), len(unimportant_alleles), len(unimportant_snp_allele)
  #         #print important_ratio, unimportant_ratio
  #         #print max_snp
  #       #print "MAX: ", max_snp
  #   snp_filter.append(max_snp)
  #   return self.check_alleles(important_alleles, unimportant, snp_filter) + self.check_alleles(remaining_important, unimportant)
  
  def check_alleles2(self, important, unimportant, snp_filter=[]):
    
    print important, "FILTER", snp_filter
    important_alleles = self.filter_alleles(important, snp_filter)
    unimportant_alleles = self.filter_alleles(unimportant, snp_filter, False)
    remaining_important = [allele for allele in important if allele not in important_alleles]
    filtered_snps = [temp_snp['snp'] for temp_snp in snp_filter]
    if not important or len(important_alleles) == 0:
      return [] 
    if len(snp_filter) > 1 and len(important_alleles) == 1 and len(unimportant_alleles) == 0:
      return [{'allele': important_alleles[0], 'snps': filtered_snps}]
      
    all_alleles = important_alleles + unimportant_alleles
    num_important = len(important_alleles)
    num_unimportant = len(unimportant_alleles)
    total = len(all_alleles)
    snp_map = defaultdict(lambda: defaultdict(list))
    for allele in all_alleles:
      seq = self.alleles[allele]
      for item in seq:
        snp_map[item]['present'].append(allele)
        
    all_snp_scores = []
    for snp, allele_map in snp_map.iteritems():
      important_snp_allele = [allele for allele in allele_map['present'] if allele in important_alleles]
      raw_score = len(important_snp_allele) / float(num_important)
      if '_' not in snp and '|' not in snp:
        all_snp_scores.append({'snp': snp, 'score': raw_score})
        
    sorted_scores = sorted(all_snp_scores, key=itemgetter('score')) 
    for to_test in sorted_scores:
      if to_test['snp'] not in filtered_snps:
        snp_filter.append(to_test)
        allele_test = self.check_alleles2(important, unimportant, snp_filter)
        snp_filter.pop()
        if allele_test:
          return allele_test + self.check_alleles2(remaining_important, unimportant, [])
    return []
  
  def does_not_contain_alleles(self, allele_combo, test_alleles):
    return not all(allele in test_alleles for allele in allele_combo)
    
  def test_allele_combinations(self, important_allele, seq):
    for i in range(len(seq)):
      for allele_combo in combinations(self.important_alleles[important_allele], i):
        if all(self.does_not_contain_alleles(allele_combo, v) for v in self.unimportant_alleles.values()):
          return allele_combo
  
  def parse_alleles(self, allele_file, wildtype):
    allele_definitions = {}
    with open(allele_file) as f:
      for i, line in enumerate(f):
        if i < 3:
          continue
        allele_def = line.split()
        allele_string = allele_def[1:]
        # if allele_def[0] == wildtype:
        #   allele_string = "".join(allele_string).translate(string.maketrans('ACGT', '____'))
        #   print allele_string
        allele_definitions[allele_def[0]] = ["%s:%s" % (i, base) for i, base in enumerate(allele_string)]
    return allele_definitions
    
  def snp_locs(self, allele1, allele2):
    return ["%s:%s" % (i, allele2[i]) for i, base in enumerate(allele1) if base != allele2[i] ]
    
  
if __name__ == "__main__":
  
  allele_diff = AlleleDiff()
  important_alleles = ["B*51240301"]
  #important_alleles = ["B*51240301", "B*55070101", "B*35700101", "B*07420101", "B*27090101"]
  unimportant_alleles = [allele for allele in allele_diff.alleles.keys() if allele not in important_alleles]
  allele_mapping = allele_diff.check_alleles(important_alleles, unimportant_alleles)
  print allele_mapping
  #print allele_diff.unimportant_alleles
