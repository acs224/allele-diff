import string
import operator
import logging
import itertools
import numpy
from collections import defaultdict, Counter
from operator import itemgetter
from optparse import OptionParser

class HaploBuilder:
  """Object to facilitate the building of a haplotype translation file."""
    
  def __init__(self, haplo_file="HLAB.txt"):
    """Construct an HaploBuilder object.
    
    Kwargs: 
      haplo_file: File to use for haplotype definitions.
      
    """
    # Number of header rows in the haplotype files
    self.HEADER_ROWS = 3
    # Special characters in haplotype file
    self.EXON = '|'
    self.NO_CALL = '*'
    self.CALLED_BASES = ['A', 'C', 'G', 'T', '_', '.']
    
    self.haplos = self.parse_haplos(haplo_file)
    self.call_cache = defaultdict(lambda: defaultdict(list))

  def parse_haplos(self, haplo_file):
    """Parse the haplo_file into a list of POS:BASE combinations.
    
    Args:
      haplo_file: File to parse
      
    Returns a dict where haplotype is the key and the value is the list of SNP definitions for that haplotype.
      
    """
    haplo_definitions = {}
    allele_list = []
    self.allele_index = []
    self.allele_hash = {}
    self.unidentifiable = []
    with open(haplo_file) as f:
      for i, line in enumerate(f):
        if i < self.HEADER_ROWS:
          continue
        allele_def = line.split()
        allele_string = allele_def[1:]
        allele_list.append(allele_string)
        string_hash = "".join(allele_string)
        if string_hash in self.allele_hash:
          #print "ERROR: Allele %s is duplicate of %s.  Cannot uniquely identify either allele." % (allele_def[0], self.allele_hash[string_hash])
          self.unidentifiable = self.unidentifiable + self.allele_hash[string_hash]
          self.unidentifiable.append(allele_def[0])
        else:
          self.allele_hash[string_hash] = []
        self.allele_hash[string_hash].append(allele_def[0])
        self.allele_index.append(allele_def[0])
    self.data_matrix = numpy.array(allele_list)
    return haplo_definitions
 
  def check_haplos(self, important, unimportant, snp_filter):
    """Recursively builds a haplotype translation table based upon the important and unimportant haplotypes.
    
    Args:
      important: List of important haplotypes.
      unimportant: List of unimportant haplotypes.
      snp_filter: List of snps to use in filtering haplotypes.
      
    Returns a list of dictonaries, one for each important haplotype.  Each dictonary is formatted:
    {haplotype: haplotype name, snps: snps to translate that haplotype}
    
    """
    data = self.check_haplos2(important, unimportant, snp_filter, {})
    print data
    return data
    
  def check_haplos2(self, important, unimportant, snp_filter, defining_snps):
    
    identifiable = [hap for hap in important if hap not in self.unidentifiable]
    if not identifiable:
      return []
    scored_snps = self.score_snps2(identifiable, unimportant, snp_filter)
    
    brand_new_snp_filter = list(snp_filter)
    added_snps = False
    for snp in scored_snps:
      new_snp_filter = list(snp_filter)
      new_snp_filter.append(snp['snp'])
      new = unimportant + defining_snps.keys()
      call = self.is_callable(identifiable, new, new_snp_filter)
      if call:
        print "callable:", call, new_snp_filter
        for hap, snps in call.iteritems():
          if hap not in defining_snps:
            defining_snps[hap] = snps
        print "DEFINING", defining_snps
        brand_new_snp_filter.append(snp['snp'])
        added_snps = True
          
    uncallable_haps = [hap for hap in identifiable if hap not in defining_snps.keys()]
    if not uncallable_haps:
      return defining_snps
    else:
      new = unimportant + defining_snps.keys()
      scored_snps = self.score_snps2(uncallable_haps, new, snp_filter)
      print "fitler", brand_new_snp_filter, uncallable_haps, scored_snps[0]
      if not added_snps:
        brand_new_snp_filter.append(scored_snps[0]['snp'])
      return self.check_haplos2(uncallable_haps, unimportant, brand_new_snp_filter, defining_snps)
    
  def filter_haps(self, hap_list, snps):
    alleles = []
    locs = []
    for snp in snps:
      loc, allele = snp.split(":")
      alleles.append(allele)
      locs.append(int(loc))
    result = self.data_matrix[:, locs] == alleles
    result_any = numpy.any(result, axis=1)
    result_all = numpy.all(result, axis=1)
    index_any = numpy.where(result_any == True)
    index_all = numpy.where(result_all == True)
    
    good_haps = []
    for index in index_all[0]:
      hap_name = self.allele_index[index]
      if hap_name in hap_list:
        good_haps.append(hap_name)
    return good_haps
    
  def call_haplotype(self, haplotype, unimportant, snp_filter, defining_snps):
    already_called = [snp for snp_list in defining_snps.values() for snp in snp_list]
    scored_snps = self.score_snps3([haplotype], unimportant, snp_filter, already_called)
    brand_new_snp_filter = list(snp_filter)
    added_snps = False
    for snp in scored_snps:
      print snp
      new_snp_filter = list(snp_filter)
      new_snp_filter.append(snp['snp'])
      call = self.is_callable2([haplotype], unimportant, new_snp_filter)
      if call:
        print "about to return:", call, new_snp_filter
        return call
      else:
        results = self.call_haplotype(haplotype, unimportant, new_snp_filter, defining_snps)
        if results:
          return results
        else:
          return []
          
    return []
  
  def call_haplotypes(self, important, unimportant, snp_filter):
    identifiable = [hap for hap in important if hap not in self.unidentifiable]
    if not identifiable:
      return []
    all_haps = {}
    for haplotype in identifiable:
      print "TRYINT TO CALL", haplotype
      new = unimportant + [hap for hap in identifiable if hap != haplotype]
      result = self.call_haplotype(haplotype, new, snp_filter, all_haps)
      print "CALLED:", result
      all_haps.update(result)
    return all_haps

  def is_callable(self, important, unimportant, snps):
    matching_alleles = self.find_matching(important, unimportant, snps)
    snp_map = defaultdict(list)
    good_snp_map = {}
    #print snps, matching_alleles
    for haplotype, matching_snps in matching_alleles.iteritems():
      for i in range(1, len(matching_snps) + 1):
        for combo in itertools.combinations(matching_snps, i):
          snp_hash = "|".join(combo)
          snp_map[snp_hash].append(haplotype)
    #print snp_map
    for snp_hash, haps in snp_map.iteritems():
      if len(haps) == 1 and haps[0] in important:
        good_snp_map[haps[0]] = snp_hash.split("|")
    #print "MAP:",  snp_map, good_snp_map
    #print snps, len(matching_alleles), matching_alleles
    if good_snp_map:
      #print "MATCHING", matching_alleles
      #print "snp_map", snp_map
      #print "returning",  good_snp_map
      return good_snp_map
    return False
    
  def is_callable2(self, important, unimportant, snps):
    matching_alleles = self.find_matching2(important, unimportant, snps)
    snp_map = defaultdict(list)
    good_snp_map = {}
    print important, snps, matching_alleles
    
    for haplotype, matching_snps in matching_alleles.iteritems():
      snp_hash = "|".join(matching_snps)
      snp_map[snp_hash].append(haplotype)
    #print important, snp_map
    for snp_hash, haps in snp_map.iteritems():
      if len(haps) == 1 and haps[0] in important:
        good_snp_map[haps[0]] = snp_hash.split("|")
    print "MAP:",  snp_map, good_snp_map
    #print snps, len(matching_alleles), matching_alleles
    if good_snp_map:
      #print "MATCHING", matching_alleles
      #print "snp_map", snp_map
      print "returning",  good_snp_map
      return good_snp_map
    return False
    
  def print_matrix(self,  sample_index, matrix):
    for i, row in enumerate(matrix):
      print self.allele_index[sample_index[i]], row
  
  def find_matching(self, important, unimportant, snps):
    
    alleles = []
    locs = []
    matching_alleles = defaultdict(list)
    for snp in snps:
      loc, allele = snp.split(":")
      alleles.append(allele)
      locs.append(int(loc))
    result = self.data_matrix[:, locs] == alleles
    result_any = numpy.any(result, axis=1)
    result_all = numpy.all(result, axis=1)
    index_any = numpy.where(result_any == True)
    index_all = numpy.where(result_all == True)
    
    for index in index_any[0]:
      hap_name = self.allele_index[index]
      compare_results = result[index]
      for i, r in enumerate(compare_results):
        if r:
          matching_alleles[hap_name].append(snps[i])
    return matching_alleles
    
  def find_matching2(self, important, unimportant, snps):

    alleles = []
    locs = []
    matching_alleles = defaultdict(list)
    for snp in snps:
      loc, allele = snp.split(":")
      alleles.append(allele)
      locs.append(int(loc))
    result = self.data_matrix[:, locs] == alleles
    result_any = numpy.any(result, axis=1)
    result_all = numpy.all(result, axis=1)
    index_any = numpy.where(result_any == True)
    index_all = numpy.where(result_all == True)

    for index in index_all[0]:
      hap_name = self.allele_index[index]
      compare_results = result[index]
      for i, r in enumerate(compare_results):
        if r:
          matching_alleles[hap_name].append(snps[i])
    return matching_alleles
    
  def score_snps2(self, important_haplos, unimportant_haplos, snp_filter):
    """Assign a score to each snp.  This score is meant to minimize the number
    of snps that should be used to define a haplotype.  The current formula is:
    
    score = (ratio of snp in important haplos + .5) / (ratio of snp in unimportant haplos + .5) + (weighted score to favor snps over reference)
    Note: .5 is added to each ratio in case the denominator is 0
    
    Args:
      important_haplos: List of important haplotypes.
      unimportant_haplos: List of unimportant haplotypes.
      snp_map: Dictionary of all snps.
      
    Returns a list of scored snps. (Higher score is better, has more chance of having useful information)
    """
    all_snp_scores = []
    unimportant_filtered = unimportant_haplos
    important_filtered = important_haplos
    if snp_filter:
      print "FILTER", snp_filter
      unimportant_filtered = self.filter_haps(unimportant_haplos, snp_filter)
      #important_filtered = self.filter_haps(important_haplos, snp_filter)
    important_indexes = [i for i, item in enumerate(self.allele_index) if item in important_filtered]
    unimportant_indexes = [i for i, item in enumerate(self.allele_index) if item in unimportant_filtered]
    
    for i, col in enumerate(self.data_matrix.T):
      important_counter = Counter(col[important_indexes,:])
      unimportant_counter = Counter(col[unimportant_indexes,:])
      for base in self.CALLED_BASES:
        if base not in important_counter and base not in unimportant_counter:
          continue
        snp = "%s:%s" % (i, base)
        if snp in snp_filter:
          continue
        important_ratio = 0
        if len(important_indexes) > 0:
          important_ratio = important_counter[base] / float(len(important_indexes))
        unimportant_ratio = 0
        if len(unimportant_indexes) > 0:
          unimportant_ratio = unimportant_counter[base] / float(len(unimportant_filtered))
        if important_ratio > 0:
          raw_score = (1 - important_ratio) + (1 - unimportant_ratio) * (.9 if "_" in snp else 1)
          if raw_score > 0:
            all_snp_scores.append({'snp': snp, 'score': raw_score, 'import_ratio': important_ratio, 'unimport_ratio': unimportant_ratio})
    
    return sorted(all_snp_scores, key=itemgetter('score'), reverse=True)
    
  def score_snps3(self, important_haplos, unimportant_haplos, snp_filter, used_snps):
    """Assign a score to each snp.  This score is meant to minimize the number
    of snps that should be used to define a haplotype.  The current formula is:

    score = (ratio of snp in important haplos + .5) / (ratio of snp in unimportant haplos + .5) + (weighted score to favor snps over reference)
    Note: .5 is added to each ratio in case the denominator is 0

    Args:
      important_haplos: List of important haplotypes.
      unimportant_haplos: List of unimportant haplotypes.
      snp_map: Dictionary of all snps.

    Returns a list of scored snps. (Higher score is better, has more chance of having useful information)
    """
    all_snp_scores = []
    unimportant_filtered = unimportant_haplos
    important_filtered = important_haplos
    if snp_filter:
      print "FILTER", snp_filter
      unimportant_filtered = self.filter_haps(unimportant_haplos, snp_filter)
      #important_filtered = self.filter_haps(important_haplos, snp_filter)
    important_indexes = [i for i, item in enumerate(self.allele_index) if item in important_filtered]
    unimportant_indexes = [i for i, item in enumerate(self.allele_index) if item in unimportant_filtered]

    for i, col in enumerate(self.data_matrix.T):
      important_counter = Counter(col[important_indexes,:])
      unimportant_counter = Counter(col[unimportant_indexes,:])
      for base in self.CALLED_BASES:
        if base not in important_counter and base not in unimportant_counter:
          continue
        snp = "%s:%s" % (i, base)
        if snp in snp_filter:
          continue
        # important_ratio = 0
        # if len(important_indexes) > 0:
        important_ratio = important_counter[base] / float(len(important_indexes))
        # unimportant_ratio = 0
        # if len(unimportant_indexes) > 0:
        unimportant_ratio = unimportant_counter[base] / float(len(unimportant_filtered))
        if important_ratio > 0:
          raw_score = (1 - important_ratio) + (1 - unimportant_ratio) * (.9 if "_" in snp else 1) * (.9 if snp in used_snps else 1)
          if raw_score > 0:
            all_snp_scores.append({'snp': snp, 'score': raw_score, 'import_ratio': important_ratio, 'unimport_ratio': unimportant_ratio})

    return sorted(all_snp_scores, key=itemgetter('score'), reverse=True)
  
    
  
if __name__ == "__main__":
  parser = OptionParser()
  parser.add_option("-f", "--file", dest="filename", help="File with good haplotypes to read.  One haplotype per line. (Required)")
  parser.add_option("-t", "--types", dest="hap_file", help="File with all haplotypes to read.")
  parser.add_option("-d", "--debug", action="store_false", dest="debug", default=False, help="Print debug messages.")
  
  (options, args) = parser.parse_args()
  if not options.filename: 
    parser.error('--file is required')
  if options.debug:
    logging.basicConfig(level=logging.DEBUG)
  haplo_builder = HaploBuilder() if not options.hap_file else HaploBuilder(options.hap_file)
  important_haplos = open(options.filename).read().split("\n")
  unimportant_haplos = [hap for hap in haplo_builder.allele_index if hap not in important_haplos]
  haplo_mapping = haplo_builder.call_haplotypes(important_haplos, unimportant_haplos, [])
  for item, snps in haplo_mapping.iteritems():
    print "%s\t%s" % (item, "\t".join(sorted(snps)))
  # 
  # important_haplos = ['B*55070101']
  # all_haplo_diff = HaploBuilder("HLAB.txt")
  # unimportant_haplos = [hap for hap in all_haplo_diff.haplos.keys() if hap not in important_haplos]
  # checked = all_haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
  # print checked
