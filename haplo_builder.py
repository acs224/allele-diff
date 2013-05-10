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
    self.REF = '_'
    
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
    with open(haplo_file) as f:
      for i, line in enumerate(f):
        if i < self.HEADER_ROWS:
          continue
        allele_def = line.split()
        allele_string = allele_def[1:]
        allele_list.append(allele_string)
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
    scored_snps = self.score_snps2(important, unimportant, snp_filter)
    
    uncallable_haps = []
    brand_new_snp_filter = list(snp_filter)
    for snp in scored_snps:
      new_snp_filter = list(snp_filter)
      new_snp_filter.append(snp['snp'])
      for hap in important:
        if hap in defining_snps:
          continue
        others = [i for i in important if i != hap and hap not in defining_snps]
        others = others + unimportant
        call = self.is_callable(hap, others, new_snp_filter)
        if call:
          print "callable:", hap, new_snp_filter
          defining_snps[hap] = call
          brand_new_snp_filter.append(snp['snp'])
          
    uncallable_haps = [hap for hap in important if hap not in defining_snps.keys()]
    if not uncallable_haps:
      return defining_snps
    else:
      scored_snps = self.score_snps2(uncallable_haps, unimportant, snp_filter)
      print "fitler", brand_new_snp_filter
      brand_new_snp_filter.append(scored_snps[0]['snp'])
      return self.check_haplos2(important, unimportant, brand_new_snp_filter, defining_snps)
    
  def is_callable(self, haplotype, others, snps):
    matching_alleles = self.find_matching(snps)
    haplotype_alleles = []
    other_alleles = []
    for hap, alleles in matching_alleles.iteritems():
      all_combos = []
      for i in range(1, len(snps) + 1):
        for combo in itertools.combinations(alleles, i):
          all_combos.append(combo)
      if hap == haplotype:
        for combo in all_combos:
          haplotype_alleles.append(combo)
      else:
        for combo in all_combos:
          other_alleles.append(combo)
    not_in_other = [combo for combo in haplotype_alleles if combo not in other_alleles]
    #print haplotype_alleles, other_alleles, not_in_other
    if not_in_other:
      return list(not_in_other)
    return False
  
  def find_matching(self, snps):
    # snp_hash = "".join(snps)
    # if snp_hash in self.call_cache:
    #   return self.call_cache[snp_hash]
    # alleles = []
    # locs = []
    # for snp in snps:
    #   loc, allele = snp.split(":")
    #   alleles.append(allele)
    #   locs.append(int(loc))
    # result = self.data_matrix[:, locs] == alleles
    # result_all = numpy.all(result, axis=1)
    # matching_index = numpy.where(result_all == True)
    # matching = [hap for i, hap in enumerate(self.allele_index) if i in matching_index[0]]
    # self.call_cache[snp_hash] = matching
    # return matching
    
    alleles = []
    locs = []
    matching_alleles = defaultdict(list)
    for snp in snps:
      loc, allele = snp.split(":")
      alleles.append(allele)
      locs.append(int(loc))
    result = self.data_matrix[:, locs] == alleles
    result_any = numpy.any(result, axis=1)
    matching_index = numpy.where(result_any == True)
    
    #print "RESULT",  result, snps
    for index in matching_index[0]:
      hap_name = self.allele_index[index]
      compare_results = result[index]
      #print "COMP", compare_results, snps, hap_name
      for i, r in enumerate(compare_results):
        if r:
          matching_alleles[hap_name].append(snps[i])
    #print matching_alleles
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
    important_indexes = [i for i, item in enumerate(self.allele_index) if item in important_haplos]
    unimportant_indexes = [i for i, item in enumerate(self.allele_index) if item in unimportant_haplos]
    for i, col in enumerate(self.data_matrix.T):
      important_counter = Counter(col[important_indexes,:])
      unimportant_counter = Counter(col[unimportant_indexes,:])
      for base in ['A', 'C', 'G', 'T', '_']:
        if base not in important_counter and base not in unimportant_counter:
          continue
        snp = "%s:%s" % (i, base)
        if snp in snp_filter:
          continue
        important_ratio = (important_counter[base] + .5) / (len(important_indexes) + .5)
        unimportant_ratio = (unimportant_counter[base] + .5) / (len(unimportant_indexes) + .5)
        if important_ratio > 0:
          raw_score = ((important_ratio + .5) / (unimportant_ratio + .5)) + (0 if base == '_' else 1)
          all_snp_scores.append({'snp': snp, 'score': raw_score})
    
    return sorted(all_snp_scores, key=itemgetter('score'), reverse=True)
    
  
if __name__ == "__main__":
  # parser = OptionParser()
  # parser.add_option("-f", "--file", dest="filename", help="File with good haplotypes to read.  One haplotype per line. (Required)")
  # parser.add_option("-t", "--types", dest="hap_file", help="File with all haplotypes to read.")
  # parser.add_option("-d", "--debug", action="store_false", dest="debug", default=False, help="Print debug messages.")
  # 
  # (options, args) = parser.parse_args()
  # if not options.filename: 
  #   parser.error('--file is required')
  # if options.debug:
  #   logging.basicConfig(level=logging.DEBUG)
  # haplo_builder = HaploBuilder() if not options.hap_file else HaploBuilder(options.hap_file)
  # important_haplos = open(options.filename).read().split("\n")
  # unimportant_haplos = [hap for hap in haplo_builder.haplos.keys() if hap not in important_haplos]
  # haplo_mapping = haplo_builder.check_haplos(important_haplos, unimportant_haplos, [])
  # for item in haplo_mapping:
  #   print "%s\t%s" % (item['haplotype'], "\t".join(sorted(item['snps'])))
  
  important_haplos = ['B*55070101']
  all_haplo_diff = HaploBuilder("HLAB.txt")
  unimportant_haplos = [hap for hap in all_haplo_diff.haplos.keys() if hap not in important_haplos]
  checked = all_haplo_diff.check_haplos(important_haplos, unimportant_haplos, [])
  print checked
