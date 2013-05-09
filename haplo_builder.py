import string
import operator
import logging
import itertools
import numpy
from collections import defaultdict
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
        # Ignore positions that are NO_CALL or EXON definitions.  
        # NO_CALL can be ignored since it does not contain any basepair definition.
        alleles = ["%s:%s" % (i, base) for i, base in enumerate(allele_string) if (self.EXON != base and self.NO_CALL != base)]
        haplo_definitions[allele_def[0]] = alleles
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
    # logging.debug("In check_haplos.  Important: %s Filter: %s" % (important, snp_filter))
    # important_haplos = self.filter_haplos(important, snp_filter)
    # unimportant_haplos = self.filter_haplos(unimportant, snp_filter, True)
    # remaining_important = [haplo for haplo in important if haplo not in important_haplos]
    # filtered_snps = [temp_snp['snp'] for temp_snp in snp_filter]
    # if not important or len(important_haplos) == 0:
    #   return []
    # logging.debug("Before recursion check.  snp_filter: %s important %s unimportant %s" % (snp_filter, important_haplos, unimportant_haplos))
    # if len(snp_filter) > 0 and len(important_haplos) == 1 and len(unimportant_haplos) == 0:
    #   logging.debug("About to recurse.  important: %s remaining_important %s" % (important_haplos[0], remaining_important))
    #   return [{'haplotype': important_haplos[0], 'snps': filtered_snps}] + self.check_haplos(remaining_important, unimportant, [])
    #   
    # all_haplos = important_haplos + unimportant_haplos
    # snp_map = defaultdict(lambda: defaultdict(list))
    # for haplo in all_haplos:
    #   seq = self.haplos[haplo]
    #   for item in seq:
    #     snp_map[item]['present'].append(haplo)
    #     
    # all_snp_scores = self.score_snps(important_haplos, unimportant_haplos, snp_map)
    # for to_test in all_snp_scores:
    #   if to_test['snp'] not in filtered_snps:
    #     snp_filter.append(to_test)
    #     return self.check_haplos(important, unimportant, snp_filter)
    # return []
    data = self.check_haplos2(important, unimportant, snp_filter, {})
    print data
    return data
    
  def check_haplos2(self, important, unimportant, snp_filter, defining_snps):
    snp_map = defaultdict(lambda: defaultdict(list))
    all_haps = (important + unimportant)
    for haplo in all_haps:
      seq = self.haplos[haplo]
      #seq = self.get_haplo(haplo)
      for item in seq:
        if item not in snp_filter:
          snp_map[item]['present'].append(haplo)
          
    scored_snps = self.score_snps(important, unimportant, snp_map)
    
    results = []
    
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
        #print "CALL:", hap, new_snp_filter, call
        if call:
          print "callable:", hap, new_snp_filter
          defining_snps[hap] = call
          brand_new_snp_filter.append(snp['snp'])
          
    uncallable_haps = [hap for hap in important if hap not in defining_snps.keys()]
    #print uncallable_haps
    if not uncallable_haps:
      return defining_snps
    else:
      scored_snps = self.score_snps(uncallable_haps, unimportant, snp_map)
      print "fitler", brand_new_snp_filter
      brand_new_snp_filter.append(scored_snps[0]['snp'])
      return self.check_haplos2(important, unimportant, brand_new_snp_filter, defining_snps)
    
  def is_callable(self, haplotype, others, snps):
    # hap_seq =  self.haplos[haplotype]
    # for i in range(1, len(snps) + 1):
    #   for combo in itertools.combinations(snps, i):
    #     if haplotype in self.call_cache and combo in self.call_cache[haplotype]:
    #       if not self.call_cache[haplotype][combo]:
    #         continue
    #       else:
    #         #print "returning cached:", haplotype, combo, self.call_cache[haplotype][combo]
    #         return self.call_cache[haplotype][combo]
    #     all_other_calls = []
    #     can_call = all([snp in hap_seq for snp in combo])
    #     #print haplotype, combo, can_call, [snp in hap_seq for snp in combo]
    #     for haplo in others:
    #       seq = self.haplos[haplo]
    #       other_call = all([snp in seq for snp in combo])
    #       all_other_calls.append(other_call)
    #     if can_call and not any(all_other_calls):
    #       self.call_cache[haplotype][combo] = list(combo)
    #       return list(combo)
    #     else:
    #       self.call_cache[haplotype][combo] = False
    # #print self.call_cache
    # return False
    for i in range(1, len(snps) + 1):
      for combo in itertools.combinations(snps, i):
        matching_alleles = self.find_matching(combo)
        #print "matches", haplotype, matching_alleles, (len(matching_alleles) == 1 and haplotype in matching_alleles), combo
        if len(matching_alleles) == 1 and haplotype in matching_alleles:
          return list(combo)
    return False
  
  def filter_haplos(self, haplos, snp_filter, include=True):
    """Filters haplotypes out based on if they have certain snps or not.
    
    Args:
      haplos: List of haplotypes to filter
      snp_filter: List of snps to use as filter criteria
      
    Kwargs:
      include:  True if criteria should be used to include haplotypes, False if it should
                exclude haplotypes
      
    Returns a list of haplotypes that have been filtered based on criteria
    
    """
    if len(snp_filter) < 1:
      return haplos
    to_filter = snp_filter
    filtered_haplos = []
    if not include:
      filtered_haplos = list(haplos)
    for haplo in haplos:
      seq = self.haplos[haplo]
      if include and all([snp in seq for snp in to_filter]):
        filtered_haplos.append(haplo)
      elif not include and not any([snp in seq for snp in to_filter]):
        filtered_haplos.remove(haplo)

    return filtered_haplos
    
  def extract(self, snps, alleles):
    indexes = [i for i, item in enumerate(self.allele_index) if item in alleles]
    return self.data_matrix[:, snps][indexes, :]
    
  def get_haplo(self, haplotype):
    indexes = self.allele_index.index(haplotype)
    return self.data_matrix[indexes, :]
    
  def find_matching(self, snps):
    snp_hash = "".join(snps)
    if snp_hash in self.call_cache:
      return self.call_cache[snp_hash]
    matching = []
    alleles = []
    locs = []
    for snp in snps:
      loc, allele = snp.split(":")
      alleles.append(allele)
      locs.append(int(loc))
    #print locs
    for i, row in enumerate(self.data_matrix[:, locs]):
      if all(row == alleles):
        matching.append(self.allele_index[i])
      #print row, alleles, all(row == alleles)
    self.call_cache[snp_hash] = matching
    return matching
    
  def score_snps(self, important_haplos, unimportant_haplos, snp_map):
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
    for snp, haplo_map in snp_map.iteritems():
      
      important_snp_haplos = [hap for hap in haplo_map['present'] if hap in important_haplos]
      unimportant_snp_haplos = [hap for hap in  haplo_map['present'] if hap in unimportant_haplos]
      important_ratio = len(important_snp_haplos) / float(len(important_haplos))
      unimportant_ratio = 0
      if len(unimportant_haplos) > 0:
        unimportant_ratio = len(unimportant_snp_haplos) / float(len(unimportant_haplos))
      if important_ratio > 0:
        raw_score = ((important_ratio + .5) / (unimportant_ratio + .5)) + (0 if "_" in snp else 1)
      #raw_score = len(important_snp_haplos)
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
