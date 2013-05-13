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
    
    self.parse_haplos(haplo_file)

  def parse_haplos(self, haplo_file):
    """Parse the haplo_file into a list of POS:BASE combinations.
    
    Args:
      haplo_file: File to parse
      
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
          self.unidentifiable = self.unidentifiable + self.allele_hash[string_hash]
          self.unidentifiable.append(allele_def[0])
        else:
          self.allele_hash[string_hash] = []
        self.allele_hash[string_hash].append(allele_def[0])
        self.allele_index.append(allele_def[0])
    self.data_matrix = numpy.array(allele_list)
  
  def filter_matrix(self, snps):
    alleles = []
    locs = []
    matching_alleles = defaultdict(list)
    for snp in snps:
      loc, allele = snp.split(":")
      alleles.append(allele)
      locs.append(int(loc))
    return self.data_matrix[:, locs] == alleles
  
  def find_matching(self, snps):

    result = self.filter_matrix(snps)
    result_all = numpy.all(result, axis=1)
    index_all = numpy.where(result_all == True)

    matching_alleles = defaultdict(list)
    for index in index_all[0]:
      hap_name = self.allele_index[index]
      compare_results = result[index]
      for i, r in enumerate(compare_results):
        if r:
          matching_alleles[hap_name].append(snps[i])
    return matching_alleles
  
  def filter_haps(self, hap_list, snps):

    result = self.filter_matrix(snps)
    result_any = numpy.any(result, axis=1)
    index_any = numpy.where(result_any == True)
    
    good_haps = []
    for index in index_any[0]:
      hap_name = self.allele_index[index]
      if hap_name in hap_list:
        good_haps.append(hap_name)
    return good_haps
    
  def call_haplotype(self, haplotype, unimportant, snp_filter, defining_snps):
    already_called = [snp for snp_list in defining_snps.values() for snp in snp_list]
    scored_snps = self.score_snps([haplotype], unimportant, snp_filter, already_called)
    for snp in scored_snps:
      new_snp_filter = list(snp_filter)
      new_snp_filter.append(snp['snp'])
      call = self.make_call(haplotype, unimportant, new_snp_filter)
      if call:
        return call
      else:
        results = self.call_haplotype(haplotype, unimportant, new_snp_filter, defining_snps)
        if results:
          return results
    return []
  
  def call_haplotypes(self, important, unimportant, snp_filter):
    """Recursively builds a haplotype translation table based upon the important and unimportant haplotypes.
    
    Args:
      important: List of important haplotypes.
      unimportant: List of unimportant haplotypes.
      snp_filter: List of snps to use in filtering haplotypes.
      
    Returns a list of dictonaries, one for each important haplotype.  Each dictonary is formatted:
    {haplotype: haplotype name, snps: snps to translate that haplotype}
    
    """
    
    identifiable = [hap for hap in important if hap not in self.unidentifiable]
    unidentifiable = [hap for hap in important if hap not in identifiable]
    if not identifiable:
      return []
    all_haps = {}
    for haplotype in identifiable:
      others = unimportant + [hap for hap in identifiable if hap != haplotype]
      result = self.call_haplotype(haplotype, others, snp_filter, all_haps)
      all_haps.update(result)
      
    for hap in unidentifiable:
      all_haps[hap] = ['Unidentifiable']
    return all_haps

  def make_call(self, haplotype, unimportant, snps):
    matching_alleles = self.find_matching(snps)
    snp_map = defaultdict(list)
    good_snp_map = {}
    for haplotype, matching_snps in matching_alleles.iteritems():
      snp_hash = "|".join(matching_snps)
      snp_map[snp_hash].append(haplotype)
    for snp_hash, haps in snp_map.iteritems():
      if len(haps) == 1 and haps[0] == haplotype:
        good_snp_map[haps[0]] = snp_hash.split("|")
    if good_snp_map:
      return good_snp_map
    return False
    
  def score_snps(self, important_haplos, unimportant_haplos, snp_filter, used_snps):
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
      unimportant_filtered = self.filter_haps(unimportant_haplos, snp_filter)
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
        important_ratio = important_counter[base] / float(len(important_indexes))
        unimportant_ratio = unimportant_counter[base] / float(len(unimportant_filtered))
        if important_ratio > 0:
          raw_score = (1 - important_ratio) + (1 - unimportant_ratio) * (.9 if "_" in snp else 1)
          if raw_score > 0:
            all_snp_scores.append({'snp': snp, 'score': raw_score})

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
