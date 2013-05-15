import string
import operator
import logging
import itertools
import numpy
from collections import defaultdict, Counter
from operator import itemgetter
from optparse import OptionParser

class HaploBuilder:
  """
  Object to facilitate the building of a haplotype translation file.  The HaploBuilder attempts to find 
  the minimal amount of SNPs needed to call a set of good haplotypes.  It scores a set of snps, then uses the
  highest scoring snp to attempt to call all snps via recursion.  The complexity of this algorithm in the worst 
  case will be O((nk)^2) where n is the number of haplotypes and k is the length of the longest haplotype.  For each
  snp in k, it must be added to the filter, scored and then tested to see if it completes a haplotype.
  
  """
    
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
    self.BASES = ['A', 'C', 'G', 'T', '_', '.']
    self.SNP_BASES = ['A', 'C', 'G', 'T', '.']
    
    self.parse_haplos(haplo_file)

  def parse_haplos(self, haplo_file):
    """
    Parse the haplo_file into a NumPy 2D array where the columns are snps and the rows are haplotypes.  
    Also looks for duplicate haplotypes and holds them in the synonyms dictionary.
    
    Args:
      haplo_file: File to parse
      
    """
    self.allele_index = []
    self.synonyms = {}
    
    haplo_definitions = {}
    allele_list = []
    allele_hash = {}
    with open(haplo_file) as f:
      for i, line in enumerate(f):
        if i < self.HEADER_ROWS:
          continue
        allele_def = line.split()
        hap_name = allele_def[0]
        allele_string = allele_def[1:]
        string_hash = "".join(allele_string)
        if string_hash in allele_hash:
          self.synonyms[hap_name] = allele_hash[string_hash][0]
          if allele_hash[string_hash][0] not in self.synonyms:
            self.synonyms[allele_hash[string_hash][0]] = allele_hash[string_hash][0]
        else:
          allele_hash[string_hash] = []
          allele_list.append(allele_string)
          self.allele_index.append(hap_name)
        allele_hash[string_hash].append(hap_name)
    
    self.data_matrix = numpy.array(allele_list)
  
  def filter_matrix(self, locs):
    """
    Filters the NumPy matrix by snp locations.

    Args:
      locs: List of snp locations.  i.e. [2, 5, 7]

    Returns a submatrix of just the provided locations.
    """
    
    return self.data_matrix[:, locs]
  
  def find_matching(self, snps):
    """
    Finds all the alleles for a set of snps for each haplotype.
    
    Args:
      snps: List of snp locations.  i.e. [2, 3]

    Returns a dictionary where the haplotype name is the key and the values are the list
    of snps that haplotype contains.
    """
    
    result = self.filter_matrix(snps)
    matching_alleles = defaultdict(list)
    for index, hap_name in enumerate(self.allele_index):
      compare_results = result[index]
      for i, r in enumerate(compare_results):
        allele_snps = []
        if r:
          for j in range(0, len(result[index])):
            allele_snps.append("%s:%s" % (snps[j], result[index][j]))
        matching_alleles[hap_name] = allele_snps
    return matching_alleles
    
  def match(self, snps):
    """
    Attempts to find a match for the provided set of snps and alleles.

    Args:
      snps: List of snp and alleles.  i.e. [2:C, 4:G]

    Returns a of all haplotypes that match those snp/allele combinations.
    """
    
    alleles = []
    locs = []
    matching = []
    for snp in snps:
      loc, allele = snp.split(":")
      alleles.append(allele)
      locs.append(int(loc))
    result = self.data_matrix[:, locs] == alleles
    result_all = numpy.all(result, axis=1)
    index_all = numpy.where(result_all == True)

    for index in index_all[0]:
      hap_name = self.allele_index[index]
      compare_results = result[index]
      for i, r in enumerate(compare_results):
        if r:
          matching.append(hap_name)
    return list(set(matching))
    
  
  def filter_haps(self, hap_list, snps):
    """
    Returns all the haplotypes that have a snp at any of the provided locations.

    Args:
      hap_list: List of haplotypes.
      snps: List of snps locations to filter by.

    Returns a list of haplotypes that have at least one of the snps.
    """
    
    result = self.filter_matrix(snps)
    result_any = numpy.any(result != self.REF, axis=1)
    index_any = numpy.where(result_any == True)
    
    good_haps = []
    for index in index_any[0]:
      hap_name = self.allele_index[index]
      if hap_name in hap_list:
        good_haps.append(hap_name)
    return good_haps

  def check_haplos(self, important, unimportant):
    """
    Wraps the recursive method, allows for collapsing of returned haplotypes to get the minimal set.
    Any duplicate haplotypes will be marked as 'Unidentifiable' since they cannot be uniquely called.

    Args:
      important: List of important haplotypes.
      unimportant: List of unimportant haplotypes.

    Returns a dictionary of haplotypes and the snps used to call them.
    """
    identifiable = [hap for hap in important if hap not in self.synonyms]
    unidentifiable = [hap for hap in important if hap not in identifiable]
    data = self.check_all(identifiable, unimportant, [], {})
    collapsed = self.collapse_haps(data)
    for hap in unidentifiable:
      collapsed[hap] = 'Unidentifiable'
    return collapsed
    
  def check_all(self, important, unimportant, snp_filter, defining_snps):
    """
    Recursively tries to call each important haplotype with the minimal amount of snps.  
    This method will search for any haplotypes that can be called simply with the addition of one
    more snp, if not it chooses the highest scoring snp possible and will recurse using that snp to search for
    more calls.  Any duplicate haplotypes will be marked as 'Unidentifiable' since they cannot be uniquely called.

    Args:
      important: List of important haplotypes.
      unimportant: List of unimportant haplotypes.
      snp_filter: List of snps to use for the filtering process.
      defining_snps: Dictionary of all haplotypes that have been called.
      
    Returns a dictionary of haplotypes and the snps used to call them.
    """
    
    identifiable = [hap for hap in important if hap not in self.synonyms]
    unidentifiable = [hap for hap in important if hap not in identifiable]
    scored_snps = self.score_snps(identifiable, unimportant, snp_filter)
    brand_new_snp_filter = list(snp_filter)
    added_snps = False
    for snp in scored_snps:
      new_snp_filter = list(snp_filter)
      new_snp_filter.append(snp['snp'])
      call = self.is_callable(identifiable, unimportant, new_snp_filter)
      if call:
        for hap, snps in call.iteritems():
          if hap not in defining_snps:
            new_snp_filter = new_snp_filter + [snp.split(":")[0] for snp in snps]
            new_snp_filter = list(set(new_snp_filter))
            defining_snps[hap] = snps
        added_snps = True

    uncallable_haps = [hap for hap in identifiable if hap not in defining_snps.keys()]
    if not uncallable_haps:
      return defining_snps
    else:
      scored_snps = self.score_snps(uncallable_haps, unimportant, snp_filter)
      for snp in scored_snps:
        new_snp_filter = list(brand_new_snp_filter)
        if not added_snps:
          new_snp_filter.append(snp['snp'])
        result = self.check_all(uncallable_haps, unimportant, new_snp_filter, defining_snps)
        if result:
          defining_snps.update(result)
          return defining_snps
    return []
      
  def collapse_haps(self, hap_defs):
    """
    Method that will collapse a set of snps that can be used to call a haplotype into the
    minimal set that is unique to that haplotype.

    Args:
      hap_defs: Dictionary of haplotypes and snp lists.

    Returns a dictionary of haplotype and the minimal amount of snps that can be used to call that haplotype.
    """
    
    collapsed_haps = {}
    all_snps = list(set([snp.split(":")[0] for snp_list in hap_defs.values() for snp in snp_list]))
    for haplotype, matching_snps in hap_defs.iteritems():
      for i in range(1, len(matching_snps) + 1):
        for combo in itertools.combinations(matching_snps, i):
          matches = self.match(combo)
          if len(matches) == 1 and matches[0] in haplotype:
            if (haplotype in collapsed_haps and len(combo) < len(collapsed_haps[haplotype])) or haplotype not in collapsed_haps:
              collapsed_haps[haplotype] = list(combo)
    return collapsed_haps 
          
  def is_callable(self, important, unimportant, snps):
    """
    Determines if a set of haplotypes is callable based upon the current set of snps.  If the snps
    are sufficient to call a haplotype, the haplotype will be returned in a dictionary with the list
    of snps used to call the haplotype.
    
    Args:
      important: List of important haplotypes.
      unimportant: List of unimportant haplotypes.
      snps: List of snps.

    Returns a dictionary of haplotypes and their corresponding list of snps that can be used to call them.
    """
    
    matching_alleles = self.find_matching(snps)
    snp_map = defaultdict(list)
    good_snp_map = {}
    for haplotype, matching_snps in matching_alleles.iteritems():
      snp_hash = "|".join(sorted(matching_snps))
      snp_map[snp_hash].append(haplotype)
    for snp_hash, haps in snp_map.iteritems():
      if len(haps) == 1 and haps[0] in important:
        good_snp_map[haps[0]] = snp_hash.split("|")
    if good_snp_map:
      return good_snp_map
    return False
      
  def score_snps(self, important_haplos, unimportant_haplos, snp_filter):
    """
    Assign a score to each snp.  This score is meant to minimize the number
    of snps that should be used to define a haplotype.  The score is meant to maximize the
    selectivity of a certain snp and to find those snps that can define a haplotype most easily.

    Args:
      important_haplos: List of important haplotypes.
      unimportant_haplos: List of unimportant haplotypes.
      snp_filter: List of all snps.

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
      snp = i
      if snp in snp_filter:
        continue
      important_score = 0
      for base, count in important_counter.iteritems():
        ratio = count / float(len(important_indexes))
        u_ratio = unimportant_counter[base] / float(len(unimportant_indexes))
        base_score = (1 - ratio) + (1 - u_ratio)
        important_score = important_score + base_score
      all_snp_scores.append({'snp': snp, 'score': important_score})

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
  haplo_mapping = haplo_builder.check_haplos(important_haplos, unimportant_haplos)
  for item, snps in haplo_mapping.iteritems():
    print "%s\t%s" % (item, "\t".join(sorted(snps)))
