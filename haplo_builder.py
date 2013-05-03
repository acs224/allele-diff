import string
import operator
import logging
from collections import defaultdict
from operator import itemgetter

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
    
    self.haplos = self.parse_haplos(haplo_file)

  def parse_haplos(self, haplo_file):
    """Parse the haplo_file into a list of POS:BASE combinations.
    
    Args:
      haplo_file: File to parse
      
    Returns a dict where haplotype is the key and the value is the list of SNP definitions for that haplotype.
      
    """
    haplo_definitions = {}
    with open(haplo_file) as f:
      for i, line in enumerate(f):
        if i < self.HEADER_ROWS:
          continue
        allele_def = line.split()
        allele_string = allele_def[1:]
        # Ignore positions that are NO_CALL or EXON definitions.  
        # NO_CALL can be ignored since it does not contain any basepair definition.
        haplo_definitions[allele_def[0]] = ["%s:%s" % (i, base) for i, base in enumerate(allele_string) if (self.EXON != base and self.NO_CALL != base)]
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
    logging.debug("In check_haplos.  Important: %s Filter: %s" % (important, snp_filter))
    important_haplos = self.filter_haplos(important, snp_filter)
    unimportant_haplos = self.filter_haplos(unimportant, snp_filter, True)
    remaining_important = [haplo for haplo in important if haplo not in important_haplos]
    filtered_snps = [temp_snp['snp'] for temp_snp in snp_filter]
    if not important or len(important_haplos) == 0:
      return []
    logging.debug("Before recursion check.  snp_filter: %s important %s unimportant %s" % (snp_filter, important_haplos, unimportant_haplos))
    if len(snp_filter) > 0 and len(important_haplos) == 1 and len(unimportant_haplos) == 0:
      logging.debug("About to recurse.  important: %s remaining_important %s" % (important_haplos[0], remaining_important))
      return [{'haplotype': important_haplos[0], 'snps': filtered_snps}] + self.check_haplos(remaining_important, unimportant, [])
      
    all_haplos = important_haplos + unimportant_haplos
    snp_map = defaultdict(lambda: defaultdict(list))
    for haplo in all_haplos:
      seq = self.haplos[haplo]
      for item in seq:
        snp_map[item]['present'].append(haplo)
        
    all_snp_scores = self.score_snps(important_haplos, unimportant_haplos, snp_map)
    for to_test in all_snp_scores:
      if to_test['snp'] not in filtered_snps:
        snp_filter.append(to_test)
        return self.check_haplos(important, unimportant, snp_filter)
    return []
  
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
    to_filter = [snp['snp'] for snp in snp_filter]
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
    
  def score_snps(self, important_haplos, unimportant_haplos, snp_map):
    """Assign a score to each snp.  This score is meant to minimize the number
    of snps that should be used to define a haplotype.  The current formula is:
    
    score = (ratio of snp in important haplos) - (ratio of snp in unimportant haplos) + (weighted score to favor snps over reference)
    
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
        raw_score = important_ratio - unimportant_ratio + (0 if "_" in snp else 1)
        all_snp_scores.append({'snp': snp, 'score': raw_score})
        
    return sorted(all_snp_scores, key=itemgetter('score'), reverse=True)
  
if __name__ == "__main__":
  #logging.basicConfig(level=logging.DEBUG)
  haplo_builder = HaploBuilder()
  important_haplos = ['B*51240301', 'B*55070101']
  unimportant_haplos = [hap for hap in haplo_builder.haplos.keys() if hap not in important_haplos]
  haplo_mapping = haplo_builder.check_haplos(important_haplos, unimportant_haplos, [])
  print haplo_mapping
