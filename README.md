# HaploBuilder

This project will read a list of known haplotypes, and a list of "good" haplotypes.  It will then create a translation of the minimum SNPs required to evaluate each "good" haplotype.

```
Usage: haplo_builder.py [options]

Options:
  -h, --help            show this help message and exit
  -f FILENAME, --file=FILENAME
                        File with good haplotypes to read.  One haplotype per
                        line. (Required)
  -t HAP_FILE, --types=HAP_FILE
                        File with all haplotypes to read.
  -d, --debug           Print debug messages.

```