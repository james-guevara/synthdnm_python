__version__="0.1.0.0"
__usage__="""
 _______  __   __  __    _  _______  __   __  ______   __    _  __   __ 
|       ||  | |  ||  |  | ||       ||  | |  ||      | |  |  | ||  |_|  |
|  _____||  |_|  ||   |_| ||_     _||  |_|  ||  _    ||   |_| ||       |
| |_____ |       ||       |  |   |  |       || | |   ||       ||       |
|_____  ||_     _||  _    |  |   |  |       || |_|   ||  _    ||       |
 _____| |  |   |  | | |   |  |   |  |   _   ||       || | |   || ||_|| |
|_______|  |___|  |_|  |__|  |___|  |__| |__||______| |_|  |__||_|   |_|


Version {}    Authors: Danny Antaki, Aojie Lian, James Guevara    
                   Contact: j3guevar@ucsd.health.edu
---------------------------------------------------------------------------------
    synthdnm  -f <in.fam>  -v  <in.vcf.gz>  [-oLgkVh]
    
necessary arguments:
  
  -v, -vcf    PATH    VCF file
  -f, -fam    PATH    PLINK pedigree (.fam/.ped) file
  
optional arguments:

  -g  --gen    STR     human reference genome version [default: hg38]
  
  -h, --help           show this message and exit
     
""".format(__version__)

import argparse

parser = argparse.ArgumentParser(usage=__usage__)

# Necessary arguments
parser.add_argument("-v","--vcf",required=True)
parser.add_argument("-f","--fam",required=True)
# Optional arguments
parser.add_argument("-g","--gen",required=False,default="hg38",choices=["hg19","hg38"])
args  = parser.parse_args()

vcf_filepath = args.vcf
fam_filepath = args.fam
gen = args.gen

from vcf import parse
parse(vcf_filepath, fam_filepath)
