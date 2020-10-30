__version__="0.1.0.1"
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

def run_synthdnm():
    import argparse
    
    parser = argparse.ArgumentParser(usage=__usage__)
    
    # Necessary arguments
    parser.add_argument("-v","--vcf",required=True)
    parser.add_argument("-f","--fam",required=True)
    # Optional arguments
    parser.add_argument("-g","--gen",required=False,default="hg38",choices=["hg19","hg38"])
    parser.add_argument("-i","--info",required=False)
    args  = parser.parse_args()
    
    vcf_filepath = args.vcf
    fam_filepath = args.fam
    gen = args.gen
    
    info_keys = []
    if args.info: 
        f = open(args.info,"r")
        for line in f:
           info_keys.append(line.rstrip()) 
    else: info_keys = ["VQSLOD","ClippingRankSum","BaseQRankSum","FS","SOR","MQ","MQRankSum","QD","ReadPosRankSum"]

    import sys
    from pathlib import Path
    # Gets the stem of the filename (removes .vcf.gz or .vcf extension)
    vcf_stem = Path(vcf_filepath).stem
    if vcf_stem.endswith(".vcf"):
        vcf_stem = Path(vcf_stem).stem
    vcf_parent = str(Path(vcf_filepath).parent) + "/"
    vcf_parent_stem = vcf_parent + vcf_stem
    fout = open(vcf_parent_stem + ".synthdnm.features.txt","w")

    from vcf import parse
    parse(vcf_filepath, fam_filepath, info_keys=info_keys, fout=fout)
    

if __name__=="__main__":
    run_synthdnm()
