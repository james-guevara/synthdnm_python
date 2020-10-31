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
  
  -v, --vcf    PATH    VCF file
  -f, --fam    PATH    PLINK pedigree (.fam/.ped) file
  
optional arguments:

  -e, --extract_features                             flag that disable classification (if you only want to extract features)
  -s, --snp_classifier                       PATH    path to snp classifier joblib file [default is pretrained classifier] 
  -l, --indel_classifier                     PATH    path to indel classifier joblib file [default is pretrained classifier]
  -g  --gen                                  STR     human reference genome version [default: hg38]
  
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

    parser.add_argument("-e", "--extract_features", action="store_true")
    parser.add_argument("-s","--snp_classifier",required=False)
    parser.add_argument("-l","--indel_classifier",required=False)

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
    feature_filename = vcf_parent_stem + ".synthdnm.features.txt"
    fout = open(feature_filename,"w")

    from .vcf import parse
    parse(vcf_filepath, fam_filepath, info_keys=info_keys, fout=fout)
    fout.close()

    if args.extract_features: return

    if args.snp_classifier:
        snv_clf_filename = args.s
    else: snv_clf_filename = "snp_100-12-10-2-1-0.0-100.joblib"
    if args.indel_classifier:
        indel_clf_filename = args.l
    else: indel_clf_filename = "indel_1000-12-25-2-1-0.0-100.joblib"

    from .clf import classify 
    classify(feature_table = feature_filename, fam_fh = fam_filepath, clf_snv = snv_clf_filename, clf_indel = indel_clf_filename)

if __name__=="__main__":
    run_synthdnm()
