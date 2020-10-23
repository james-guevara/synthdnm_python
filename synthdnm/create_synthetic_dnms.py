import sys
# Pipeline to create synethic de novo variants from private, inherited variants.
# python create_synthetic_dnms.py <in.vcf.gz>
vcf_filepath = sys.argv[1]

from pathlib import Path
vcf_stem = Path(vcf_filepath).stem
if vcf_stem.endswith(".vcf"):
    vcf_stem = Path(vcf_stem).stem

annotated_vcf_filename = "{}.{}".format(vcf_stem, "annotated.vcf.gz")

ped_filepath = sys.argv[2]
ped_stem = Path(ped_filepath).stem
upid = "{}.{}".format(ped_stem,"update.ids")
uppa = "{}.{}".format(ped_stem,"update.parents")
upsx = "{}.{}".format(ped_stem,"update.sex")
phen = "{}.{}".format(ped_stem,"update.tdt_all_case")


import subprocess
'''
annotate = subprocess.Popen(["bcftools", "annotate", "-I", "%CHROM:%POS0:%END:%REF:%ALT", vcf_filepath], stdout=subprocess.PIPE)
view = subprocess.check_call(["bcftools", "view", "-e", "N_ALT>1", "-Oz", "-o", annotated_vcf_filename], stdin=annotate.stdout)
# annotate.stdout.close()


double_ids = subprocess.check_call(["plink", "--vcf", annotated_vcf_filename, "--double-id", "--make-bed", "--out", vcf_stem, "--vcf-half-call", "m", "--allow-extra-chr"])
update_ids = subprocess.check_call(["plink", "--bfile", vcf_stem, "--allow-extra-chr", "--update-ids", upid, "--make-bed", "--out", vcf_stem + ".up"])
update_sex_parents = subprocess.check_call(["plink", "--bfile", vcf_stem + ".up", "--allow-extra-chr", "--update-sex", upsx, "--update-parents", uppa, "--make-bed", "--out", vcf_stem + ".fin"])

freq = subprocess.check_call(["plink", "--bfile", vcf_stem + ".fin", "--allow-extra-chr", "--freq", "counts", "--out", vcf_stem])
tdt = subprocess.check_call(["plink", "--bfile", vcf_stem + ".fin", "--allow-extra-chr", "--tdt", "poo", "--pheno", phen, "--out", vcf_stem])
'''
# from make_private_VCF import make_private_vcf
# make_private_vcf(annotated_vcf_filename, vcf_stem)

from vcf import parse_private_inherited
parse_private_inherited(vcf_stem + ".annotated.private.inherited.vcf", ped_filepath) 
