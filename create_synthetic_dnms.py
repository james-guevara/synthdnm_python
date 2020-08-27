import sys
# Pipeline to create synethic de novo variants from private, inherited variants.
# python create_synthetic_dnms.py <in.vcf.gz>
vcf_filepath = sys.argv[1]
ped_filepath = sys.argv[2]

from pathlib import Path
vcf_stem = Path(vcf_filepath).stem
annotated_vcf_filename = "{}.{}".format(vcf_stem, "annotated.vcf.gz")


import subprocess
# annotate = subprocess.Popen(["bcftools", "annotate", "-I", "%CHROM:%POS0:%END:%REF:%ALT", vcf_filepath], stdout=subprocess.PIPE)
# view = subprocess.Popen(["bcftools", "view", "-e", "N_ALT>1", "-Oz", "-o", annotated_vcf_filename], stdin=annotate.stdout)
# annotate.stdout.close()

ped_filepath = sys.argv[2]


up = "{}."

double_ids = subprocess.Popen(["plink", "--vcf", annotated_vcf_filename, "--double-id", "--make-bed", "--out", vcf_stem])
update_ids = subprocess.Popen["plink", vcf_stem, "--update-ids", , "--make-bed", "--out", vcf_stem])
update_sex_parents = subprocess.Popen["plink", vcf_stem, "--update-ids", , "--make-bed", "--out", vcf_stem])

