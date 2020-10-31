import sys
from .backend import tabbit

def get_private_IDs(frq_counts_fh):
    private_IDs = set()
    f = open(frq_counts_fh,"r")
    f.readline()
    for line in f:
        ID = tabbit(line).rstrip().split("\t")[1]
        C1 = int(tabbit(line).rstrip().split("\t")[4])
        if C1 == 1: private_IDs.add(ID)
    return private_IDs

def get_transmitted_variants(tdt_poo_fh):
    trans_IDs = set()
    f = open(tdt_poo_fh,"r")
    f.readline()
    for line in f:
        linesplit = tabbit(line).rstrip().split("\t")
        ID = linesplit[1] 
        T_U_PAT = linesplit[3].split(":")
        T_U_MAT = linesplit[6].split(":")
        if "." in T_U_PAT[0] or "." in T_U_PAT[1] or "." in T_U_MAT[0] or "." in T_U_MAT[1]: continue
        T = int(T_U_PAT[0]) + int(T_U_MAT[0]) # transmitted
        if T == 1: trans_IDs.add(ID)
    return trans_IDs

def make_private_vcf(annotated_vcf_fh, vcf_prev_stem):
    private_IDs = get_private_IDs(vcf_prev_stem + ".frq.counts")
    trans_IDs = get_transmitted_variants(vcf_prev_stem + ".tdt.poo")
    # Intersect to get private, transmitted variants?
    private_trans_IDs = private_IDs.intersection(trans_IDs)

    from pathlib import Path
    vcf_stem = Path(annotated_vcf_fh).stem
    if vcf_stem.endswith(".vcf"):
        vcf_stem = Path(vcf_stem).stem

    vcf_parent = str(Path(annotated_vcf_fh).parent) + "/"
    vcf_parent_stem = vcf_parent + vcf_stem

    priv_inh_vcf_filename = vcf_parent + vcf_stem + ".private.inherited.vcf"

    # fout = open(vcf_stem + ".private.vcf", "w")
    fout = open(priv_inh_vcf_filename, "w")
    if annotated_vcf_fh.endswith(".gz"):
        import gzip
        annotated_vcf_fh = gzip.open(annotated_vcf_fh,"rt",9)
    else:
        annotated_vcf_fh = open(annotated_vcf_fh,"r")
    for line in annotated_vcf_fh:
        if line.startswith("#"): fout.write(line) # Header
        else:
            ID = line.rstrip().split("\t")[2]
            if ID in private_trans_IDs: fout.write(line) # Private variants
    return priv_inh_vcf_filename
