def build_synthdnm():
    import sys
    # Pipeline to create synthetic de novo variants from private, inherited variants.
    # python create_synthetic_dnms.py <in.vcf.gz> <in.ped>
    vcf_filepath = sys.argv[1]
    
    from pathlib import Path
    # Gets the stem of the filename (removes .vcf.gz or .vcf extension)
    vcf_stem = Path(vcf_filepath).stem
    if vcf_stem.endswith(".vcf"):
        vcf_stem = Path(vcf_stem).stem
    
    
    ped_filepath = sys.argv[2]
    print(ped_filepath)
    
    # Create plink files (using ped_stem)
    ped_stem = Path(ped_filepath).stem
    ped_parent = str(Path(ped_filepath).parent) + "/"
    
    upid = "{}.{}".format(ped_parent + ped_stem,"update.ids")
    uppa = "{}.{}".format(ped_parent + ped_stem,"update.parents")
    upsx = "{}.{}".format(ped_parent + ped_stem,"update.sex")
    phen = "{}.{}".format(ped_parent + ped_stem,"update.tdt_all_case")


    def preprocess_ped(ped_filepath,upid,uppa,upsx,phen):
    
        fout_upid = open(upid,"w")
        fout_uppa = open(uppa,"w")
        fout_upsx = open(upsx,"w")
        fout_phen = open(phen,"w")
    
        f = open(ped_filepath,"r")
        for line in f:
            linesplit = line.rstrip().split("\t")
            fid,iid,father_iid,mother_iid,sex,phen = linesplit[0],linesplit[1],linesplit[2],linesplit[3],linesplit[4],linesplit[5]
            fout_upid.write("{} {} {} {}\n".format(iid,iid,fid,iid))
            fout_uppa.write("{} {} {} {}\n".format(fid,iid,father_iid,mother_iid))
            fout_upsx.write("{} {} {}\n".format(fid,iid,sex))
            # Generate phen file
            if not (father_iid == "0" and mother_iid == "0"): # Skip parents
                fout_phen.write("{} {} 2\n".format(fid,iid)) # In order to obtain transmission counts from the controls, we treat controls as cases (the PLINK TDT seems to require this)
        f.close()
        fout_upid.close()
        fout_uppa.close()
        fout_upsx.close()
        fout_phen.close()
    
    preprocess_ped(ped_filepath,upid,uppa,upsx,phen)
    
    
    vcf_parent = str(Path(vcf_filepath).parent) + "/"
    vcf_parent_stem = vcf_parent + vcf_stem
    annotated_vcf_filename = "{}.{}".format(vcf_parent + vcf_stem, "annotated.vcf.gz")
    
    import subprocess
    
    annotate = subprocess.Popen(["bcftools", "annotate", "-I", "%CHROM:%POS0:%END:%REF:%ALT", vcf_filepath], stdout=subprocess.PIPE)
    view = subprocess.check_call(["bcftools", "view", "-e", "N_ALT>1", "-Oz", "-o", annotated_vcf_filename], stdin=annotate.stdout)
    # annotate.stdout.close()
    
    
    double_ids = subprocess.check_call(["plink", "--vcf", annotated_vcf_filename, "--double-id", "--make-bed", "--out", vcf_parent_stem, "--vcf-half-call", "m", "--allow-extra-chr"])
    update_ids = subprocess.check_call(["plink", "--bfile", vcf_parent_stem, "--allow-extra-chr", "--update-ids", upid, "--make-bed", "--out", vcf_parent_stem + ".up"])
    update_sex_parents = subprocess.check_call(["plink", "--bfile", vcf_parent_stem + ".up", "--allow-extra-chr", "--update-sex", upsx, "--update-parents", uppa, "--make-bed", "--out", vcf_parent_stem + ".fin"])
    
    freq = subprocess.check_call(["plink", "--bfile", vcf_parent_stem + ".fin", "--allow-extra-chr", "--freq", "counts", "--out", vcf_parent_stem])
    tdt = subprocess.check_call(["plink", "--bfile", vcf_parent_stem + ".fin", "--allow-extra-chr", "--tdt", "poo", "--pheno", phen, "--out", vcf_parent_stem])
    
    from make_private_VCF import make_private_vcf
    priv_inh_vcf_filename = make_private_vcf(annotated_vcf_filename, vcf_parent_stem)
    
    bgzip = subprocess.check_call(["bgzip", priv_inh_vcf_filename])
    tabix = subprocess.check_call(["tabix", priv_inh_vcf_filename + ".gz"])
    
    from swap import swap_ped
    swapped_ped_absolute_path = swap_ped(ped_filepath)
    
    # Make output file for training set
    private_inherited_vcf_absolute_path = vcf_parent_stem + ".annotated.private.inherited.vcf.gz"
    
    fout = open(vcf_parent_stem + ".training_set.txt", "w")
    from vcf import parse
    info_keys = ["VQSLOD","BaseQRankSum","FS","SOR","MQ","MQRankSum","QD","ReadPosRankSum"]
    # Get positive training examples
    parse(private_inherited_vcf_absolute_path, swapped_ped_absolute_path, info_keys,fout=fout,training_examples="1")
    # parse(private_inherited_vcf_absolute_path, swapped_ped_absolute_path, info_keys)
    # Get negative training examples
    parse(vcf_filepath, ped_filepath, info_keys, fout=fout, training_examples="0")

if __name__=="__main__":
    build_synthdnm()

