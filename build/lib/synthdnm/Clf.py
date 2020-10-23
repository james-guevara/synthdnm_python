import pandas as pd
from sklearn.externals import joblib
from Backend import get_path
import os,sys
import numpy as np

def classify_dataframe(df, clf, ofh,pyDNM_header=False, mode="a",keep_fp=True):
    pd.options.mode.chained_assignment = None
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(axis=0,subset=df.columns[12:36])
    if df.empty: 
        # print("Empty dataframe.")
        return 0
    X = df[df.columns[12:36]].values
    df["pred"] = clf.predict(X)
    df["prob"] = clf.predict_proba(X)[:,1]
    if keep_fp == False:
        df = df.loc[df["pred"] == 1]
    with open(ofh, mode) as f:
        df.to_csv(f, sep="\t",header=pyDNM_header, index=False)

def get_sex(fam_fh):
    fam = open(fam_fh, "r")
    fam_dict = {}
    for line in fam:
        linesplit = line.rstrip().split("\t")
        iid = linesplit[1]
        sex = linesplit[4]
        fam_dict[iid] = sex 
    df = pd.Series(fam_dict).to_frame("sex")
    df["iid"] = df.index
    df.reset_index(inplace=True)
    df.drop(columns=["index"],inplace=True)
    return df

def classify(ofh_tmp=None,ofh=None,keep_fp=None,pseud=None,vcf=None,make_bed=True,make_vcf=True,fam_fh=None):
    ofh_new = ofh
    
    # fam 
    #fam_fh = "/home/a1lian/reach_ssc1-4.fam"
    df_fam = get_sex(fam_fh)
    pseud_chrX = pseud["chrX"]
    pseud_chrX_interval_one = pseud_chrX[0]
    pseud_chrX_interval_two = pseud_chrX[1]
    pseud_chrY = pseud["chrY"]
    pseud_chrY_interval_one = pseud_chrY[0]
    pseud_chrY_interval_two = pseud_chrY[1]
    # Get classifiers
    # OLD CLFS:
    snv_clf = get_path()+'/pydnm.snv.clf.joblib'
    indels_clf = get_path()+'/pydnm.indels.clf.joblib'

    # NEW CLFS:
    # snv_clf = get_path()+'/ssc1-1-snp-clf-2020-08-09.joblib'
    # indels_clf = get_path()+'/ssc1-1-indel-clf-2020-08-09.joblib'

    snv_chrX_clf=get_path()+'/chrX_training_snps.joblib'
    snv_chrY_clf= get_path()+'/chrY_training_snps.joblib'
    indels_chrX_chrY_clf= get_path()+'/chrX_chrY_training_indels.joblib'
    if not os.path.isfile(snv_clf):
            sys.stderr.write('FATAL ERROR: {} CLASSIFIER NOT FOUND\n'.format(snv_clf))
            sys.exit(1)
            
    clf = joblib.load(snv_clf)
    clf_indels = joblib.load(indels_clf)
    clf_chrX_snps = joblib.load(snv_chrX_clf)
    clf_chrY_snps = joblib.load(snv_chrY_clf)
    clf_chrX_chrY_indels = joblib.load(indels_chrX_chrY_clf)
    # Make dataframe from input pydnm file
    df = pd.read_csv(ofh_tmp,sep="\t",dtype={"chrom": str})
    
    # Filter original dataframe
    # df = df.loc[(df["offspring_gt"] != "0/0")] 
    
    df['iid']=df['iid'].astype(str)
    df_fam['iid']=df_fam['iid'].astype(str)
    df = pd.merge(df, df_fam, on="iid")
    df=df[~df["chrom"].str.contains("GL*")]
    df["chrom"]=df["chrom"].astype(str)
    df["chrom"] = df["chrom"].apply(lambda s: "chr" + s if not s.startswith("chr") else s)
    df_autosomal_SNV = df.loc[(df["chrom"] != "chrY") & (df["chrom"] != "chrX") &(df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0') & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)]
    df_autosomal_indel = df.loc[(df["chrom"] != "chrY") & (df["chrom"] != "chrX") &(df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1))]
    df_female_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "2") & (df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0') & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)]
    df_female_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "2") &(df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0') & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1))]
    df_male_nonPAR_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["offspring_gt"]=='1/1') & (df["mother_gt"]=='0/0') & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & ~(df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    df_male_nonPAR_Y_SNV = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") & (df["offspring_gt"]=='1/1')&(df["father_gt"]=='0/0') & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & ~(df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]
    df_male_nonPAR_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["offspring_gt"]=='1/1') & (df["mother_gt"]=='0/0') & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & ~(df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    df_male_nonPAR_Y_indel = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") & (df["offspring_gt"]=='1/1') &(df["father_gt"]=='0/0')& ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & ~(df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]
    df_male_PAR_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & (df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    df_male_PAR_Y_SNV = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") &(df["offspring_gt"]=='0/1') &(df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & (df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]
    df_male_PAR_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & (df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    df_male_PAR_Y_indel = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") &(df["offspring_gt"]=='0/1') &(df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & (df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]

    classify_dataframe(df_autosomal_SNV,clf,ofh_new,False,"w")
    classify_dataframe(df_autosomal_indel,clf_indels,ofh_new)
    classify_dataframe(df_female_X_SNV,clf,ofh_new)
    classify_dataframe(df_female_X_indel,clf_indels,ofh_new)
    classify_dataframe(df_male_nonPAR_X_SNV,clf_chrX_snps,ofh_new)
    classify_dataframe(df_male_nonPAR_Y_SNV,clf_chrY_snps,ofh_new)    
    classify_dataframe(df_male_nonPAR_X_indel,clf_chrX_chrY_indels,ofh_new)    
    classify_dataframe(df_male_nonPAR_Y_indel,clf_chrX_chrY_indels,ofh_new)
    classify_dataframe(df_male_PAR_X_SNV,clf,ofh_new)
    classify_dataframe(df_male_PAR_Y_SNV,clf,ofh_new)
    classify_dataframe(df_male_PAR_X_indel,clf_indels,ofh_new)
    classify_dataframe(df_male_PAR_Y_indel,clf_indels,ofh_new)
   
    df = pd.read_csv(ofh_new,sep="\t",header=None)
    df.columns=['chrom','pos','id','ref','alt','iid','offspring_gt','father_gt','mother_gt','nalt','filter','qual','parent_ar_max','parent_ar_min','offspring_ar','parent_dp_max','parent_dp_min','offspring_dp','parent_dnm_pl_max','parent_dnm_pl_min','parent_inh_pl_max','parent_inh_pl_min','offspring_dnm_pl','offspring_inh_pl','parent_gq_max','parent_gq_min','offspring_gq','VQSLOD','ClippingRankSum','BaseQRankSum','FS','SOR','MQ','MQRankSum','QD','ReadPosRankSum','sex','pred','prob']
    with open(ofh_new,'w') as f:
        df.to_csv(f, sep="\t", index=False)
    if make_bed: 
        ofb = make_output_bed(ofh_new)

    
def make_output_bed(ofh_new):
    ofb =ofh_new+".bed"
    fout = open(ofb,"w")
    f = open(ofh_new,"r")
    f.readline()
    dnm_bed = []
    for line in f:
        linesplit = line.rstrip().split("\t")
        chrom = linesplit[0]
        pos = linesplit[1]
        ref = linesplit[3]
        alt = linesplit[4]
        ###SNPs
        if len(ref) == 1 and len(alt) == 1: 
            iid = linesplit[5]
            pred = linesplit[-2]
            prob = linesplit[-1]
            pos_0 = str(int(pos)-1)
            pos_1 = pos
            ID_col = "{}:{}:{}:{}:{}:{}:{}".format(chrom,pos,ref,alt,iid,pred,prob)
            newline = "{}\t{}\t{}\t{}\n".format(chrom,pos_0,pos_1,ID_col)
            fout.write(newline)
            dnm_bed.append(newline)
        ###INDELs
        else:
            iid = linesplit[5]
            pred = linesplit[-2]
            prob = linesplit[-1]
            pos_0 = str(int(pos)-1)
            pos_1= str(int(pos) + len(ref) - 1)
            ID_col = "{}:{}:{}:{}:{}:{}:{}".format(chrom,pos,ref,alt,iid,pred,prob)
            newline = "{}\t{}\t{}\t{}\n".format(chrom,pos_0,pos_1,ID_col)
            fout.write(newline)
            dnm_bed.append(newline)

    return ofb 
 
