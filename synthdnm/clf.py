import pandas as pd
# from sklearn.externals import joblib
import joblib
import os,sys
import numpy as np

def classify_dataframe(df = None, clf = None, ofh = None, mode = "a", keep_fp = True):
    pd.options.mode.chained_assignment = None
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(axis=0,subset=df.columns[12:36])

    # ClippingRankSum (temporary solution)
    df["ClippingRankSum"] = 0

    if df.empty: 
        # print("Empty dataframe.")
        return 0
    X = df[df.columns[12:36]].values
    df["pred"] = clf.predict(X)
    df["prob"] = clf.predict_proba(X)[:,1]

    if keep_fp == False:
        df = df.loc[df["pred"] == 1]

    with open(ofh, mode) as f:
        df.to_csv(f, sep="\t", header = False, index=False)

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

# def classify(ofh_tmp=None,ofh=None,keep_fp=None,pseud=None,vcf=None,make_bed=True,make_vcf=True,fam_fh=None):
def classify(feature_table=None,keep_fp=True,pseud=None,fam_fh=None,snv_clf="snp_100-12-10-2-1-0.0-100.joblib",clf_indel="ssc-jg.half.indel_100-auto-None.joblib"):
    # Get classifiers
    clf = joblib.load(snv_clf)
    clf_indels = joblib.load(clf_indel)
    
    # Make dataframe from input pydnm file
    df = pd.read_csv(feature_table,sep="\t",dtype={"chrom": str})
    df_fam = get_sex(fam_fh)

    # pseud_chrX = pseud["chrX"]
    # pseud_chrX_interval_one = pseud_chrX[0]
    # pseud_chrX_interval_two = pseud_chrX[1]
    # pseud_chrY = pseud["chrY"]
    # pseud_chrY_interval_one = pseud_chrY[0]
    # pseud_chrY_interval_two = pseud_chrY[1]

    
    # def classify_dataframe(df, clf, ofh,pyDNM_header=False, mode="a",keep_fp=True):
    from pathlib import Path
    feature_file_stem = Path(feature_filename).stem
    feature_file_parent = str(Path(feature_filename).parent) + "/"
    feature_file_parent_stem = feature_file_parent + feature_file_stem 
    ofh = feature_file_parent_stem + ".preds.txt"

    with open(feature_filename) as f:
        header_list = f.readline().rstrip().split("\t")
    header_list = header_list + ["sex","pred","prob"]

    df_out = pd.DataFrame(columns = header_list)
    df_out.to_csv(ofh, sep = "\t", index = False)
    
    df['iid']=df['iid'].astype(str)
    df_fam['iid']=df_fam['iid'].astype(str)
    df = pd.merge(df, df_fam, on="iid")
    df=df[~df["chrom"].str.contains("GL*")]
    df["chrom"]=df["chrom"].astype(str)
    df["chrom"] = df["chrom"].apply(lambda s: "chr" + s if not s.startswith("chr") else s)


    df_autosomal_SNV = df.loc[(df["chrom"] != "chrY") & (df["chrom"] != "chrX") &(df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0') & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)]
    df_autosomal_indel = df.loc[(df["chrom"] != "chrY") & (df["chrom"] != "chrX") &(df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1))]


    # df_female_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "2") & (df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0') & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1)]
    # df_female_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "2") &(df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0') & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1))]
    # df_male_nonPAR_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["offspring_gt"]=='1/1') & (df["mother_gt"]=='0/0') & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & ~(df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    # df_male_nonPAR_Y_SNV = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") & (df["offspring_gt"]=='1/1')&(df["father_gt"]=='0/0') & (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & ~(df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]
    # df_male_nonPAR_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["offspring_gt"]=='1/1') & (df["mother_gt"]=='0/0') & ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & ~(df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    # df_male_nonPAR_Y_indel = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") & (df["offspring_gt"]=='1/1') &(df["father_gt"]=='0/0')& ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & ~(df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]
    # df_male_PAR_X_SNV = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & (df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    # df_male_PAR_Y_SNV = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") &(df["offspring_gt"]=='0/1') &(df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& (df["ref"].str.len() == 1) & (df["alt"].str.len() == 1) & (df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]
    # df_male_PAR_X_indel = df.loc[(df["chrom"] == "chrX") & (df["sex"] == "1") & (df["offspring_gt"]=='0/1') & (df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & (df["pos"].between(pseud_chrX_interval_one[0],pseud_chrX_interval_one[1]) | df["pos"].between(pseud_chrX_interval_two[0], pseud_chrX_interval_two[1]))]
    # df_male_PAR_Y_indel = df.loc[(df["chrom"] == "chrY") & (df["sex"] == "1") &(df["offspring_gt"]=='0/1') &(df["mother_gt"]=='0/0') &(df["father_gt"]=='0/0')& ((df["ref"].str.len() != 1) | (df["alt"].str.len() != 1)) & (df["pos"].between(pseud_chrY_interval_one[0],pseud_chrY_interval_one[1]) | df["pos"].between(pseud_chrY_interval_two[0], pseud_chrY_interval_two[1]))]



    classify_dataframe(df = df_autosomal_SNV, clf = clf, ofh = ofh)
    classify_dataframe(df = df_autosomal_indel,clf = clf_indels, ofh = ofh)


    # classify_dataframe(df_female_X_SNV,clf,ofh_new)
    # classify_dataframe(df_female_X_indel,clf_indels,ofh_new)
    # classify_dataframe(df_male_nonPAR_X_SNV,clf_chrX_snps,ofh_new)
    # classify_dataframe(df_male_nonPAR_Y_SNV,clf_chrY_snps,ofh_new)    
    # classify_dataframe(df_male_nonPAR_X_indel,clf_chrX_chrY_indels,ofh_new)    
    # classify_dataframe(df_male_nonPAR_Y_indel,clf_chrX_chrY_indels,ofh_new)
    # classify_dataframe(df_male_PAR_X_SNV,clf,ofh_new)
    # classify_dataframe(df_male_PAR_Y_SNV,clf,ofh_new)
    # classify_dataframe(df_male_PAR_X_indel,clf_indels,ofh_new)
    # classify_dataframe(df_male_PAR_Y_indel,clf_indels,ofh_new)
   
    ofb = feature_file_parent_stem + ".preds.bed"
    fout = open(ofb,"w")
    f = open(ofh,"r")
    make_output_bed(f = f, fout = fout)
    
def make_output_bed(f = None, fout = None):
    f.readline()
    for line in f:
        linesplit = line.rstrip().split("\t")
        chrom,pos,ref,alt,iid,pred,prob = linesplit[0],linesplit[1],linesplit[3],linesplit[4],linesplit[5],linesplit[-2],linesplit[-1]
        pos_0 = str(int(pos)-1)
        pos_1= str(int(pos) + len(ref) - 1)
        ID_column = "{}:{}:{}:{}:{}:{}:{}".format(chrom,pos,ref,alt,iid,pred,prob)
        row = "{}\t{}\t{}\t{}\n".format(chrom,pos_0,pos_1,ID_column)
        fout.write(row)

if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    feature_filename = sys.argv[1]
    ped_filename = sys.argv[2]

    # testing

    # snv_clf_filename = sys.argv[3]
    # indel_clf_filename = sys.argv[4]

    classify(feature_table = feature_filename, fam_fh = ped_filename)
