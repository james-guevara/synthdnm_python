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
    synthdnm-classify  -f <in.fam>  -d  <in.features.txt>
    
necessary arguments:
  
  -f, --fam         PATH    PLINK pedigree (.fam/.ped) file
  -d, --features    PATH    feature file
  
optional arguments:

  -s, --snp_classifier                       PATH    path to snp classifier joblib file
  -l, --indel_classifier                     PATH    path to indel classifier joblib file 
  -p, --keep_all_putative_dnms                       flag that retains all putative dnms (and their scores) in the output files 

  
  -h, --help           show this message and exit
     
""".format(__version__)



import pandas as pd
# from sklearn.externals import joblib
import joblib
import os,sys
import numpy as np

def classify_dataframe(df = None, clf = None, ofh = None, mode = "a", keep_fp = False, features = None):
    pd.options.mode.chained_assignment = None
    df = df.replace([np.inf, -np.inf], np.nan)
    df = df.dropna(axis=0,subset=df.columns[12:36])

    # ClippingRankSum (temporary solution)
    df["ClippingRankSum"] = 0

    if df.empty: 
        # print("Empty dataframe.")
        return 0

    X = df[features].to_numpy()

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

def classify(feature_table=None,keep_fp=False,pseud=None,fam_fh=None,clf_snv="snp_100-12-10-2-1-0.0-100.joblib",clf_indel="indel_1000-12-25-2-1-0.0-100.joblib"):
    # Get classifiers
    clf = joblib.load(clf_snv)
    clf_indels = joblib.load(clf_indel)
    
    # Make dataframe from input pydnm file
    df = pd.read_csv(feature_table,sep="\t",dtype={"chrom": str})
    # Get the list of features
    columns = list(df.columns)
    non_features = ['chrom', 'pos', 'ID', 'ref', 'alt', 'iid', 'offspring_gt', 'father_gt', 'mother_gt', 'nalt', 'filter', 'qual']
    features = [elem for elem in columns if elem not in non_features]

    df_fam = get_sex(fam_fh)

    # pseud_chrX = pseud["chrX"]
    # pseud_chrX_interval_one = pseud_chrX[0]
    # pseud_chrX_interval_two = pseud_chrX[1]
    # pseud_chrY = pseud["chrY"]
    # pseud_chrY_interval_one = pseud_chrY[0]
    # pseud_chrY_interval_two = pseud_chrY[1]

    
    from pathlib import Path
    feature_filename = feature_table
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



    classify_dataframe(df = df_autosomal_SNV, clf = clf, ofh = ofh, features = features, keep_fp = keep_fp)
    classify_dataframe(df = df_autosomal_indel,clf = clf_indels, ofh = ofh, features = features, keep_fp = keep_fp)


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
    import argparse
    parser = argparse.ArgumentParser(usage=__usage__)
    # Necessary arguments
    parser.add_argument("-d","--features",required=True)
    parser.add_argument("-f","--fam",required=True)
    # Optional arguments
    parser.add_argument("-s","--snp_classifier",required=False)
    parser.add_argument("-i","--indel_classifier",required=False)
    parser.add_argument('-p',"--keep_all_putative_dnms", action='store_true')    
    args  = parser.parse_args()


    # feature_filename = sys.argv[1]
    feature_filename = args.features
    # ped_filename = sys.argv[2]
    ped_filename = args.fam

    if args.snp_classifier:
        snv_clf_filename = args.s
    else: snv_clf_filename = "snp_100-12-10-2-1-0.0-100.joblib"
    if args.indel_classifier:
        indel_clf_filename = args.i
    else: indel_clf_filename = "indel_1000-12-25-2-1-0.0-100.joblib"

    keep_fp = args.keep_all_putative_dnms

    classify(feature_table = feature_filename, fam_fh = ped_filename, clf_snv = snv_clf_filename, clf_indel = indel_clf_filename, keep_fp = keep_fp)
