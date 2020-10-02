import re
# Splits a line by spaces or tabs, returns a list object
def tokenize(line):
    linesplit = line.rstrip().split("\t")
    if len(linesplit) == 1: linesplit = line.rstrip().split(" ")
    else: return linesplit

# Returns dictionary containing information of pedigree file
def process_ped(fam_filepath):
    ped_dictionary = {} # ped_dictionary[iid] = (fid, iid, father_iid, mother_iid, sex, phen)
    with open(fam_filepath) as f:
        for line in f:
            linesplit = tokenize(line)
            fid,iid,father_iid,mother_iid,sex,phen = linesplit[0],linesplit[1],linesplit[2],linesplit[3],linesplit[4],linesplit[5]
            if sex != "1" and sex != "2": continue
            if father_iid == "0" and mother_iid == "0": continue
            ped_dictionary[iid] = (fid,iid,father_iid,mother_iid,sex,phen)
    return ped_dictionary

# Convert spaces to tabs
def tabbit(line):
    linesplit = re.split("\s+", line.rstrip())
    return "\t".join(linesplit[1:])
