import sys
# Uses the pedigree (.ped) file to create files that will be used by PLINK. 
ped_filepath = sys.argv[1]
# Gets the stem of the pedigree file (retaining only the base filename and removing the file extension.
from pathlib import Path
ped_stem = Path(ped_filepath).stem

upid = "{}.{}".format(ped_stem,"update.ids")
uppa = "{}.{}".format(ped_stem,"update.parents")
upsx = "{}.{}".format(ped_stem,"update.sex")
phen = "{}.{}".format(ped_stem,"update.tdt_all_case")

fout_upid = open(upid,"w")
fout_uppa = open(upid,"w")
fout_upsx = open(upid,"w")
fout_phen = open(upid,"w")
f = open(ped_filepath,"r")
for line in f:
    linesplit = line.rstrip().split("\t")
    fid,iid,father_iid,mother_iid,sex,phen = linesplit[0],linesplit[1],linesplit[2],linesplit[3],linesplit[4]
    fout_upid.write("{} {} {} {}".format(iid,iid,fid,iid))
    fout_uppa.write("{} {} {} {}".format(fid,iid,father_iid,mother_iid))
    fout_upsx.write("{} {} {}".format(fid,iid,sex))
    # Generate phen file
    if not (father_iid == "0" and mother_iid == "0"): # Skip parents
        fout_phen.write("{} {} 2".format(fid,iid)) # In order to obtain transmission counts from the controls, we treat controls as cases (the PLINK TDT seems to require this)
f.close()
fout_upid.close()
fout_uppa.close()
fout_upsx.close()
fout_phen.close()
# Return the names of the created files
return [upid,uppa,upsx,phen]
