from collections import OrderedDict

num_parents = {}
# skip families that don't have both parents
def skip_families(ped_filename):
    num_parents = {}
    f = open(ped_filename,"r")
    for line in f:
        linesplit = line.rstrip().split("\t")
        fid, iid, iid_father, iid_mother, sex = linesplit[0],linesplit[1],linesplit[2],linesplit[3],linesplit[4]
        if iid_father == "0" and iid_mother == "0": # this means it's one of the parents
            if fid not in num_parents: num_parents[fid] = 1
            else: num_parents[fid] += 1
    return num_parents


# fam/ped file: FID, IID, IID_Father, IID_Mother, Sex, Phenotype
def make_family_ordered_dict(ped_filename, num_parents):
    fids_od = OrderedDict()
    f = open(ped_filename,"r")
    for line in f:
        linesplit = line.rstrip().split("\t")
        fid, iid, iid_father, iid_mother, sex = linesplit[0],linesplit[1],linesplit[2],linesplit[3],linesplit[4]
        if iid_father == "0" and iid_mother == "0": # this means it's one of the parents
            if num_parents[fid] != 2: continue # don't use this family if there aren't 2 parents
            if fid not in fids_od: fids_od[fid] = [None]*2
            if sex == "1": # male (father)
                fids_od[fid][0] = iid
            elif sex == "2": # female (mother)
                fids_od[fid][1] = iid
    return fids_od
    
import sys
ped_filename = sys.argv[1]
num_parents = skip_families(ped_filename)


fids_od = make_family_ordered_dict(ped_filename, num_parents)
fid1_key = None
fid1_val = None
fid2_key = None
fid2_val = None
swapped_parents = []
new_families_set = set() # the set of families used in the swap
new_families_dict = {}

for i, (key,val) in enumerate(fids_od.items()):
    if not fid1_key: # if fid1_key is None, make this family fid1
        fid1_key = key
        fid1_val = val
        continue
    elif not fid2_key: # if fid2_key is None, make this family fid2
        fid2_key = key
        fid2_val = val
    if fid1_key and fid2_key: # swap the parents of these 2 families and then makes the keys none
        swapped_parents.append("{}\t{}\t{}\t{}\t{}\t{}".format(fid1_key,fid2_val[0], "0", "0", "1", "0"))
        swapped_parents.append("{}\t{}\t{}\t{}\t{}\t{}".format(fid1_key,fid2_val[1], "0", "0", "2", "0"))
        swapped_parents.append("{}\t{}\t{}\t{}\t{}\t{}".format(fid2_key,fid1_val[0], "0", "0", "1", "0"))
        swapped_parents.append("{}\t{}\t{}\t{}\t{}\t{}".format(fid2_key,fid1_val[1], "0", "0", "2", "0"))
        new_families_set.add(fid1_key)
        new_families_set.add(fid2_key)
        new_families_dict[fid1_key] = fid2_val
        new_families_dict[fid2_key] = fid1_val
        fid1_key = None
        fid2_key = None

def print_new_ped(ped_filename, swapped_parents,new_families_set):
    f = open(ped_filename,"r")
    for line in f:
        linesplit = line.rstrip().split("\t")
        fid, iid, iid_father, iid_mother, sex, phen = linesplit[0],linesplit[1],linesplit[2],linesplit[3],linesplit[4],linesplit[5]
        if fid not in new_families_set: continue 
        if not (iid_father == "0" and iid_mother == "0"): # this means it's *not* one of the parents
            iid_father_new, iid_mother_new = new_families_dict[fid][0], new_families_dict[fid][1]
            print("{}\t{}\t{}\t{}\t{}\t{}".format(fid,iid,iid_father_new,iid_mother_new,sex,phen))
    for elem in swapped_parents:
        print(elem)

print_new_ped(ped_filename,swapped_parents,new_families_set)
