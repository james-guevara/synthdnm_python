from backend import tokenize
def index_samples(CHROM_header_line, ped_dictionary):
    """Processes CHROM header line and stores column indices with sample ID"""
    iid_indices = {}  # iid_indices[iid] = index
    linesplit = tokenize(CHROM_header_line) 
    for i in range(9, len(linesplit)):
        iid = linesplit[i] # The sample ID (IID) at the ith column
        iid_indices[iid] = i
    return iid_indices

def index_format(FORMAT):
    """Determines indices of FORMAT values"""
    format_indices = {}
    fields = FORMAT.split(":")
    for i in range(0, len(fields)): 
        field = fields[i]
        format_indices[field] = i
    return format_indices 

def get_feature(iid_feature_index, format_index):
    return iid_feature_index.split(":")[format_index]

# Make more general
def load_genotypes(linesplit, iid_indices, format_indices):
    """Stores genotypes"""
    iid_genotypes = {} # iid_genotypes[iid] = "0/1" (example)
    GT_index = format_indices["GT"]
    for iid,iid_index in iid_indices.items():
        iid_genotypes[iid] = linesplit[iid_index].split(":")[GT_index]
    return iid_genotypes 

'''
def allele_ratio(offspring_AD):
    """
    Returns the allele ratio of the DNM and inherited allele
    The log2 coverage ratio (relative to median)
    """
    _buffer = 1 # Add one to avoid 0 values
    _ad = str(offspring_AD).split(":")
    pass
'''


def parse_variant(line,ped_dict,iid_indices):
    linesplit = tokenize(line)
    variant = tuple(map(str,linesplit[0:5]))
    chrom,pos,ref,alt = variant[0],variant[1],variant[3],variant[4]
    # Zero-based positions
    start,end = int(pos)-1,int(pos)+len(ref)-1
    # chrom:start:end:ref:alt
    ID = "{}:{}:{}:{}".format(chrom,pos,ref,alt)

    # Skip multiallelic variants
    n_alt = len(alt.split(","))
    if n_alt > 1: return

    # Skip X and Y variants for now
    if "X" in chrom or "Y" in chrom: return


    FORMAT = linesplit[8] # Genotype FORMAT field (contains GT, DP, etc.)
    format_indices = index_format(FORMAT)
    # According to doc, none of the fields in the FORMAT column are required 
    # We should require GT, GQ, PL, AD, DP?
    # iid_genotypes = load_genotypes(linesplit,iid_indices,format_indices)

    for offspring_iid in ped_dict:
        father_iid, mother_iid = ped_dict[offspring_iid][2],ped_dict[offspring_iid][3]
        if offspring_iid not in iid_indices or father_iid not in iid_indices or mother_iid not in iid_indices: continue
        offspring_GT = get_feature(linesplit[iid_indices[offspring_iid]],format_indices["GT"])
        offspring_AD = get_feature(linesplit[iid_indices[offspring_iid]],format_indices["AD"])
        offspring_DP = get_feature(linesplit[iid_indices[offspring_iid]],format_indices["DP"])
        offspring_GQ = get_feature(linesplit[iid_indices[offspring_iid]],format_indices["GQ"])
        offspring_PL = get_feature(linesplit[iid_indices[offspring_iid]],format_indices["PL"])

        # (Maybe later allow for 0/1,1/1,1/1 genotypes)
        offspring_alleles, father_alleles, mother_alleles = offspring_genotype.replace("|","/").split("/"), father_genotype.replace("|","/").split("/"), mother_genotype.replace("|","/").split("/")
        # Get genotype features  

        # Print putative DNM for that offspring


def parse(vcf_filepath, ped_filepath):
    from backend import process_ped
    ped_dict = process_ped(ped_filepath)

    if vcf_filepath.endswith(".gz"):
        import gzip
        f = gzip.open(vcf_filepath,"rt",9)
    else: f = open(vcf_filepath,"r")

    for line in f:
        if line.startswith("#CHROM"): iid_indices = index_samples(line, ped_dict) # Use the CHROM header line to index the samples
        if line.startswith("#"): continue

        parse_variant(line, ped_dict, iid_indices)
        # print(line.rstrip())
        
