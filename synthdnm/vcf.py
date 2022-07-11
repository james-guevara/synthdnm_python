from .backend import tokenize
import numpy as np
import sys


def index_samples(CHROM_header_line, ped_dictionary):
    """Processes CHROM header line and stores column indices with sample ID"""
    iid_indices = {}  # iid_indices[iid] = index
    linesplit = tokenize(CHROM_header_line)
    for i in range(9, len(linesplit)):
        iid = linesplit[i]  # The sample ID (IID) at the ith column
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
    """Gets feature value at format_index"""
    if format_index >= len(iid_feature_index.split(":")):
        return None
    return iid_feature_index.split(":")[format_index]


def get_info_features(info_keys, info_column):
    """Get info features for variant"""
    info_features = {
        key: None for key in info_keys
    }  # Initialize each key's value to None
    for field in info_column.split(";"):
        if "=" not in field:
            continue
        key, val = field.split("=")
        if key not in info_keys:
            continue
        else:
            info_features[key] = val
    return info_features


def parse_variant(
    line, ped_dict, iid_indices, info_keys, fout=None, training_examples=None
):
    """Parse line for variant"""
    linesplit = tokenize(line)
    variant = tuple(map(str, linesplit[0:5]))
    chrom, pos, ref, alt = variant[0], variant[1], variant[3], variant[4]
    start, end = int(pos) - 1, int(pos) + len(ref) - 1  # Zero-based positions
    ID = "{}:{}:{}:{}".format(chrom, start, end, ref, alt)  # chrom:start:end:ref:alt

    # Skip multiallelic variants
    n_alt = len(alt.split(","))
    if n_alt > 1:
        return

    # Skip X and Y variants for now
    if "X" in chrom or "Y" in chrom:
        return

    INFO = linesplit[7]
    info_features = get_info_features(info_keys, INFO)

    format_indices = index_format(
        FORMAT=linesplit[8]
    )  # Genotype FORMAT field (contains GT, DP, etc.)
    # According to doc, none of the fields in the FORMAT column are required
    # We should require GT, GQ, PL, AD, DP?
    required_format_fields = ["GT", "GQ", "PL", "AD", "DP"]
    if not all(field in format_indices for field in required_format_fields):
        print(
            "Variant is missing a required field; skipping [{}]".format(ID),
            file=sys.stderr,
        )
        return

    for offspring_iid in ped_dict:
        father_iid, mother_iid = ped_dict[offspring_iid][2], ped_dict[offspring_iid][3]
        if (
            offspring_iid not in iid_indices
            or father_iid not in iid_indices
            or mother_iid not in iid_indices
        ):
            continue
        offspring_GT, father_GT, mother_GT = (
            get_feature(linesplit[iid_indices[offspring_iid]], format_indices["GT"]),
            get_feature(linesplit[iid_indices[father_iid]], format_indices["GT"]),
            get_feature(linesplit[iid_indices[mother_iid]], format_indices["GT"]),
        )
        offspring_AD, father_AD, mother_AD = (
            get_feature(linesplit[iid_indices[offspring_iid]], format_indices["AD"]),
            get_feature(linesplit[iid_indices[father_iid]], format_indices["AD"]),
            get_feature(linesplit[iid_indices[mother_iid]], format_indices["AD"]),
        )
        offspring_DP, father_DP, mother_DP = (
            get_feature(linesplit[iid_indices[offspring_iid]], format_indices["DP"]),
            get_feature(linesplit[iid_indices[father_iid]], format_indices["DP"]),
            get_feature(linesplit[iid_indices[mother_iid]], format_indices["DP"]),
        )
        offspring_GQ, father_GQ, mother_GQ = (
            get_feature(linesplit[iid_indices[offspring_iid]], format_indices["GQ"]),
            get_feature(linesplit[iid_indices[father_iid]], format_indices["GQ"]),
            get_feature(linesplit[iid_indices[mother_iid]], format_indices["GQ"]),
        )
        offspring_PL, father_PL, mother_PL = (
            get_feature(linesplit[iid_indices[offspring_iid]], format_indices["PL"]),
            get_feature(linesplit[iid_indices[father_iid]], format_indices["PL"]),
            get_feature(linesplit[iid_indices[mother_iid]], format_indices["PL"]),
        )

        # (Maybe later allow for 0/1,1/1,1/1 genotypes)
        if not (
            offspring_GT.replace("|", "/") == "0/1"
            and father_GT.replace("|", "/") == "0/0"
            and mother_GT.replace("|", "/") == "0/0"
        ):
            continue
        if offspring_PL is None or father_PL is None or mother_PL is None:
            continue
        # Get genotype features
        offspring_PL, father_PL, mother_PL = (
            offspring_PL.split(","),
            father_PL.split(","),
            mother_PL.split(","),
        )
        (offspring_hom_ref_pl, offspring_het_pl, offspring_hom_alt_pl) = (
            offspring_PL[0],
            offspring_PL[1],
            offspring_PL[2],
        )
        (father_hom_ref_pl, father_het_pl, father_hom_alt_pl) = (
            father_PL[0],
            father_PL[1],
            father_PL[2],
        )
        (mother_hom_ref_pl, mother_het_pl, mother_hom_alt_pl) = (
            mother_PL[0],
            mother_PL[1],
            mother_PL[2],
        )

        (parent_min_DNM_PL, parent_max_DNM_PL) = sorted(
            np.array([father_het_pl, mother_het_pl], dtype=np.float64)
        )
        (parent_min_HOM_REF_PL, parent_max_HOM_REF_PL) = sorted(
            np.array([father_hom_ref_pl, mother_hom_ref_pl], dtype=np.float64)
        )

        offspring_ar = get_offspring_ar(offspring_AD)
        offspring_DP_, father_DP_, mother_DP_ = (
            get_log2_coverage_ratio(offspring_AD),
            get_log2_coverage_ratio(father_AD),
            get_log2_coverage_ratio(mother_AD),
        )
        (parent_min_DP, parent_max_DP) = sorted(
            np.array([father_DP_, mother_DP_], dtype=np.float64)
        )
        (parent_min_GQ, parent_max_GQ) = sorted(
            np.array([father_GQ, mother_GQ], dtype=np.float64)
        )

        filter_ = None
        qual = None

        parent_ar_max = 0
        parent_ar_min = 0

        # fout line with features
        feats = [
            offspring_GT,
            father_GT,
            mother_GT,
            n_alt,
            filter_,
            qual,
            parent_ar_max,
            parent_ar_min,
            offspring_ar,
            parent_max_DP,
            parent_min_DP,
            offspring_DP_,
            parent_max_DNM_PL,
            parent_min_DNM_PL,
            parent_max_HOM_REF_PL,
            parent_min_HOM_REF_PL,
            offspring_het_pl,
            offspring_hom_ref_pl,
            parent_max_GQ,
            parent_min_GQ,
            offspring_GQ,
        ] + list(info_features.values())

        out = (
            "{}\t{}\t{}\t{}\t{}\t{}".format(chrom, pos, ID, ref, alt, offspring_iid)
            + "\t"
            + "\t".join([str(elem) if elem is not None else "" for elem in feats])
        )
        if training_examples:
            out = out + "\t" + training_examples
        if fout:
            print(out, file=fout)
        else:
            print(out)


def get_offspring_ar(offspring_AD):
    # Add one to avoid 0 values in the denominator
    _buff = 1
    ref_AD, alt_AD = offspring_AD.split(",")[0], offspring_AD.split(",")[1]
    if ref_AD is None or alt_AD is None:
        return None
    offspring_AR = float(alt_AD) / (float(ref_AD) + _buff)
    return offspring_AR


def get_log2_coverage_ratio(AD):
    _buff = 1
    ref_AD, alt_AD = float(AD.split(",")[0]), float(AD.split(",")[1])
    if ref_AD is None or alt_AD is None:
        return None
    coverage_ratio = np.log2(
        (ref_AD + alt_AD + _buff) / (np.median([ref_AD, alt_AD]) + _buff)
    )
    return coverage_ratio


def parse(
    vcf_filepath,
    ped_filepath,
    info_keys,
    variants_to_keep=None,
    fout=None,
    training_examples=None,
):
    from .backend import process_ped

    ped_dict = process_ped(ped_filepath)

    # The INFO column depends on the VCF type
    header = "\t".join(
        ["chrom", "pos", "ID", "ref", "alt", "iid"]
        + [
            "offspring_gt",
            "father_gt",
            "mother_gt",
            "nalt",
            "filter",
            "qual",
            "parent_ar_max",
            "parent_ar_min",
            "offspring_ar",
            "parent_dp_max",
            "parent_dp_min",
            "offspring_dp",
            "parent_dnm_pl_max",
            "parent_dnm_pl_min",
            "parent_inh_pl_max",
            "parent_inh_pl_min",
            "offspring_dnm_pl",
            "offspring_inh_pl",
            "parent_gq_max",
            "parent_gq_min",
            "offspring_gq",
        ]
        + info_keys
    )
    if training_examples:
        header = header + "\t" + "TRUTH"
    if fout:
        if fout.tell() == 0:
            print(header, file=fout)
    else:
        print(header)

    if vcf_filepath.endswith(".gz"):
        import gzip

        f = gzip.open(vcf_filepath, "rt", 9)
    else:
        f = open(vcf_filepath, "r")

    for line in f:
        if line.startswith("#CHROM"):
            iid_indices = index_samples(
                line, ped_dict
            )  # Use the CHROM header line to index the samples
        if line.startswith("#"):
            continue
        parse_variant(
            line,
            ped_dict,
            iid_indices,
            info_keys,
            fout=fout,
            training_examples=training_examples,
        )


###
def parse_private_inherited_variant(line, ped_dict, iid_indices, fout):
    """Parse line for variant"""
    linesplit = tokenize(line)
    variant = tuple(map(str, linesplit[0:5]))
    chrom, pos, ref, alt = variant[0], variant[1], variant[3], variant[4]
    start, end = int(pos) - 1, int(pos) + len(ref) - 1  # Zero-based positions
    ID = "{}:{}:{}:{}".format(chrom, start, end, ref, alt)  # chrom:start:end:ref:alt
    format_indices = index_format(
        FORMAT=linesplit[8]
    )  # Genotype FORMAT field (contains GT, DP, etc.)
    for offspring_iid in ped_dict:
        father_iid, mother_iid = ped_dict[offspring_iid][2], ped_dict[offspring_iid][3]
        if (
            offspring_iid not in iid_indices
            or father_iid not in iid_indices
            or mother_iid not in iid_indices
        ):
            continue
        offspring_GT, father_GT, mother_GT = (
            get_feature(linesplit[iid_indices[offspring_iid]], format_indices["GT"]),
            get_feature(linesplit[iid_indices[father_iid]], format_indices["GT"]),
            get_feature(linesplit[iid_indices[mother_iid]], format_indices["GT"]),
        )
        # (Maybe later allow for 0/1,1/1,1/1 genotypes)
        if (
            offspring_GT.replace("|", "/") == "0/1"
            and father_GT.replace("|", "/") == "0/1"
            and mother_GT.replace("|", "/") == "0/0"
        ) or (
            offspring_GT.replace("|", "/") == "0/1"
            and father_GT.replace("|", "/") == "0/0"
            and mother_GT.replace("|", "/") == "0/1"
        ):
            # print("{}\t{}\t{}\t{}".format(ID,offspring_iid,father_iid,mother_iid))
            fout.write(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    ID,
                    offspring_iid,
                    father_iid,
                    mother_iid,
                    offspring_GT,
                    father_GT,
                    mother_GT,
                )
            )


def parse_private_inherited(vcf_filepath, ped_filepath):
    from .backend import process_ped

    ped_dict = process_ped(ped_filepath)
    if vcf_filepath.endswith(".gz"):
        import gzip

        f = gzip.open(vcf_filepath, "rt", 9)
    else:
        f = open(vcf_filepath, "r")

    fout = open("private_transmitted_IDs.txt", "w")
    for line in f:
        if line.startswith("#CHROM"):
            iid_indices = index_samples(
                line, ped_dict
            )  # Use the CHROM header line to index the samples
        if line.startswith("#"):
            continue
        parse_private_inherited_variant(line, ped_dict, iid_indices, fout)
