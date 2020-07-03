from __future__ import print_function

from cyvcf2 import VCF
import sys
import logging

logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__file__)

__version__ = "0.1"

args = sys.argv

vcf_file = str(args[1])
caller = str(args[2])
max_gnomad_af = float(args[3])
max_internal_af = 0.2
vep_impact_terms = ["HIGH", "MODERATE", "MODIFIER"]

LOG.info("running SV-filter version %s with parameters", __version__)
LOG.info("vcf-file: %s", vcf_file)
LOG.info("caller: %s", caller)
LOG.info("Max GnomAD AF: %s", max_gnomad_af)
LOG.info("Max within samples AF: %s", max_internal_af)
LOG.info("VEP IMPACT terms: %s", vep_impact_terms)



if caller == "cnvnator":
    min_len = 5000
elif caller == "delly":
    min_len = 50
else:
    LOG.info("Enter delly, or cnvnator as caller")
    sys.exit()


def CSQ_to_dict(csq_keys, csq_string):
    csq_dict = dict()
    csq_values = csq_string.split("|")
    for key, value in zip(csq_keys, csq_values):
        if value == "":
            csq_value = None
        else:
            csq_value = value
        csq_dict[key] = csq_value

    return csq_dict

def parse_genotype(vcf_record, samples):

    genotype_dict = dict()
    for sample, genotype in zip(samples, vcf_record.genotypes):
        sample_name = sample.split("/")[-1]
        if 0 in genotype[0:2] and 1 in genotype[0:2]:
            genotype_dict[sample_name] = "HET"

        elif 1 in genotype[0:2]:
            genotype_dict[sample_name] = "ALT_HOM"

        else:
            genotype_dict[sample_name] = "REF_HOM"

    return genotype_dict


def get_alt_samples(genotype_dict):

    het_list = []
    hom_list = []

    for key, value in genotype_dict.items():
        if value == "HET":
            het_list.append(key)
        if value == "ALT_HOM":
            hom_list.append(key)

    return het_list, hom_list


vcf = VCF(vcf_file)

samples = vcf.samples

tot_alleles = 2*len(samples)

CSQ_description = vcf.get_header_type(key="CSQ").get("Description")
CSQ_format = CSQ_description.split("Format: ")[-1]
CSQ_keys = CSQ_format.split("|")

print(vcf.raw_header, end="")

for record in vcf:

    svlen = record.INFO.get("SVLEN")
    if svlen is None:
        svlen = record.INFO["END"] - record.POS


    genotypes = parse_genotype(record, samples)
    het_list, hom_list = get_alt_samples(genotypes)

    num_alleles = len(het_list) + 2*len(hom_list)
    af = num_alleles / tot_alleles

    consequences = record.INFO.get("CSQ").split(",")
    first_csq = CSQ_to_dict(csq_keys=CSQ_keys, csq_string=consequences[0])

    if first_csq["gnomad_sv_AF"] not in [None, ""]:
        gnomad_af = float(first_csq["gnomad_sv_AF"].split("&")[0])
        if gnomad_af > max_gnomad_af:
            continue
    if svlen < min_len:
        continue
    if af > max_internal_af:
        continue

    for consequence in consequences:

        csq_dict = CSQ_to_dict(csq_keys=CSQ_keys, csq_string=consequence)

        if csq_dict["IMPACT"] in vep_impact_terms:
            print(record, end="")
            break
