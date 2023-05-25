### main.py is at the very end

###
### preprocessor.py
###
# from sample import generate_samples, find_variants, FamilyMember
# from filter import VariantFilter, filter_variants_without_log
# from variantFields import parse_variant_header
# from analysisDriver import AnalysisType, AffectedStatus
# reads and processes input files: gets VariantFields, samples to include in analysis, and proband's variants to start analysis
class Preprocessor:
    def __init__(self, analysis_type, pedigree, annotation_files, log):
        self.proband_start_vars = None
        self.log = log
        self.member_ids = [FamilyMember.proband.value] # initialize with only proband
        self.header_string = find_header_string(annotation_files[0])
        self.variant_fields = parse_variant_header(self.header_string)
        self.samples = self.make_samples(analysis_type, pedigree, annotation_files)
    def make_samples(self, analysis_type, ped_file, annotation_files):
        samples_all = generate_samples(ped_file, annotation_files)
        samples_included = choose_samples(samples_all, analysis_type)
        return self.add_variants(samples_included)
    # add each sample's variants to its Sample object
    # proband handled differently from other samples
    def add_variants(self, samples_empty):
        samples_vars = []
        for sample in samples_empty:
            if sample.family_member_id == FamilyMember.proband.value:
                with_vars = self.find_proband_vars(sample)
            else:
                with_vars = self.find_other_vars(sample)
            samples_vars.append(with_vars)
        return samples_vars
    # finds proband vars: adds steps to log, assigns class variable proband initial vars
    def find_proband_vars(self, proband):
        proband_all_vars = find_variants(proband)
        self.proband_start_vars = len(proband_all_vars)
        proband_filter = VariantFilter(proband_all_vars, self.variant_fields, self.log)
        proband_variants = proband_filter.filter()
        proband.set_variants(proband_variants)
        return proband
    # finds variants for non-proband samples: does not add steps to log
    def find_other_vars(self, sample):
        sample_all_vars = find_variants(sample)
        sample_variants = filter_variants_without_log(sample_all_vars, self.variant_fields)
        sample.set_variants(sample_variants)
        self.member_ids.append(sample.family_member_id)
        return sample
    def get_all_samples(self):
        return self.samples
    # returns int
    def get_proband_start_vars(self):
        return self.proband_start_vars
    def get_member_ids(self):
        return self.member_ids
    def get_header_str(self):
        return self.header_string
    def get_variant_fields(self):
        return self.variant_fields
    ############ END CLASS METHODS ############
def choose_samples(samples_all, analysis_type):
    samples_selected = [samples_all[0]]    # always include proband
    # AD_NM excludes all except proband and parent(s)
    if analysis_type == AnalysisType.AD_NM:
        for sample in samples_all[1:]:    # iterate over non-proband samples
            i = sample.family_member_id
            if i == FamilyMember.father.value or i == FamilyMember.mother.value:
                samples_selected.append(sample)
        return samples_selected
    # AR_H excludes affected parents
    elif analysis_type == AnalysisType.AR_H:
        for sample in samples_all[1:]:    # iterate over non-proband samples
            i = sample.family_member_id
            if i == FamilyMember.mother.value or i == FamilyMember.father.value:
                if sample.get_affected_status() != AffectedStatus.Affected.value:
                    samples_selected.append(sample)
            else:
                samples_selected.append(sample)
        return samples_selected
    else:
        return samples_all
# def find_header_string(annotation_file):
#     header_string = ''
#     with open(annotation_file, 'r') as f:
#         for line in f.readlines():
#             if line.startswith('Chromosome') or line.startswith('Chr'):
#                 header_string = line
#                 break
#     if header_string == '':
#         raise Exception('Header not found in file:', annotation_file)
#     return header_string
def find_header_string(annotation_file):
    header_string = ''
    with open(annotation_file, 'r') as f:
        for line in f.readlines():
            if line.startswith('#CHROM') or line.startswith('Chr'):
                header_string = line
                break
    if header_string == '':
        raise Exception('Header not found in file:', annotation_file)
    return header_string




###
### analysisDriver.py
###
# from compoundAnalyzer import *
from enum import Enum
# from analyzer import Analyzer, FORMAT_START, GENOTYPE_HOM, GENOTYPE_HET
# from sample import FamilyMember
# AFFECTED STATUS ACCORDING TO GATK PED
class AffectedStatus(Enum):
    Affected = 2
    Unaffected = 1
    Missing = 0
# ANALYSIS_TYPE_AR_CH: ALL VARIANTS APPEARING AT LEAST TWICE IN A GENE IN THE SAME PATTERN
#     IN ALL AFFECTED SAMPLES, AND NOT FOUND IN THAT PATTERN IN ANY UNAFFECTED SAMPLE
# ANALYSIS_TYPE_AR_H: ALL VARIANTS FOUND ON ALL AFFECTED INDIVIDUALS AS HOMOZYGOUS AND NOT FOUND ON UNAFFECTED
#     INDIVIDUALS AS HOMOZYGOUS
# ANALYSIS_TYPE_AD_NM: ALL VARIANTS FOUND ON PROBAND NOT PRESENT ON PARENTS (ONLY APPLIED FOR TRIO ANALYSIS)
# ANALYSIS_TYPE_AD_IM: ALL VARIANTS FOUND ON AFFECTED SAMPLES AND NOT FOUND ON UNAFFECTED SAMPLES (APPLIED FOR
#     MULTIPLES AFFECTED INDIVIDUALS IN THE FAMILY)
# ANALYSIS_TYPE_AD_V: ALL HETEROZYGOUS VARIANTS PRESENT ON PROBAND
class AnalysisType(Enum):
    AR_CH = 0
    AR_H = 1
    AD_NM = 2
    AD_IM = 3
    AD_V = 4
class Analysis(object):
    def __init__(self, variants_fields, log_text, all_samples):
        self.format_var = FORMAT_START
        self.variants_fields = variants_fields
        self.analyzer = Analyzer(variants_fields, self.format_var)
        self.log_text = log_text
        self.all_samples = all_samples
        self.affected = []  # list of affected samples, including proband
        self.unaffected = []  # list of unaffected samples
        self.get_samples_status()  # populates the above two samples lists
        self.proband_vars = self.all_samples[0].get_variants()  # proband is always first
    def analyze_variants(self, analysis_type):
        if analysis_type is AnalysisType.AD_V:
            return self.run_AD_V()
        elif analysis_type == AnalysisType.AD_IM:
            return self.run_AD_IM()
        elif analysis_type == AnalysisType.AD_NM:
            return self.run_AD_NM()
        elif analysis_type == AnalysisType.AR_H:
            return self.run_AR_H()
        elif analysis_type == AnalysisType.AR_CH:
            return self.run_AR_CH()
        else:
            raise Exception('Invalid analysis type:', analysis_type)
    #
    # ---------------- 1. AD_V ----------------
    # only uses proband
    # returns all proband's heterozygous variants
    def run_AD_V(self):
        proband_het_vars = self.analyzer.get_heterozygous_vars(self.proband_vars, self.log_text)
        return proband_het_vars
    #
    # ---------------- 2. AD_IM ----------------
    # requires proband and at least one relative, who can be either affected or not
    # returns variants shared by only affected people
    def run_AD_IM(self):
        proband_het_vars = self.analyzer.get_heterozygous_vars(self.proband_vars, self.log_text)
        inherited_vars = []
        # -------- 2.a. get variants in both proband and other affected samples --------
        for var in proband_het_vars:
            proband_var = self.analyzer.format_variant(var)
            inherited = True
            for sample in self.affected[1:]:  # start at 1 to skip proband
                other_vars = sample.get_variants()
                # every proband variant must be in every affected sample's variants to be appended
                if not self.analyzer.match_variant(proband_var, other_vars):
                    inherited = False
                    break
            if inherited:
                inherited_vars.append(var)
        self.log_text.append_text(
            'Including if start position is in all affected members with het Genotype {0} variants -> {1} variants'.format(
                str(len(proband_het_vars)), str(len(inherited_vars))))
        # -------- 2.b. exclude variants in an unaffected sample --------
        return self.analyzer.exclude_unaffected(inherited_vars, self.unaffected, self.log_text)
    #
    # ---------------- 3. AD_NM ----------------
    # requires unaffected parents, and at least 1 parent sample
    # returns variants (het and hom) that only the proband has
    # all samples but proband and parent(s) excluded by Preprocessor
    def run_AD_NM(self):
        proband_het_hom = []
        # -------- 3.a. get all proband's heterozygous and homozygous variants
        for var in self.proband_vars:
            genotype = var.split("\t")[self.variants_fields.get_genotype()]
            if genotype == GENOTYPE_HET or genotype == GENOTYPE_HOM:
                proband_het_hom.append(var)
        self.log_text.append_text('Including {0}, {1} in Genotype, {2} variants -> {3} variants'.format(
            GENOTYPE_HET, GENOTYPE_HOM, str(len(self.proband_vars)), str(len(proband_het_hom))))
        # -------- 3.b. exclude variants in unaffected parent samples --------
        # the only others in all_samples are unaffected parent(s)
        result_variants = proband_het_hom.copy()
        for sample in self.all_samples[1:]:  # start at 1 to skip proband
            start_count = len(result_variants)
            other_vars = sample.get_variants()
            for var in proband_het_hom:
                variant = self.analyzer.format_variant(var)
                if self.analyzer.match_variant(variant, other_vars):
                    try:
                        result_variants.remove(var)
                    except ValueError:  # do not show error if already removed
                        pass
            parent = FamilyMember(sample.get_family_member_id()).name
            self.log_text.append_text(
                'Excluding inherited from {0} {1} variants -> {2} variants'.format(
                    parent, str(start_count), str(len(result_variants))))
        return result_variants
    #
    # ---------------- 4. AR_H ----------------
    # requires proband only
    # returns variants for which only affected non-parents are homozygous
    # affected parent samples were removed by Preprocessor
    def run_AR_H(self):
        # -------- 4.a. start with all proband's homozygous variants --------
        proband_hom_vars = self.analyzer.get_homozygous_vars(self.proband_vars, True, self.log_text)
        # -------- 4.b. keep only variants homozygous in all other affected samples --------
        vars_hom_affected = []
        inherit_affected = []
        for sample in self.affected[1:]:  # start at 1 to skip proband
            vars_affected = sample.get_variants()
            vars_hom_affected.append(self.analyzer.get_homozygous_vars(vars_affected, False, self.log_text))
        for var in proband_hom_vars:
            variant = self.analyzer.format_variant(var)
            inherited = True
            for var_affected in vars_hom_affected:
                if not self.analyzer.match_variant(variant, var_affected):
                    inherited = False
                    break
            if inherited:
                inherit_affected.append(var)
        self.log_text.append_text(
            'Including if start position is in all affected members with hom Genotype, {0} variants -> {1} variants'.format(
                str(len(proband_hom_vars)), str(len(inherit_affected))))
        # -------- 4.c. exclude variants homozygous in any unaffected sample --------
        vars_hom_unaffected = []
        result_variants = []
        for sample in self.unaffected:
            vars_unaffected = sample.get_variants()
            vars_hom_unaffected.append(self.analyzer.get_homozygous_vars(vars_unaffected, False, self.log_text))
        for var in inherit_affected:
            variant = self.analyzer.format_variant(var)
            inherited = True
            for var_unaffected in vars_hom_unaffected:
                if self.analyzer.match_variant(variant, var_unaffected):
                    inherited = False
                    break
            if inherited:
                result_variants.append(var)
        self.log_text.append_text(
            'Excluding if start position is in any unaffected members with hom Genotype, '
            '{0} variants -> {1} variants'.format(
                str(len(inherit_affected)), str(len(result_variants))))
        return result_variants
    #
    # ---------------- 5. AR_CH ----------------
    # Requires proband only
    # Returns combination of >= 2 variants for a gene which only affected individuals have
    def run_AR_CH(self):
        # -------- 5.a. exclude variants in any unaffected sample --------
        # proband_im = self.analyze_variants(AnalysisType.AD_IM, log_text, affected)
        vars_affected_only = self.analyzer.exclude_unaffected(self.proband_vars, self.unaffected, self.log_text)
        # -------- 5.b. exclude genes with only a single variant --------
        c_analyzer = CompoundAnalyzer(self.analyzer)
        genes_all = c_analyzer.split_vars_into_genes(vars_affected_only)
        genes = get_genes_multiple_vars(genes_all)
        self.log_text.append_text(
            'Excluding variants in a gene with single variants, {0} variants -> {1} variants'.format(
                str(len(vars_affected_only)), str(count_variants_genes(genes))))
        # -------- 5.c. if two parent samples: exclude variants from only one parent --------
        two_parent_samples = check_two_parents(self.all_samples)
        inherited_genes = []
        if not two_parent_samples:
            inherited_genes = genes
        else:
            parents_genes = []
            for sample in self.all_samples:
                family_member = sample.get_family_member_id()
                if family_member == FamilyMember.father.value or family_member == FamilyMember.mother.value:
                    split_vars = c_analyzer.split_vars_into_genes(sample.get_variants())
                    parents_genes.append(split_vars)
            for gene in genes:
                find_in_single_parent = False
                for parent_genes in parents_genes:
                    if c_analyzer.match_genes_ch(gene, parent_genes):
                        find_in_single_parent = True
                        break
                if not find_in_single_parent:
                    inherited_genes.append(gene)
            self.log_text.append_text(
                'Excluding variants with single parent inheritance, {0} variants -> {1} variants'.format(
                    str(count_variants_genes(genes)),
                    str(count_variants_genes(inherited_genes))))
    # -------- 5.d. exclude genes in unaffected samples --------
        genes_unaffected_variants = []
        analyzed_genes = []
        result_variants = []
        for sample in self.unaffected:
            vars_unaffected = sample.get_variants()
            genes_unaffected_variants.append(c_analyzer.split_vars_into_genes(vars_unaffected))
        for inherited_gene in inherited_genes:
            inherited = True
            for genes_unaffected in genes_unaffected_variants:
                if c_analyzer.match_genes_ch(inherited_gene, genes_unaffected):
                    inherited = False
                    break
            if inherited:
                analyzed_genes.append(inherited_gene)
        for gene in analyzed_genes:
            gene_vars = gene.get_variants()
            for var in gene_vars:
                result_variants.append(var)
        self.log_text.append_text(
            'Excluding gene if all start positions are in any unaffected members for each gene, '
            '{0} variants -> {1} variants'.format(
                str(count_variants_genes(inherited_genes)), str(len(result_variants))))
        return result_variants
    def get_samples_status(self):
        for sample in self.all_samples:
            if sample.get_affected_status() == AffectedStatus.Affected.value:
                self.affected.append(sample)
            else:
                self.unaffected.append(sample)




###
### log.py
###
from datetime import datetime
from pytz import timezone
class Log:
    def __init__(self):
        self.text = ""  # single string for analysis descriptions, appended throughout app
        self.time = get_time()
    def append_text(self, text):
        self.text = self.text + text + "\n"
    # write analysis descriptions to file
    def write_log(self, analysis_type):
        out_name = 'Log_' + analysis_type.name + '_' + self.time + '.txt'
        with open(out_name, 'w') as f:
            f.write(self.text + '\n\n')
    # write final variants to csv file
    def write_variants(self, analysis_type, header, variants):
        out_name = 'PhenoDB_Analysis_' + analysis_type.name + '_' + self.time + '.tsv'
        with open(out_name, 'w') as f:
            f.write(header + '\n')
            for var in variants:
                f.write(var + '\n')
    def print_summary(self, analysis_type):
        print('Summary of analysis', analysis_type, '\n')
        print(self.text + '\n\n')
def get_time():
    est = datetime.now(timezone('US/Eastern'))
    return est.strftime('%Y_%m_%d_%H-%M-%S')




###
### sample.py
###
from enum import Enum
# MEMBER_ID ACCORDING TO PHENODB
# REMAINING SAMPLES SHOULD BE NUMBERED FROM 4 TO N, REGARDLESS THE PRESENCE OF MOTHER OR FATHER
class FamilyMember(Enum):
    proband = 1
    mother = 2
    father = 3
# -----------------------------------------------------------------------------
#   Class:          Sample
#   Description:    Class to collect information from PED file and header from ANNOVAR file
#   Methods:        getters(self): returns the content of the field
#                   setters(self, param): set the content of the field according to the parameter
class Sample:
    def __init__(self, family_member_id, family_id, sample_id, paternal_id, maternal_id, sex, affected_status,
                 file_path):
        self.family_member_id = family_member_id
        self.family_id = family_id
        self.sample_id = sample_id
        self.paternal_id = paternal_id
        self.maternal_id = maternal_id
        self.sex = sex
        self.affected_status = affected_status
        self.variants = []
        self.file_path = file_path
    def get_variants(self):
        return self.variants
    def get_family_member_id(self):
        return self.family_member_id
    def get_family_id(self):
        return self.family_id
    def get_sample_id(self):
        return self.sample_id
    def get_paternal_id(self):
        return self.paternal_id
    def get_maternal_id(self):
        return self.maternal_id
    def get_sex(self):
        return self.sex
    def get_affected_status(self):
        return self.affected_status
    def get_file_path(self):
        return self.file_path
    def set_variants(self, variants):
        self.variants = variants
    def __str__(self):
        return f'Sample({self.family_member_id},{self.family_id},{self.sample_id},{self.paternal_id},' \
               f',{self.maternal_id}{self.sex},{self.affected_status},{self.file_path},{len(self.variants)})'
# -----------------------------------------------------------------------------
#   Method:         read_sample(sample):
#   Description:    read variants on ANNOVAR file from sample
#   Parameters:     sample              sample object
#   Return:         variants            list of variants
def find_variants(sample):
    try:
        with open(sample.get_file_path(), 'r') as reader:
            variants = []
            for line in reader.readlines():
                if line[0] != "#" and not line.startswith('Chromosome'):
                    variants.append(line[:-1])
            return variants
    except:
        raise Exception('File not found:', sample.get_file_path())
# -----------------------------------------------------------------------------
#   Method:         read_ped_file(file_ped, file_path):
#   Description:    create a list of Sample() objects from PED file and ANNOVAR paths
#   Parameters:     file_ped            ped file containing information of samples
#                   file_path[]         list of files containing the ANNOVAR variants
#   Return:         sample[]            list of Sample() objects
# TODO relatedness temporarily determined from ped file, prod setup will be different
# in ped file: first line must be proband
# second line is either mother or word 'missing', third line either father or word 'missing'
def generate_samples(ped_file, annotation_files):
    with open(ped_file, 'r') as reader:
        family_member_id_count = 0
        sample_count = 0
        samples = []
        for line in reader.readlines():
            family_member_id_count += 1    # proband is first line, and family member id 1
            # do not create sample for family ids 2 or 3 (mother, father) if sample unavailable
            if line.startswith('missing'):
                continue
            family_id = line.split("\t")[0]
            sample_id = line.split("\t")[1]
            paternal_id = line.split("\t")[2]
            maternal_id = line.split("\t")[3]
            sex = int(line.split("\t")[4])
            affected_status = int(line.split("\t")[5].replace("\n", ""))
            samples.append(
                Sample(family_member_id_count, family_id, sample_id, paternal_id, maternal_id, sex, affected_status,
                       annotation_files[sample_count]))
            sample_count += 1
    return samples




###
### filter.py
###
# from log import Log
# THIS IS THE DEFAULT SELECTION THAT WAS IMPLEMENTED FOR ALL ANALYSIS
# PHENODB ALLOWS USERS TO FREELY SELECT THE DESIRED COMBINATION OF OPTIONS
# REFGENE_GENE_LOCATION_DEFAULT_SELECTION = ("exonic", "exonic;splicing", "splicing")
REFGENE_GENE_LOCATION_DEFAULT_SELECTION = ( # from https://useast.ensembl.org/info/genome/variation/prediction/predicted_data.html
    "transcript_ablation",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "stop_gained",
    "frameshift_variant",
    "stop_lost",
    "start_lost",
    "transcript_amplification",
    "inframe_insertion",
    "inframe_deletion",
    "missense_variant",
    "protein_altering_variant",
    "splice_region_variant",
    "incomplete_terminal_codon_variant",
    "start_retained_variant",
    "stop_retained_variant",
    "coding_sequence_variant")
# THIS IS THE DEFAULT MAF CUTOFF, SHOULD BE CHANGED ACCORDING TO THE NEEDS
MAF_CUTOFF = 0.01
# exclude some variants from analysis
class VariantFilter(object):
    def __init__(self, variants, columns, log_text):
        self.variants = variants
        self.columns = columns
        self.maf_cols = columns.get_cols_maf()  # col numbers with MAF from ANNOVAR file
        self.log_text = log_text
    def filter(self):
        temp0 = self.filter_refgene_gene_location(self.variants)
        # temp1 = self.filter_refgene_exon_function(temp0)
        # temp2 = self.filter_maf(temp1)
        temp2 = self.filter_maf(temp0)
        return temp2
    # -----------------------------------------------------------------------------
    #   Method:         filter_refgene_gene_location(variants):
    #   Description:    filter variants based on refGene gene location according to selection
    #   Parameters:     variants        list of variants
    #   Return:         variants_tmp    list of variants filtered variants
    def filter_refgene_gene_location(self, vs):
        variants_tmp = []
        # SELECT THE DESIRED REFGENE GENE LOCATIONS FROM THE OPTIONS DESCRIBED ON
        # REFGENE_GENE_LOCATION_OPTIONS, STANDARD PHENODB SELECTION IMPLEMENTED
        refgene_gene_location_selection = REFGENE_GENE_LOCATION_DEFAULT_SELECTION
        for var in vs:
            if self.check_refgene_gene_location(var, refgene_gene_location_selection):
                variants_tmp.append(var)
        self.log_text.append_text(
            'Including {0} in RefGeneLocation, {1} variants -> {2} variants'.format(refgene_gene_location_selection,
                                                                                    str(len(vs)),
                                                                                    str(len(variants_tmp))))
        return variants_tmp
    # -----------------------------------------------------------------------------
    #   Method:         check_refgene_gene_location(line, refgene_location_selection):
    #   Description:    check if a variant is at desired location on refGene gene location
    #   Parameters:     line                            ANNOVAR line with a variant
    #                   refgene_location_selection      selection of desired locations
    #   Return:         True              qualified variant
    #                   False             not a qualified variant
    def check_refgene_gene_location(self, line, refgene_location_selection):
        refgene_gene_location = line.split("\t")[self.columns.get_gene_location()]
        for location in refgene_location_selection:
            # if refgene_gene_location == location:
            if location in refgene_gene_location:
                return True
        return False
    # -----------------------------------------------------------------------------
    #   Method:         filter_refgene_exon_function(variants):
    #   Description:    filter synonymous variants except for exonic:splicing
    #   Parameters:     variants        list of variants
    #   Return:         variants_tmp    list of variants filtered variants
    # def filter_refgene_exon_function(self, vs):
    #     variants_tmp = []
    #     for var in vs:
    #         if self.check_refgene_exon_function(var):
    #             variants_tmp.append(var)
    #     self.log_text.append_text(
    #         'Excluding synonymous SNV in RefGeneExonFunction, except for exonic:splicing variants {0} variants -> {1} variants'.format(
    #             str(len(vs)), str(len(variants_tmp))))
    #     return variants_tmp
    # -----------------------------------------------------------------------------
    #   Method:         check_refgene_exon_function(line):
    #   Description:    check if a synonymous is a qualified variant (synonymous & exonic:splicing)
    #   Parameters:     line              ANNOVAR line with a variant
    #   Return:         True              qualified variant
    #                   False             not a qualified variant
    # def check_refgene_exon_function(self, line):
    #     refgene_gene_location = line.split("\t")[self.columns.get_gene_location()]
    #     refgene_exon_function = line.split("\t")[self.columns.get_exon_func()]
    #     if refgene_exon_function == "synonymous SNV":
    #         if refgene_gene_location == "exonic;splicing":
    #             return True
    #         else:
    #             return False
    #     else:
    #         return True
    # filters variants based on maf
    def filter_maf(self, vs):
        variants_tmp = []
        for var in vs:
            print_var = True
            for n in self.maf_cols:
                if maf_format(var.split("\t")[n]) > MAF_CUTOFF:
                    print_var = False
                    break
            if print_var:
                variants_tmp.append(var)
        self.log_text.append_text(
            'Excluding if value is greater than, {0} in selected MAF projects, {1} variants -> {2} variants'.format(
                MAF_CUTOFF, str(len(vs)), str(len(variants_tmp))))
        return variants_tmp
# formats MAF to float for filtering comparison
def maf_format(m):
    if m == ".":
        return 0.0
    elif len(m) == 0:
        return 0.0
    elif m == "\t":
        return 0.0
    elif m == "\n":
        return 0.0
    elif m == ".\n":
        return 0.0
    else:
        return float(m)
# filters variants but without adding to the log used by the rest of the analysis
def filter_variants_without_log(sample_variants, columns):
    log = Log()
    filter_variants = VariantFilter(sample_variants, columns, log)
    variants = filter_variants.filter()
    return variants




###
### variantFields.py
###
# fields of ANNOVAR files to be used for analysis
# getters: return the index number of the field
class VariantFields(object):
    # def __init__(self, chromosome, start_position, end_position, ref_allele, alt_allele, genotype,
    #              gene_location, gene_name, exon_func, cols_maf):
    def __init__(self, chromosome, start_position, ref_allele, alt_allele, genotype,
                 gene_location, gene_name, cols_maf):
        self.chromosome = chromosome
        self.start_position = start_position
        # self.end_position = end_position
        self.end_position = start_position + length(ref_allele) - 1
        self.ref_allele = ref_allele
        self.alt_allele = alt_allele
        self.genotype = genotype
        self.gene_location = gene_location
        self.gene_name = gene_name
        # self.exon_func = exon_func
        self.cols_maf = cols_maf
    def get_chromosome(self):
        return self.chromosome
    def get_start_position(self):
        return self.start_position
    def get_end_position(self):
        return self.end_position
    def get_ref_allele(self):
        return self.ref_allele
    def get_alt_allele(self):
        return self.alt_allele
    def get_genotype(self):
        return self.genotype
    def get_gene_location(self):
        return self.gene_location
    def get_gene_name(self):
        return self.gene_name
    # def get_exon_func(self):
    #     return self.exon_func
    def get_cols_maf(self):
        return self.cols_maf
    def __str__(self):
        # return f'({self.chromosome},{self.start_position},{self.end_position},{self.ref_allele},' \
        #         f'{self.alt_allele},{self.genotype},{self.gene_location},{self.gene_name},{self.exon_func})'
        return f'({self.chromosome},{self.start_position},{self.end_position},{self.ref_allele},' \
               f'{self.alt_allele},{self.genotype},{self.gene_location},{self.gene_name})'
# read the tab-separated header of an ANNOVAR file to find the index numbers of its fields
# def parse_variant_header(header):
#     fields = header.split("\t")
#     n = 0
#     cols_maf = []
#     for field in fields:
#         if field == "Chromosome" or field == "Chr":
#             chromosome = n
#         elif field == "StartPosition" or field == "Start":
#             start_position = n
#         elif field == "EndPosition" or field == "End":
#             end_position = n
#         elif field == "ReferenceAllele" or field == "Ref":
#             ref_allele = n
#         elif field == "AlternativeAllele" or field == "Alt":
#             alt_allele = n
#         # elif field == "Genotype" or field == "GT":
#         elif field == "Otherinfo1":
#             genotype = n
#         # elif field == "RefgeneGeneLocation":
#         elif field == "Func.refGene":
#             gene_location = n
#         # elif field == "RefgeneGeneName":
#         elif field == "Gene.refGene":
#             gene_name = n
#         # elif field == "RefgeneExonFunction":
#         elif field == "ExonicFunc.refGene":
#             exon_func = n
#         elif field[:5] == "ExAC_" or field[:12] == "esp6500siv2_" \
#                 or field[:13] == "1000g2014oct_" or field[:13] == "1000g2015aug_":
#             cols_maf.append(n)
#         n = n + 1
#     return VariantFields(chromosome, start_position, end_position, ref_allele, alt_allele, genotype,
#                    gene_location, gene_name, exon_func, cols_maf)
def parse_variant_header(header):
    header = header.replace("\n", "")
    fields = header.split("\t")
    n = 0
    cols_maf = []
    for field in fields:
        if field == "#CHROM" or field == "#CHR":
            chromosome = n
        elif field == "POS":
            start_position = n
        # elif field == "EndPosition" or field == "End":
        #     end_position = n
        #     there is no end_position in vcf converted tsv files
        elif field == "REF":
            ref_allele = n
        elif field == "ALT":
            alt_allele = n
        elif field.endswith("_AdjustedGenotype"):
            genotype = n
        # elif field == "RefgeneGeneLocation":
        elif field == "CSQ_Consequence":
            gene_location = n
        # elif field == "RefgeneGeneName":
        elif field == "CSQ_SYMBOL":
            gene_name = n
        # elif field == "RefgeneExonFunction":
        # elif field == "ExonicFunc.refGene":
        #     exon_func = n
        #     there is no exon_func in vcf converted tsv files
        elif field[:20] == "INFO_gnomad_3_1_1_AF" or field[:13] == "CSQ_gnomAD_AF":
            cols_maf.append(n)
        n = n + 1
    return VariantFields(chromosome, start_position, ref_allele, alt_allele, genotype,
                   gene_location, gene_name, cols_maf)




###
### compoundAnalyzer.py
###
# from sample import FamilyMember
class CompoundAnalyzer(object):
    def __init__(self, analyzer):
        self.variant_fields = analyzer.get_variant_fields()
        self.format_var = analyzer.get_format_var()
        self.analyzer = analyzer
    # -----------------------------------------------------------------------------
    #   Method:         match_genes_ch(gene, match_genes):
    #   Description:    Check whether all variants in a given gene match to all variants in a given list of genes
    #   Parameters:     gene                Gene()
    #                   match_genes[]       list of Gene() to match
    #   Return:         True                If all variants in given gene matched to at all variants in list of genes
    #                   False               If all variants in given gene did not match to all variants in list of genes
    def match_genes_ch(self, gene, match_genes):
        for match_gene in match_genes:
            if gene.get_name() == match_gene.get_name():
                # print("Proband Gene " + gene.get_name() + " Parent Gene " + match_gene.get_name())
                if self.match_gene_variants(gene.get_variants(), match_gene.get_variants()):
                    return True
                else:
                    return False
        return False
    # -----------------------------------------------------------------------------
    #   Method:         match_gene_variants(variants, match_variants):
    #   Description:    remove gene with only one variant
    #   Parameters:     variants[]          list of variants
    #                   match_variants[]    list of variants to match
    #   Return:         True                If all variants matched between two sets of variants
    #   Return:         False               If all variants did not match between two sets of variants
    def match_gene_variants(self, variants, match_variants):
        match = True
        for var in variants:
            variant = self.analyzer.format_variant(var)
            if not self.analyzer.match_variant(variant, match_variants):
                match = False
        return match
    # -----------------------------------------------------------------------------
    #   Method:         split_vars_into_genes(variants):
    #   Description:    Split a list of variants into a list of genes with the given variants
    #   Parameters:     variants[]          list of ANNOVAR lines to be split into genes
    #   Return:         genes[]           list of genes containing ANNOVAR lines specific for each gene
    def split_vars_into_genes(self, variants):
        genes = []
        for variant in variants:
            gene_name = variant.split("\t")[self.variant_fields.get_gene_name()]
            index = find_gene(genes, gene_name)
            if index == -1:
                index = len(genes)
                genes.append(Gene(gene_name))
            genes.__getitem__(index).add_variant(variant)
        return genes
# -----------------------------------------------------------------------------
#   Method:         count_variants_genes(genes):
#   Description:    count number of variants in a list of genes
#   Parameters:     genes[]     list of genes
#   Return:         int         count of variants
def count_variants_genes(genes):
    analyzed_variants = []
    for gene in genes:
        variants = gene.get_variants()
        for var in variants:
            analyzed_variants.append(var)
    return len(analyzed_variants)
# -----------------------------------------------------------------------------
#   Method:         clear_genes_ch(genes):
#   Description:    remove gene with only one variant
#   Parameters:     genes[]       list of genes with >1 variant
def get_genes_multiple_vars(genes):
    genes_tmp = []
    for index in range(0, len(genes)):
        if genes[index].get_count() > 1:
            genes_tmp.append(genes[index])
    return genes_tmp
# returns True if there are two parent samples
def check_two_parents(all_samples):
    parents = 0
    for sample in all_samples:
        family_member = sample.get_family_member_id()
        if family_member == FamilyMember.father.value or family_member == FamilyMember.mother.value:
            parents += 1
    return parents == 2
# -----------------------------------------------------------------------------
#   Method:         find_gene(genes, gene_name):
#   Description:    Find a gene in a list and return its index
#   Parameters:     genes[]       list of genes
#                   gene_name     name of the gene to search on list of genes
#   Return:         int: -1 if not found, or the index if found
def find_gene(genes, gene_name):
    for index in range(0, len(genes)):
        if genes[index].get_name() == gene_name:
            return index
    return -1
# -----------------------------------------------------------------------------
#   Class:          Gene
#   Description:    Class to collect information regarding all variants in each gene for filtering steps
#   Methods:        getters(self): returns the content of the field
#                   setters(self, param): set the content of the field according to the parameter
class Gene:
    def __init__(self, name):
        self.name = name
        self.variants = []
        self.count = 0
    def get_name(self):
        return self.name
    def get_variants(self):
        return self.variants
    def get_count(self):
        return self.count
    def name(self, name):
        self.name = name
    def set_variants(self, variants):
        self.variants = variants
    def add_variant(self, variant):
        self.variants.append(variant)
        self.count = self.count + 1




###
### analyzer.py
###
# OUTPUT FORMAT FOR VARIANTS
# FORMAT_FULL {CHR}:{START}-{END}{REF}>{ALT}
# FORMAT_START {CHR}:{START}
FORMAT_FULL = 1
FORMAT_START = 2
# GENOTYPE ANNOTATION ACCORDING TO ANNOVAR
# GENOTYPE_HET = 'het'
# GENOTYPE_HOM = 'hom'
# GENOTYPE FROM VCF
GENOTYPE_HET = '0/1'
GENOTYPE_HOM = '1/1'
class Analyzer(object):
    def __init__(self, variant_fields, format_var):
        self.variant_fields = variant_fields
        self.format_var = format_var
    # -----------------------------------------------------------------------------
    #   Method:         match_variant(variant, sample_variants):
    #   Description:    match a given variant to a list of variants
    #   Parameters:     variant               formatted proband variant to match
    #                   sample_variants[]     list of variants to be matched
    #   Return:         True if variant match and
    #                   False if variant does not match
    def match_variant(self, variant, sample_variants):
        for var in sample_variants:
            matching_var = self.format_variant(var)
            if variant == matching_var:
                return True
        return False
    # -----------------------------------------------------------------------------
    #   Method:         format_variant(var, out_format):
    #   Description:    format variant as a single String
    #   Parameters:     var                 variant
    #                   format              output format full (1:10001-10001A>C) start (1:10001)
    #   Return:         String with variant in the chosen format
    def format_variant(self, var):
        try:
            chromosome = var.split("\t")[self.variant_fields.get_chromosome()]
            start = var.split("\t")[self.variant_fields.get_start_position()]
            end = var.split("\t")[self.variant_fields.get_end_position()]
            ref = var.split("\t")[self.variant_fields.get_ref_allele()]
            alt = var.split("\t")[self.variant_fields.get_alt_allele()]
            chr_num = chromosome
            if chromosome.capitalize().startswith("CHR"):
                chr_num = chromosome.capitalize().replace("CHR", "")
            formatted_variant = chr_num + ":" + start
            if self.format_var == FORMAT_FULL:
                formatted_variant = chr_num + ":" + start + "-" + end + ref + ">" + alt
            return formatted_variant
        except:
            raise Exception('Unable to format variant:', var)
    def parse_genotype(self, var):
        try:
            return var.split("\t")[self.variant_fields.get_genotype()]
        except:
            raise Exception('Unable to parse variant to get genotype:', var)
    def get_heterozygous_vars(self, variants, log_text):
        analyzed = []
        for var in variants:
            if self.parse_genotype(var) == GENOTYPE_HET:
                analyzed.append(var)
        log_text.append_text(
            'Including {0} in Genotype {1} variants -> {2} variants'.format(
                GENOTYPE_HET, str(len(variants)), str(len(analyzed))))
        return analyzed
    def get_homozygous_vars(self, variants, append_log, log_text):
        hom_variants = []
        for var in variants:
            if self.parse_genotype(var) == GENOTYPE_HOM:
                hom_variants.append(var)
        if append_log:
            log_text.append_text(
                'Including {0} in Genotype {1} variants -> {2} variants'.format(
                    GENOTYPE_HOM, str(len(variants)), str(len(hom_variants))))
        return hom_variants
    def exclude_unaffected(self, start_vars, unaffected, log_text):
        result_variants = []
        for var in start_vars:
            variant = self.format_variant(var)
            inherited = True
            for sample in unaffected:
                vars_unaffected = sample.get_variants()
                if self.match_variant(variant, vars_unaffected):
                    inherited = False
                    break
            if inherited:
                result_variants.append(var)
        if len(unaffected) > 0:
            log_text.append_text(
                'Excluding if start position is in any unaffected samples, {0} variants -> {1} variants'.format(
                    str(len(start_vars)), str(len(result_variants))))
        return result_variants
    def get_variant_fields(self):
        return self.variant_fields
    def get_format_var(self):
        return self.format_var




###
### main.py
###
import sys
# from analysisDriver import Analysis, AnalysisType
# from log import Log
# from preprocessor import Preprocessor
# -----------------------------------------------------------------------------
#   Method:         def main():
#   Description:    main method
#   Parameters:     selected_analysis_type     One of constants for ANALYSIS_TYPE
#                   file_ped                   Single ped file for the dataset
#                   file_path                  Array of annovar files to be analyzed, proband must be at [0]
def main(analysis_type, ped, annotation_files):
    log = Log()
    # 1. parse inputs into samples
    preprocessor = Preprocessor(analysis_type, ped, annotation_files, log)
    # 2. analyze samples
    samples = preprocessor.get_all_samples()
    variant_fields = preprocessor.get_variant_fields()
    analysis_job = Analysis(variant_fields, log, samples)
    final_variants = analysis_job.analyze_variants(analysis_type)
    # 3. append initial and final variant counts to log
    proband_start_vars = preprocessor.get_proband_start_vars()
    log.append_text(
        'Initial proband count = {0} / Final count = {1}'.format(str(proband_start_vars), str(len(final_variants))))
    # 4. write analysis log and final variants to separate files
    log.write_log(analysis_type)
    header = preprocessor.get_header_str()
    log.write_variants(analysis_type, header, final_variants)
# -----------------------------------------------------------------------------
# The user inputs four items into Cavatica: [analysis_types], ped, proband_annotation, [other_annotations]
# Data types: ped and proband_annotation are single files; analysis_types is enum array len >= 1; other_annotations is file array len >=0
# The args are always sent in this order, so proband is the first '.txt' file
# args are received in a one-dimensional array, ie the arrays have each item added individually
# eg sys.argv[1:] = [AR_H, AD_V, path/samples.ped, path/proband.txt, path/annot_1.txt, path/annot_2.txt, path/annot_3.txt]
def get_enums(analyses_strs):
    analyses = []
    for a in analyses_strs:
        if a == 'AR_CH':
            analyses.append(AnalysisType.AR_CH)
        elif a == 'AR_H':
            analyses.append(AnalysisType.AR_H)
        elif a == 'AD_V':
            analyses.append(AnalysisType.AD_V)
        elif a == 'AD_IM':
            analyses.append(AnalysisType.AD_IM)
        elif a == 'AD_NM':
            analyses.append(AnalysisType.AD_NM)
        else:
            raise Exception('Invalid analysis type:', a)
    return analyses
if __name__ == '__main__':
    ped = None
    annotations = []
    analyses_strs = []
    for arg in sys.argv[1:]:
        if arg.endswith('.ped'):
            ped = arg
        elif arg.endswith('.txt'):  # proband is always first text file, becomes annotations[0]
            annotations.append(arg)
        else:
            analyses_strs.append(arg)
    analyses = get_enums(analyses_strs)
    for analysis in analyses:
        main(analysis, ped, annotations)
# -----------------------------------------------------------------------------
# Use this to run phenodb analysis directly with a python command
# Parameters:
#     analysis_types: list of strings, must be options in get_enums()
#     ped: is a filepath to a pedigree
#     annotation_files: list of filepaths to annotation files. Proband must be at index [0]
def main_direct(analyses_strs, ped, annotations):
    analyses = get_enums(analyses_strs)
    for analysis in analyses:
        main(analysis, ped, annotations)
