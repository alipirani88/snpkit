__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from sys import platform as _platform
from modules.log_modules import keep_logging
from modules.logging_subprocess import *

def vcfstats(final_raw_vcf, out_path, analysis):
    base_tabix = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['tabix_bin']
    base_vcftools_perl = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_perl_bin']
    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    # bgzip_cmd = "%s/bgzip %s" % (base_tabix, final_raw_vcf)
    # os.system(bgzip_cmd)
    # tabix_cmd = "%s/tabix %s.gz" % (base_tabix, final_raw_vcf)
    # os.system(tabix_cmd)
    vcfstat_cmd = "%s/vcf-stats %s.gz > %s/%s_vcf_stats" % (base_vcftools_perl, final_raw_vcf, out_path, analysis)
    print "\nRunning:\n [%s] \n" % vcfstat_cmd
    os.system(vcfstat_cmd)
    vcf_stats_file = "%s/%s_vcf_stats" % (out_path, analysis)
    return vcf_stats_file



def only_snp_filter1_vcf(gatk_filter1_final_vcf, out_path, analysis, reference):
    base_tabix = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['tabix_bin']
    base_vcftools_perl = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_perl_bin']
    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    onlysnp_filter1_cmd = "%s/vcftools --vcf %s --remove-indels --recode --recode-INFO-all --out %s/%s_filter1_onlysnp.vcf" % (base_vcftools_bin, gatk_filter1_final_vcf, out_path, analysis)
    #print onlysnp_cmd
    print "Running Command: [%s]" % onlysnp_filter1_cmd
    os.system(onlysnp_filter1_cmd)
    only_snp_filter1_vcf_file = "%s/%s_filter1_onlysnp.vcf.recode.vcf" % (out_path, analysis)
    return only_snp_filter1_vcf_file

def only_snp_filter2_vcf(gatk_filter2_final_vcf, out_path, analysis, reference):
    base_tabix = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['tabix_bin']
    base_vcftools_perl = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_perl_bin']
    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    onlysnp_filter2_cmd = "%s/vcftools --vcf %s --remove-indels --recode --recode-INFO-all --out %s/%s_filter2_onlysnp.vcf" % (base_vcftools_bin, gatk_filter2_final_vcf, out_path, analysis)
    #print onlysnp_cmd
    print "Running Command: [%s]" % onlysnp_filter2_cmd
    os.system(onlysnp_filter2_cmd)
    only_snp_filter2_vcf_file = "%s/%s_filter2_onlysnp.vcf.recode.vcf" % (out_path, analysis)
    return only_snp_filter2_vcf_file


def vcftools_vcf2fasta_filter1(only_snp_filter1_vcf_file, out_path, analysis, reference):
    base_tabix = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['tabix_bin']
    base_vcftools_perl = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_perl_bin']
    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    #print only_snp_filter1_vcf_file
    bgzip_cmd = "%s/bgzip -f %s" % (base_tabix, only_snp_filter1_vcf_file)
    os.system(bgzip_cmd)
    tabix_cmd = "%s/tabix %s.gz" % (base_tabix, only_snp_filter1_vcf_file)
    os.system(tabix_cmd)
    vcftools_vcf2fasta_filter1_cmd = "cat %s | %s/vcf-consensus %s.gz > %s/%s_filter1_consensus.fa" % (reference, base_vcftools_bin, only_snp_filter1_vcf_file, out_path, analysis)
    print "Running Command: [%s]" % vcftools_vcf2fasta_filter1_cmd
    #os.system(vcftools_vcf2fasta_filter1_cmd)
    #print "Consensus fasta file created"
    if _platform == "darwin":
        change_header_cmd = "sed -i '' 's/>.*/>%s/g' %s/%s_filter1_consensus.fa" % (analysis, out_path, analysis)
        os.system(change_header_cmd)
    else:
        change_header_cmd = "sed -i 's/>.*/>%s/g' %s/%s_filter1_consensus.fa" % (analysis, out_path, analysis)
        os.system(change_header_cmd)


def vcftools_vcf2fasta_filter2(only_snp_filter2_vcf, out_path, analysis, reference, logger, Config):
    base_tabix = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['tabix_bin']
    base_vcftools_perl = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_perl_bin']
    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    bgzip_cmd = "%s/bgzip -f %s" % (base_tabix, only_snp_filter2_vcf)
    os.system(bgzip_cmd)
    tabix_cmd = "%s/tabix %s.gz" % (base_tabix, only_snp_filter2_vcf)
    os.system(tabix_cmd)
    vcftools_vcf2fasta_filter2_cmd = "cat %s | %s/vcf-consensus %s.gz > %s_filter2_consensus.fa" % (reference, base_vcftools_bin, only_snp_filter2_vcf, only_snp_filter2_vcf)
    keep_logging(vcftools_vcf2fasta_filter2_cmd, vcftools_vcf2fasta_filter2_cmd, logger, 'debug')
    try:
        call(vcftools_vcf2fasta_filter2_cmd, logger)
        #print ""
    except sp.CalledProcessError:
        keep_logging('Error in vcftools vcf 2 fasta step. Exiting.', 'Error in vcftools vcf 2 fasta step. Exiting.', logger, 'exception')
        sys.exit(1)
    bash_script_file = "%s.sh" % (only_snp_filter2_vcf)
    f1=open(bash_script_file, 'w+')
    f1.write(vcftools_vcf2fasta_filter2_cmd)
    bash_command = "bash %s" % bash_script_file
    os.system(bash_command)
    if _platform == "darwin":
        change_header_cmd = "sed -i '' 's/>.*/>%s/g' %s_filter2_consensus.fa" % (analysis, only_snp_filter2_vcf)
        os.system(change_header_cmd)
    else:
        change_header_cmd = "sed -i 's/>.*/>%s/g' %s_filter2_consensus.fa" % (analysis, only_snp_filter2_vcf)
        os.system(change_header_cmd)

#base_tabix = ConfigSectionMap("bin_path")['binbase'] + "/" + ConfigSectionMap("vcftools")['tabix_bin']
#base_vcftools = ConfigSectionMap("bin_path")['binbase'] + "/" + ConfigSectionMap("vcftools")['vcftools_bin']


def vcfstats(final_raw_vcf, out_path, analysis, logger, Config):
    base_tabix = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['tabix_bin']
    base_vcftools_perl = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_perl_bin']
    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    bgzip_cmd = "%s/bgzip -f %s" % (base_tabix, final_raw_vcf)
    os.system(bgzip_cmd)
    tabix_cmd = "%s/tabix %s.gz" % (base_tabix, final_raw_vcf)
    os.system(tabix_cmd)
    vcfstat_cmd = "%s/vcf-stats %s.gz > %s/%s_vcf_stats" % (base_vcftools_bin, final_raw_vcf, out_path, analysis)
    keep_logging(vcfstat_cmd, vcfstat_cmd, logger, 'debug')
    try:
        call(vcfstat_cmd, logger)
        #print ""
    except sp.CalledProcessError:
        keep_logging('Error in vcftools vcf stats step. Exiting.', 'Error in vcftools vcf stats step. Exiting.', logger, 'exception')
        sys.exit(1)
    vcf_stats_file = "%s/%s_vcf_stats" % (out_path, analysis)
    keep_logging('VCF Stats file: {}'.format(vcf_stats_file), 'VCF Stats file: {}'.format(vcf_stats_file), logger, 'debug')
    return vcf_stats_file



def only_snp_raw_vcf(final_raw_vcf, out_path, analysis, reference):
    base_tabix = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['tabix_bin']
    base_vcftools_perl = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_perl_bin']
    base_vcftools_bin = ConfigSectionMap("bin_path", Config)['binbase'] + "/" + ConfigSectionMap("vcftools", Config)['vcftools_bin']
    onlysnp_raw_cmd = "%s/vcftools --vcf %s --remove-indels --recode --recode-INFO-all --out %s/%s_raw_onlysnp.vcf" % (base_vcftools_bin, final_raw_vcf, out_path, analysis)
    #print onlysnp_cmd
    print "Running Command: [%s]" % onlysnp_raw_cmd
    os.system(onlysnp_raw_cmd)
    only_snp_raw_vcf_file = "%s/%s_raw_onlysnp.vcf.recode.vcf" % (out_path, analysis)
    return only_snp_raw_vcf_file




