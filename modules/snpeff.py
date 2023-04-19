__author__ = 'alipirani'
import os
from config_settings import ConfigSectionMap
from modules.log_modules import keep_logging
from modules.logging_subprocess import *
from sys import platform as _platform
import subprocess as subprocess
import errno
from modules.tabix import *

def make_sure_path_exists(out_path):
    """This function checks if the args out_path exists and generates an empty directory if it doesn't.

    :param:
        out_path: Directory path to check or create a new directory.

    :return: null/exception

    """

    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('- Errors in output folder path! please change the output path or analysis name',
                         '- Errors in output folder path! please change the output path or analysis name', logger,
                         'info')
            exit()

def prepare_snpEff_db(reference, vc_logs_folder, logger, Config):
    
    keep_logging('- Preparing snpEff database requirements.', '- Preparing snpEff database requirements.', logger, 'info')

    reference_basename = (os.path.basename(reference)).split(".")

    global bin_dir
    proc = subprocess.Popen(["which snpEff"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    bin_dir = os.path.dirname(out.decode("utf-8"))
    
    ## Great Lakes Changes
    proc = subprocess.Popen(["find $CONDA_PREFIX/share/ -name snpEff.config"], stdout=subprocess.PIPE, shell=True)
    (out2, err2) = proc.communicate()
    
    if out2:
        snpeff_config = (out2.decode("utf-8")).strip()
    else:
        print ("- Unable to find snpEff config file in conda Environment share directory")
        exit()
    
    # os.system("cp %s $CONDA_PREFIX/bin/" % snpeff_config)
    os.system("cp %s %s" % (snpeff_config, bin_dir))
    
    if os.path.isfile("%s/snpEff.config" % bin_dir):
        # os.system("cp %s/%s/snpEff.config %s" % (ConfigSectionMap("bin_path", Config)['binbase'], ConfigSectionMap("snpeff", Config)['snpeff_bin'], args.filter2_only_snp_vcf_dir))
        # keep_logging("cp %s/snpEff.config %s" % (bin_dir, vc_logs_folder),
        #              "cp %s/snpEff.config %s" % (bin_dir, vc_logs_folder), logger, 'debug')
        call("cp %s/snpEff.config %s" % (bin_dir, vc_logs_folder), logger)
    else:
        keep_logging("- Error: %s/snpEff.config doesn't exists.\nExiting..." % bin_dir,
                     "- Error: %s/snpEff.config doesn't exists.\nExiting..." % bin_dir, logger, 'exception')
        exit()
    make_sure_path_exists("%s/data/%s" % (bin_dir, reference_basename[0]))
    make_sure_path_exists("%s/data/genomes/" % bin_dir)
    
    # keep_logging("cp %s %s/data/genomes/%s.fa" % (reference, bin_dir, reference_basename[0]),
    #              "cp %s %s/data/genomes/" % (reference, bin_dir), logger, 'debug')
    call("cp %s %s/data/genomes/%s.fa" % (reference, bin_dir, reference_basename[0]), logger)
    with open("%s/snpEff.config" % vc_logs_folder, "a") as conf_file:
        conf_file.write(
            "\n\n##Building Custom Database###\n%s.genome\t: %s\n\n" % (reference_basename[0], reference_basename[0]))
    conf_file.close()
    # get the gff name from config file
    if os.path.isfile("%s/%s.gff" % (os.path.dirname(reference), reference_basename[0])):
        # keep_logging("cp %s/%s.gff %s/data/%s/genes.gff" % (
        #     os.path.dirname(reference), reference_basename[0], bin_dir, reference_basename[0]),
        #              "cp %s/%s.gff %s/data/%s/genes.gff" % (os.path.dirname(reference), reference_basename[0],
        #                                                     bin_dir,
        #                                                     reference_basename[0]), logger, 'debug')
        # keep_logging("cp %s/%s.gb* %s/data/%s/genes.gbk" % (
        #     os.path.dirname(reference), reference_basename[0], bin_dir, reference_basename[0]),
        #              "cp %s/%s.gff %s/data/%s/genes.gff" % (os.path.dirname(reference), reference_basename[0],
        #                                                     bin_dir,
        #                                                     reference_basename[0]), logger, 'debug')
        call("cp %s/%s.gff %s/data/%s/genes.gff" % (
            os.path.dirname(reference), reference_basename[0], bin_dir, reference_basename[0]), logger)
        call("cp %s/%s.gb* %s/data/%s/genes.gbk" % (
            os.path.dirname(reference), reference_basename[0], bin_dir, reference_basename[0]), logger)
    else:
        keep_logging(
        "- Warning: %s/%s.gff file doesn't exists. Make sure the GFF file has the same prefix as reference fasta file.\n- Warning: This will cause issues in subsequent analyses." % (
        os.path.dirname(reference), reference_basename[0]),
        "- Warning: %s/%s.gff file doesn't exists. Make sure the GFF file has the same prefix as reference fasta file\n- Warning: This will cause issues in subsequent analyses." % (
        os.path.dirname(reference), reference_basename[0]), logger, 'exception')
        #exit()
    
    if os.path.isfile("%s/%s.gbf" % (os.path.dirname(reference), reference_basename[0])):
        # keep_logging("cp %s/%s.gbf %s/data/%s/genes.gbk" % (
        #     os.path.dirname(reference), reference_basename[0], bin_dir, reference_basename[0]),
        #              "cp %s/%s.gbf %s/data/%s/genes.gbk" % (os.path.dirname(reference), reference_basename[0],
        #                                                     bin_dir,
        #                                                     reference_basename[0]), logger, 'debug')
        call("cp %s/%s.gbf %s/data/%s/genes.gbk" % (
            os.path.dirname(reference), reference_basename[0], bin_dir, reference_basename[0]), logger)
    else:
        keep_logging(
            "- Warning: %s/%s.gbf file doesn't exists. Make sure the Genbank file has the same prefix as reference fasta file\n- Warning: This will cause issues in subsequent analyses." % (
            os.path.dirname(reference), reference_basename[0]),
            "- Warning: %s/%s.gff file doesn't exists. Make sure the Genbank file has the same prefix as reference fasta file\n- Warning: This will cause issues in subsequent analyses." % (
            os.path.dirname(reference), reference_basename[0]), logger, 'exception')
        #exit()
    
    call("%s build -genbank -v %s -c %s/snpEff.config -dataDir %s/data" % (
    ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], vc_logs_folder, bin_dir),
         logger)
    keep_logging("%s build -genbank -v %s -c %s/snpEff.config -dataDir %s/data" % (ConfigSectionMap("snpeff", Config)['base_cmd'], reference_basename[0], vc_logs_folder, bin_dir), "", logger, 'debug')
     
def variant_annotation(vcf_file, reference, vc_logs_folder, Config, logger):
    global bin_dir
    proc = subprocess.Popen(["which snpEff"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    bin_dir = os.path.dirname(out.decode("utf-8"))

    reference_basename = (os.path.basename(reference)).split(".")

    if ConfigSectionMap("snpeff", Config)['prebuild'] == "yes":
        if ConfigSectionMap("snpeff", Config)['db']:
            print ("- Using pre-built snpEff database: %s" % ConfigSectionMap("snpeff", Config)['db'])
            ## Great Lakes Changes
            proc = subprocess.Popen(["%s databases | grep %s" % (
            ConfigSectionMap("snpeff", Config)['base_cmd'], ConfigSectionMap("snpeff", Config)['db'])],
                                    stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            if out2:
                snpeffdb = ConfigSectionMap("snpeff", Config)['db']
            else:
                print ("- The database name %s provided was not found. Check the name and try again" % \
                      ConfigSectionMap("snpeff", Config)['db'])
                exit()
        else:
            print ("- snpEff db section is not set in config file")
            exit()
    
    else:
        snpeffdb = reference_basename[0]
    
    
    raw_vcf = vcf_file.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_aln_mpileup_raw.vcf')
    annotate_vcf_cmd = "%s -csvStats %s_ANN.csv -dataDir %s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                        (ConfigSectionMap("snpeff", Config)['base_cmd'], raw_vcf, bin_dir,
                        ConfigSectionMap("snpeff", Config)['snpeff_parameters'], vc_logs_folder,
                        snpeffdb, raw_vcf, raw_vcf)
    
    final_vcf = vcf_file
    annotate_final_vcf_cmd = "%s -csvStats %s_ANN.csv -dataDir %s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                                (ConfigSectionMap("snpeff", Config)['base_cmd'], final_vcf, bin_dir,
                                ConfigSectionMap("snpeff", Config)['snpeff_parameters'],
                                vc_logs_folder, snpeffdb, final_vcf, final_vcf)
    print (annotate_vcf_cmd)
    print (annotate_final_vcf_cmd)
    call(annotate_vcf_cmd, logger)
    call(annotate_final_vcf_cmd, logger)
    files_for_tabix = ['%s_ANN.vcf' % raw_vcf, '%s_ANN.vcf' % final_vcf]
    tabix(files_for_tabix, "vcf", logger, Config)

def indel_annotation(vcf_file, reference, vc_logs_folder, Config, logger):
    global bin_dir
    proc = subprocess.Popen(["which snpEff"], stdout=subprocess.PIPE, shell=True)
    (out, err) = proc.communicate()
    bin_dir = os.path.dirname(out.decode("utf-8"))

    reference_basename = (os.path.basename(reference)).split(".")

    if ConfigSectionMap("snpeff", Config)['prebuild'] == "yes":
        if ConfigSectionMap("snpeff", Config)['db']:
            print ("- Using pre-built snpEff database: %s" % ConfigSectionMap("snpeff", Config)['db'])
            proc = subprocess.Popen(["%s databases | grep %s" % (
            ConfigSectionMap("snpeff", Config)['base_cmd'], ConfigSectionMap("snpeff", Config)['db'])],
                                    stdout=subprocess.PIPE, shell=True)
            (out2, err2) = proc.communicate()
            if out2:
                snpeffdb = ConfigSectionMap("snpeff", Config)['db']
            else:
                print ("- The database name %s provided was not found. Check the name and try again" % \
                      ConfigSectionMap("snpeff", Config)['db'])
                exit()
        else:
            print ("- snpEff db section is not set in config file")
            exit()
    else:
        snpeffdb = reference_basename[0]

    final_vcf = vcf_file.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf')
    annotate_final_vcf_cmd = "%s -csvStats %s_ANN.csv -dataDir %s/data/ %s -c %s/snpEff.config %s %s > %s_ANN.vcf" % \
                                (ConfigSectionMap("snpeff", Config)['base_cmd'], final_vcf, bin_dir,
                                ConfigSectionMap("snpeff", Config)['snpeff_parameters'],
                                vc_logs_folder, snpeffdb, final_vcf, final_vcf)
    print (annotate_final_vcf_cmd)
    call(annotate_final_vcf_cmd, logger)

    files_for_tabix = ['%s_ANN.vcf' % final_vcf]
    tabix(files_for_tabix, "vcf", logger, Config)
