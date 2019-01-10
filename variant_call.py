__author__ = 'alipirani'

import sys
import os
import argparse
import errno
import ConfigParser
import glob
import re
import multiprocessing
import subprocess as subprocess
from datetime import datetime
from joblib import Parallel, delayed
from config_settings import ConfigSectionMap
from modules.logging_subprocess import *
from modules.log_modules import *
from argparse import RawTextHelpFormatter
from modules.phage_detection import *
from modules.variant_diagnostics.find_repeats import *
from modules.variant_diagnostics.mask_regions import *

""" Command Line Argument Parsing """
def parser():
    parser = argparse.ArgumentParser(description='\nVariant Calling pipeline for Illumina PE/SE data.\n', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-type', action='store', dest="type", help='Type of reads: SE or PE', required=True)
    required.add_argument('-readsdir', action='store', dest="dir", help='Path to Sequencing Reads Data directory. NOTE: Provide full/absolute path.', required=True)
    required.add_argument('-outdir', action='store', dest="output_folder", help='Output Folder Path ending with output directory name to save the results. Creates a new output directory path if it doesn\'t exist. NOTE: Provide full/absolute path.', required=True)
    required.add_argument('-index', action='store', dest="index", help='Reference Index Name. Most Frequently used reference genomes index options: KPNIH1 | MRSA_USA_300 | MRSA_USA_100 | CDIFF_630 | paris Make sure the paths are properly set in config file', required=True)
    required.add_argument('-steps', action='store', dest="steps", help='Variant Calling Steps in sequential order.\n'
                                                                     '1.   All: Run all variant calling  steps starting from trimming the reads, mapping, post-processing the alignments and calling variants;\n'
                                                                     '2.   core_All: Extract core snps and generate different types of alignments, SNP/Indel Matrices and diagnostics plots.')
    required.add_argument('-analysis', action='store', dest="analysis_name", help='Unique analysis name that will be used as prefix to saving results and log files.', required=True)
    optional.add_argument('-config', action='store', dest="config", help='Path to Config file, Make sure to check config settings before running pipeline', required=False)
    optional.add_argument('-suffix', action='store', dest="suffix", help='Fastq reads suffix such as fastq, fastq.gz, fq.gz, fq; Default: fastq.gz', required=False)
    optional.add_argument('-filenames', action='store', dest="filenames", help='fastq filenames with one single-end filename per line. \nIf the type is set to PE, it will detect the second paired-end filename with the suffix from first filename. \nUseful for running variant calling pipeline on selected files in a reads directory or extracting core snps for selected samples in input reads directory. \nOtherwise the pipeline will consider all the samples available in reads directory.', required=False)
    optional.add_argument('-cluster', action='store', dest='cluster', help='Run variant calling pipeline in one of the four modes. Default: local. Suggested mode for core snp is \"cluster\" that will run all the steps in parallel with the available cores. Make sure to provide a large memory node for this option\nThe possible modes are: cluster/parallel-local/local\nSet your specific hpc cluster parameters in config file under the [scheduler] section. Supports only PBS scheduling system. ')
    optional.add_argument('-clean', action="store_true", help='clean up intermediate files. Default: OFF')
    optional.add_argument('-extract_unmapped', action='store', dest="extract_unmapped", help='Extract unmapped reads, assemble it and detect AMR genes using ariba')
    optional.add_argument('-datadir', action='store', dest="datadir", help='Path to snpEff data directory')
    optional.add_argument('-snpeff_db', action='store', dest="snpeff_db", help='Name of pre-build snpEff database to use for Annotation')
    optional.add_argument('-debug_mode', action='store', dest="debug_mode", help='yes/no for debug mode')
    optional.add_argument('-gubbins', action='store', dest="gubbins", help='yes/no for running gubbins')
    return parser


""" Sanity checks and directory structure maintenance methods """
def file_exists(path1):
    if not os.path.isfile(path1):
        file_basename = os.path.basename(path1)
        keep_logging('The file {} does not exists. Please provide another file with full path or check the files path.\n'.format(file_basename), 'The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
        exit()

def make_sure_path_exists(out_path):
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
            exit()

def get_filenames(dir, type, filenames, analysis, suffix):
    if not filenames:
        if not suffix:
            suffix = ".fastq.gz"
        try:
            list_of_files = glob.glob("%s/*%s" % (dir, suffix))
            if len(list_of_files) < 1:
                keep_logging('No fastq files with suffix %s found in reads directory %s' % (suffix, dir), 'No fastq files with suffix %s found in reads directory %s' % (suffix, dir), logger, 'info')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                keep_logging('Error while listing files in reads directory.', 'Error while listing files in reads directory.', logger, 'exception')
                exit()
    else:
        list_of_files = []
        with open(filenames) as fp:
            for line in fp:
                line = line.strip()
                line = dir + "/" + line
                list_of_files.append(line)
    return list_of_files

""" Methods to generate jobs for various pipeline tasks """
def create_varcall_jobs(filenames_array, type, output_folder, reference, steps, config_file, logger):
    jobs_temp_dir = "%s/temp_jobs" % output_folder
    make_sure_path_exists(jobs_temp_dir)
    keep_logging('Generating cluster jobs in temporary directory %s' % jobs_temp_dir, 'Generating cluster jobs in temporary directory %s' % jobs_temp_dir, logger, 'exception')
    ##Change the email address; nodes and processor requirements accordingly

    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])

    for file in filenames_array:
        filename_base = os.path.basename(file)
        #print filename_base
        if "R1_001_final.fastq.gz" in filename_base or "R1.fastq.gz" in filename_base or "1_combine.fastq.gz" in filename_base or "1_sequence.fastq.gz" in filename_base or "_forward.fastq.gz" in filename_base or "R1_001.fastq.gz" in filename_base or "_1.fastq.gz" in filename_base or ".1.fastq.gz" in filename_base or "_R1.fastq.gz" in filename_base or "_L001_R1_001.fastq.gz" in filename_base:
            # Forward reads file name and get analysis name from its name
            first_file = file
            # Get the name of reverse reads files
            if "R1_001_final.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "R1_001.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "1_combine.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "1_sequence.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_forward.fastq.gz" in filename_base:
                second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "R1_001.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif "_1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
            elif ".1.fastq.gz" in filename_base:
                second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
                first_part_split = filename_base.split('.1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)

            """ Have a standard filename preparation step"""
            # else:
            #     print "Using Standard second file naming convention"
            #     second_part = filename_base.replace("_R1_", "_R2_")
            #     first_part_split = filename_base.split('_R1.fastq.gz')
            #     first_part = first_part_split[0].replace('_L001', '')
            #     first_part = re.sub("_S.*_", "", first_part)
            second_file = args.dir + "/" + second_part
            job_name = jobs_temp_dir + "/" + first_part + ".pbs"
            if not steps:
                steps == "All"
            if type == "SE":
                if args.clean:
                    command = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/pipeline.py -PE1 %s -o %s/%s -analysis %s -index %s -type SE -config %s -steps %s -clean" % (
                    os.path.dirname(os.path.abspath(__file__)), first_file, output_folder, first_part, first_part,
                    reference, config_file, steps)
                else:
                    command = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/pipeline.py -PE1 %s -o %s/%s -analysis %s -index %s -type SE -config %s -steps %s" % (os.path.dirname(os.path.abspath(__file__)), first_file, output_folder, first_part, first_part, reference, config_file, steps)

            else:
                if args.clean:
                    command = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/pipeline.py -PE1 %s -PE2 %s -o %s/%s -analysis %s -index %s -type PE -config %s -steps %s -clean" % (
                    os.path.dirname(os.path.abspath(__file__)), first_file, second_file, output_folder, first_part,
                    first_part, reference, config_file, steps)
                else:
                    command = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/pipeline.py -PE1 %s -PE2 %s -o %s/%s -analysis %s -index %s -type PE -config %s -steps %s" % (os.path.dirname(os.path.abspath(__file__)), first_file, second_file, output_folder, first_part, first_part, reference, config_file, steps)

            with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                #out.write(cd_command+'\n') ## changed it to automatically change to PBS working directory
                out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
                out.write("echo $PBS_O_WORKDIR" + '\n')
                out.write("cd %s/temp_jobs" % output_folder + '\n')
                out.write(command+'\n')
        elif "R2_001_final.fastq.gz" in filename_base or "R2.fastq.gz" in filename_base or "2_combine.fastq.gz" in filename_base or "2_sequence.fastq.gz" in filename_base or "_reverse.fastq.gz" in filename_base or "R2_001.fastq.gz" in filename_base or "_2.fastq.gz" in filename_base or ".2.fastq.gz" in filename_base or "_R2.fastq.gz" in filename_base:
            continue
        else:
            keep_logging('Error while generating cluster jobs. Make sure the fastq filenames ends with one of these suffix: R1_001_final.fastq.gz, R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, R1_001.fastq.gz, _1.fastq.gz, .1.fastq.gz, _R1.fastq.gz', 'Error while generating cluster jobs. Make sure the fastq filenames ends with one of these suffix: R1_001_final.fastq.gz, R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, R1_001.fastq.gz, _1.fastq.gz, .1.fastq.gz, _R1.fastq.gz', logger, 'exception')
            print filename_base
            exit()
    list_of_jobs = glob.glob("%s/*.pbs" % jobs_temp_dir)
    return list_of_jobs

def generate_custom_vcf_file_list(filenames_array, logger):
    keep_logging('Generating custom vcf file list for core pipeline steps', 'Generating custom vcf files list for core pipeline steps', logger, 'exception')
    list_of_vcf_files = []
    for file in filenames_array:
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base or "R1.fastq.gz" in filename_base or "1_combine.fastq.gz" in filename_base or "1_sequence.fastq.gz" in filename_base or "_forward.fastq.gz" in filename_base or "R1_001.fastq.gz" in filename_base or "_1.fastq.gz" in filename_base or ".1.fastq.gz" in filename_base or "_R1.fastq.gz" in filename_base:
            # Forward reads file name and get analysis name from its name
            first_file = file
            if "R1_001_final.fastq.gz" in filename_base:
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            elif "_R1.fastq.gz" in filename_base:
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            elif "R1.fastq.gz" in filename_base:
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            elif "1_combine.fastq.gz" in filename_base:
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            elif "1_sequence.fastq.gz" in filename_base:
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            elif "_forward.fastq.gz" in filename_base:
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            elif "R1_001.fastq.gz" in filename_base:
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            elif "_1.fastq.gz" in filename_base:
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            elif ".1.fastq.gz" in filename_base:
                first_part_split = filename_base.split('.1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S.*_", "", first_part)
                first_part = first_part + "_filter2_final.vcf_no_proximate_snp.vcf"
            list_of_vcf_files.append(first_part)
    return list_of_vcf_files

def run_command(job):
    keep_logging('Running Job: bash %s' % job, 'Running Job: bash %s' % job, logger, 'info')
    call("bash %s" % job, logger)
    done = "Job Run completed: %s" % job
    return done

def run_command_list(command):
    keep_logging('Running command: %s' % command, 'Running command: %s' % command, logger, 'info')
    call("%s" % command, logger)
    done = "Command Run completed: %s" % command
    return done

def run_varcall_jobs(list_of_jobs, cluster, log_unique_time, analysis_name, output_folder, logger):
    #Generate command list to run on single cluster or in parallel
    command_list = ""
    command_list_qsub = []
    for job in list_of_jobs:
        command_list = command_list + "bash %s\n" % job
        command_list_qsub.append(job)

    #cluster mode
    if cluster == "cluster":
        # Removing support for single cluster mode and replacing it with parallel-cluster
        # keep_logging('Running Job in cluster mode', 'Running Job in cluster mode', logger, 'info')
        # cluster_job_dir = output_folder + "/temp_jobs/cluster/"
        # make_sure_path_exists(cluster_job_dir)
        # Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
        #               % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])
        # job_name = cluster_job_dir + log_unique_time + "_" + analysis_name + ".pbs"
        # with open(job_name, 'w') as out:
        #     job_title = "#PBS -N %s_%s" % (log_unique_time, analysis_name)
        #     out.write(job_title+'\n')
        #     out.write(Pbs_model_lines+'\n')
        #     out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
        #     out.write("echo $PBS_O_WORKDIR" + '\n')
        #     out.write("cd %s/temp_jobs" % output_folder + '\n')
        #     out.write(command_list+'\n')
        # out.close()
        # keep_logging('Submitting single cluster Job: qsub %s' % job_name, 'Submitting single cluster Job: qsub %s' % job_name, logger, 'info')
        # #call("qsub %s" % job_name, logger)
        # qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        # print qid.split('.')[0]
        # keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')

        keep_logging('Running Jobs in cluster mode', 'Running Jobs in cluster mode', logger, 'info')
        for job in command_list_qsub:
            keep_logging('Submitting Job: qsub %s' % job, 'Submitting Job: qsub %s' % job, logger, 'info')
            call("qsub %s" % job, logger)
        # keep_logging('You can check the job status with: qstat -u USERNAME',
        #              'You can check the job status with: qstat -u USERNAME', logger, 'info')
    elif cluster == "parallel-cluster":
        keep_logging('Running Jobs in parallel-cluster mode', 'Running Jobs in parallel-cluster mode', logger, 'info')
        for job in command_list_qsub:
            keep_logging('Submitting Job: qsub %s' % job, 'Submitting Job: qsub %s' % job, logger, 'info')
            call("qsub %s" % job, logger)
        #keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')
    elif cluster == "parallel-local":
        keep_logging('Running Jobs in parallel-local mode', 'Running Jobs in parallel-local mode', logger, 'info')
        cluster_job_dir = output_folder + "/temp_jobs/parallel-local/"
        make_sure_path_exists(cluster_job_dir)
        job_name = cluster_job_dir + log_unique_time + "_" + analysis_name + ".sh"
        with open(job_name, 'w') as out:
            out.write(command_list+'\n')
        num_cores = multiprocessing.cpu_count()
        keep_logging('Number of cores available: %s' % num_cores, 'Number of cores available: %s' % num_cores, logger, 'info')
        results = Parallel(n_jobs=num_cores)(delayed(run_command)(i) for i in command_list_qsub)
    elif cluster == "local":
        keep_logging('Running Jobs in local mode', 'Running Jobs in local mode', logger, 'info')
        for job in command_list_qsub:
            keep_logging('Running Job: bash %s' % job, 'Running Job: bash %s' % job, logger, 'info')
            call("bash %s" % job, logger)

""" Pipeline individual task methods """
def run_core_prep_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger, config_file):
    file_exists(reference)
    if args.debug_mode == "yes":
        core_prep_pipeline = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/modules/variant_diagnostics/core_pipeline_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 1 -jobrun %s -config %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, config_file)
    else:
        core_prep_pipeline = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 1 -jobrun %s -config %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, config_file)

    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"

    # Changed on 04/11/2018
    # Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=92:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
    #                   % (ConfigSectionMap("scheduler", Config)['email'],
    #                      ConfigSectionMap("scheduler", Config)['notification'],
    #                      ConfigSectionMap("scheduler", Config)['queue'],
    #                      ConfigSectionMap("scheduler", Config)['flux_account'])

    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                      % (ConfigSectionMap("scheduler", Config)['email'],
                         ConfigSectionMap("scheduler", Config)['notification'],
                         ConfigSectionMap("scheduler", Config)['resources'],
                         ConfigSectionMap("scheduler", Config)['queue'],
                         ConfigSectionMap("scheduler", Config)['flux_account'])

    # with open(job_name, 'w') as out:
    #     job_title = "#PBS -N %s_%s_core" % (log_unique_time, analysis_name)
    #     out.write(job_title+'\n')
    #     out.write(Pbs_model_lines+'\n')
    #     out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
    #     out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
    #     out.write("cd %s" % core_temp_dir + '\n')
    #     out.write(core_prep_pipeline+'\n')
    # out.close()
    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        keep_logging('Running in parallel-local mode: bash %s' % job_name, 'Running in parallel-local mode: bash %s' % job_name, logger, 'info')
        #call("bash %s" % job_name, logger)
        os.system("bash %s" % job_name)
    elif cluster == "cluster":
        #call("qsub %s" % job_name, logger)
        keep_logging('Generating core_prep command: \n%s\n' % core_prep_pipeline, 'Generating core_prep command: \n%s\n' % core_prep_pipeline, logger, 'info')
        #qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        #print qid.split('.')[0]
    elif cluster == "parallel-cluster":
        #call("qsub %s" % job_name, logger)
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s_%s_core" % (log_unique_time, analysis_name)
            out.write(job_title + '\n')
            out.write(Pbs_model_lines + '\n')
            out.write(
                "#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
            out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
            out.write("cd %s" % core_temp_dir + '\n')
            out.write(core_prep_pipeline + '\n')
        out.close()
        keep_logging('Submitting parallel-cluster Job: qsub %s' % job_name, 'Submitting parallel-cluster Job: qsub %s' % job_name, logger, 'info')
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
    #keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')

    return core_prep_pipeline

def run_core_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger, core_results_dir, config_file):
    file_exists(reference)
    if args.debug_mode == "yes":
        core_pipeline = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/modules/variant_diagnostics/core_pipeline_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 2 -jobrun %s -results_dir %s -config %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file)
    else:
        core_pipeline = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 2 -jobrun %s -results_dir %s -config %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file)
    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"

    # Changed on 11/04/2018
    # Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,mem=47000mb,walltime=92:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
    #                   % (ConfigSectionMap("scheduler", Config)['email'],
    #                      ConfigSectionMap("scheduler", Config)['notification'],
    #                      ConfigSectionMap("scheduler", Config)['queue'],
    #                      ConfigSectionMap("scheduler", Config)['flux_account'])

    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                      % (ConfigSectionMap("scheduler", Config)['email'],
                         ConfigSectionMap("scheduler", Config)['notification'],
                         ConfigSectionMap("scheduler", Config)['large_resources'],
                         ConfigSectionMap("scheduler", Config)['queue'],
                         ConfigSectionMap("scheduler", Config)['flux_account'])
    
    # with open(job_name, 'w') as out:
    #     job_title = "#PBS -N %s_%s_core" % (log_unique_time, analysis_name)
    #     out.write(job_title+'\n')
    #     out.write(Pbs_model_lines+'\n')
    #     out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
    #     out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
    #     out.write("cd %s" % core_temp_dir + '\n')
    #     out.write(core_pipeline+'\n')
    # out.close()
    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        #call("qsub %s" % job_name, logger)
        keep_logging('Generating Core Step command: \n%s\n' % core_pipeline, 'Generating Core Step command: \n%s\n' % core_pipeline, logger, 'info')
        #qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        #print qid.split('.')[0]
    elif cluster == "parallel-cluster":
        #call("qsub %s" % job_name, logger)
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s_%s_core" % (log_unique_time, analysis_name)
            out.write(job_title + '\n')
            out.write(Pbs_model_lines + '\n')
            out.write(
                "#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
            out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
            out.write("cd %s" % core_temp_dir + '\n')
            out.write(core_pipeline + '\n')
        keep_logging('Submitting parallel-cluster Job: qsub %s' % job_name, 'Submitting parallel-cluster Job: qsub %s' % job_name, logger, 'info')
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
    #keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')
    return core_pipeline

def run_report_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger, core_results_dir, config_file):
    file_exists(reference)
    if args.debug_mode == "yes":
        core_pipeline = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/modules/variant_diagnostics/core_pipeline_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 3 -jobrun %s -results_dir %s -config %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file)
    else:
        core_pipeline = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 3 -jobrun %s -results_dir %s -config %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file)
    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"
    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=92:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])
    # with open(job_name, 'w') as out:
    #     job_title = "#PBS -N %s_%s_core" % (log_unique_time, analysis_name)
    #     out.write(job_title+'\n')
    #     out.write(Pbs_model_lines+'\n')
    #     out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
    #     out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
    #     out.write("cd %s" % core_temp_dir + '\n')
    #     out.write(core_pipeline+'\n')
    # out.close()

    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        #call("qsub %s" % job_name, logger)
        keep_logging('Generating Report Step command: \n%s\n' % core_pipeline, 'Generating Report Step command: \n%s\n' % core_pipeline, logger, 'info')
        #qid = subprocess.check_output("bash %s" % job_name, shell=True)
        #print qid.split('.')[0]
    elif cluster == "parallel-cluster":
        #call("qsub %s" % job_name, logger)
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s_%s_core" % (log_unique_time, analysis_name)
            out.write(job_title + '\n')
            out.write(Pbs_model_lines + '\n')
            out.write(
                "#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
            out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
            out.write("cd %s" % core_temp_dir + '\n')
            out.write(core_pipeline + '\n')
        out.close()
        keep_logging('Submitting parallel-cluster Job: qsub %s' % job_name, 'Submitting parallel-cluster Job: qsub %s' % job_name, logger, 'info')
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
    #keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')
    return core_pipeline

def run_tree_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger, core_results_dir, config_file):
    if args.debug_mode == "yes":
        core_pipeline = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/modules/variant_diagnostics/core_pipeline_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 4 -jobrun %s -results_dir %s -config %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file)
        if args.gubbins == "yes":
            core_pipeline = core_pipeline + " -gubbins %s" % args.gubbins
    else:
        core_pipeline = "/nfs/esnitkin/bin_group/anaconda2/bin/python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 4 -jobrun %s -results_dir %s -config %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file)
        if args.gubbins == "yes":
            core_pipeline = core_pipeline + " -gubbins %s" % args.gubbins




    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"
    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])
    # with open(job_name, 'w') as out:
    #     job_title = "#PBS -N %s_%s_core_tree" % (log_unique_time, analysis_name)
    #     out.write(job_title+'\n')
    #     out.write(Pbs_model_lines+'\n')
    #     out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
    #     out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
    #     out.write("cd %s" % core_temp_dir + '\n')
    #     out.write(core_pipeline+'\n')
    # out.close()

    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        #call("qsub %s" % job_name, logger)
        keep_logging('Generating Tree Step command: \n%s\n' % core_pipeline, 'Generating Tree Step command: \n%s\n' % core_pipeline, logger, 'info')
        #qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        #print qid.split('.')[0]
    elif cluster == "parallel-cluster":
        #call("qsub %s" % job_name, logger)
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s_%s_core_tree" % (log_unique_time, analysis_name)
            out.write(job_title + '\n')
            out.write(Pbs_model_lines + '\n')
            out.write(
                "#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
            out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
            out.write("cd %s" % core_temp_dir + '\n')
            out.write(core_pipeline + '\n')
        out.close()
        keep_logging('Submitting parallel-cluster Job: qsub %s' % job_name, 'Submitting parallel-cluster Job: qsub %s' % job_name, logger, 'info')
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
    #keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')
    return core_pipeline

""" Start of Main Method/Pipeline """
if __name__ == '__main__':

    # Set up logging modules and config file
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()

    # Pass arguments to args object
    args = parser().parse_args()

    global config_file
    global log_unique_time
    global Config
    global files_to_delete
    global logger

    if args.output_folder != '':
        args.output_folder += '/'
    make_sure_path_exists(args.output_folder)

    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"


    logs_folder = args.output_folder + "/Logs"
    make_sure_path_exists(logs_folder)
    # logger = generate_logger(logs_folder, args.analysis_name, log_unique_time)

    files_to_delete = []
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)

    # Run pipeline steps
    if "core_All" not in args.steps and "core" not in args.steps and "core_prep" not in args.steps and "report" not in args.steps and "tree" not in args.steps and "2" not in args.steps:
        """ Set Up variant calling logs folder/logger object, cluster mode and copy config files to it"""
        vc_logs_folder = logs_folder + "/variant_calling"
        make_sure_path_exists(vc_logs_folder)
        logger = generate_logger(vc_logs_folder, args.analysis_name, log_unique_time)
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, vc_logs_folder, log_unique_time, args.analysis_name), logger)
        if args.cluster:
            cluster_mode = args.cluster
        else:
            cluster_mode = "local"

        keep_logging('START: Variant Calling Pipeline', 'START: Variant Calling Pipeline', logger, 'info')

        """ Main Variant calling Methods: Generate and Run the jobs"""
        list_of_files = get_filenames(args.dir, args.type, args.filenames, args.analysis_name, args.suffix)
        list_of_jobs = create_varcall_jobs(list_of_files, args.type, args.output_folder, args.index, args.steps, config_file, logger)
        run_varcall_jobs(list_of_jobs, cluster_mode, log_unique_time, args.analysis_name, args.output_folder, logger)
        time_taken = datetime.now() - start_time_2
        keep_logging('Logs were recorded in file with extension log.txt in %s' % vc_logs_folder, 'Logs were recorded in file with extension log.txt in %s' % vc_logs_folder, logger, 'info')
        keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')
        keep_logging('End: Variant calling Pipeline', 'End: Variant calling Pipeline', logger, 'info')

    elif "core_All" in args.steps or "2" in args.steps:
        core_All_cmds = []
        """ Set Up Core Prep logs folder/logger object, cluster mode and copy config files to it"""
        core_prep_logs_folder = logs_folder + "/core_prep"
        make_sure_path_exists(core_prep_logs_folder)
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        logger = generate_logger(core_prep_logs_folder, args.analysis_name, log_unique_time)
        #keep_logging('START: Copying variant calling results to %s and running core_prep step.' % core_temp_dir, 'START: Copying variant calling results to %s and running core_prep step.' % core_temp_dir, logger, 'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, core_prep_logs_folder, log_unique_time, args.analysis_name), logger)
        make_sure_path_exists(core_temp_dir)
        keep_logging('\nCopying vcf files to %s\n' % core_temp_dir, '\nCopying vcf files to %s\n' % core_temp_dir, logger, 'info')
        cp_command = "cp %s/*/*_vcf_results/*_filter2_indel_final.vcf %s/*/*_vcf_results/*_aln_mpileup_raw.vcf %s/*/*_vcf_results/*_raw.vcf_5bp_indel_removed.vcf* %s/*/*_vcf_results/*filter2_final.vcf* %s/*/*_vcf_results/*vcf_no_proximate_snp.vcf* %s/*/*_vcf_results/*array %s/*/*unmapped.bed_positions %s/*/*_vcf_results/*_indel_gatk.vcf %s/*/*_stats_results/*_depth_* %s/*/*_stats_results/*_markduplicates_metrics %s/*/*_stats_results/*_markduplicates_metrics %s" % (args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, core_temp_dir)


        # cp_command = "cp %s/*/*_filter2_indel_final.vcf %s/*/*_aln_mpileup_raw.vcf %s/*/*_raw.vcf_5bp_indel_removed.vcf* %s/*/*filter2_final.vcf* %s/*/*vcf_no_proximate_snp.vcf* %s/*/*array %s/*/*unmapped.bed_positions %s" % (
        # args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder,
        # args.output_folder, args.output_folder, core_temp_dir)
        call(cp_command, logger)

        """ Decompress zipped files in core temp folder"""
        list_of_gzipped_files = glob.glob("%s/*.gz" % core_temp_dir)
        keep_logging('Decompressing gzipped files in %s\n' % core_temp_dir, 'Decompressing gzipped files in %s\n' % core_temp_dir, logger, 'info')
        gzipped_command_list = []
        for i in list_of_gzipped_files:
            gzipped_command_list.append("gzip -df %s" % i)

        """ Set the number of cores variable for Parallel-local mode """
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command_list)(i) for i in gzipped_command_list)
        filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']

        """ Generate a custom vcf file list to process for core prep step. If file is not provided, it will consider all the samples in output folder"""
        """ Perform Sanity check to confirm that variant calling results are present in core temp folder """
        if args.filenames:
            # fix this bug
            list_of_files = get_filenames(args.dir, args.type, args.filenames, args.analysis_name, args.suffix)
            list_of_vcf_files = generate_custom_vcf_file_list(sorted(list_of_files), logger)
            keep_logging('\nNumber of final variant call vcf files: %s\n' % len(list_of_vcf_files), '\nNumber of final variant call vcf files: %s\n' % len(list_of_vcf_files), logger, 'info')
            #keep_logging('Make sure the number of final variant call vcf files looks correct', 'Make sure the number of final variant call vcf files looks correct', logger, 'info')
            #keep_logging('Running core_prep on these files...', 'Running core_prep on these files...', logger, 'info')
            with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
                for file in list_of_files:
                    #file = file.splitlines()
                    print file
                    depth = "grep -vE '^sample|Total' %s | awk -F'\t' '{print $3}'"
                    proc = subprocess.Popen(["grep -vE '^sample|Total' %s | awk -F'\t' '{print $3}'" % file.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_depth_of_coverage.sample_summary')], stdout=subprocess.PIPE, shell=True)
                    (out2, err2) = proc.communicate()
                    print file
                    cov_depth = int(float(out2.strip()))
                    if float(out2.strip()) > float(ConfigSectionMap(filter_criteria, Config)['dp']):
                        out_fp.write(os.path.basename(file) + '\n')
            out_fp.close()
        else:
            keep_logging('\nChecking if all the variant calling results exists in %s\n' % core_temp_dir, '\nChecking if all the variant calling results exists in %s\n' % core_temp_dir, logger, 'info')
            #call("ls -1a %s/*.vcf_no_proximate_snp.vcf > %s/vcf_filenames" % (core_temp_dir, core_temp_dir), logger)
            try:
                list_cmd = "ls -1a %s/*.vcf_no_proximate_snp.vcf" % core_temp_dir
                list_of_files = subprocess.check_output(list_cmd, shell=True)

                keep_logging('Number of final variant call vcf files: %s\n' % len(list_of_files.splitlines()), 'Number of final variant call vcf files: %s\n' % len(list_of_files.splitlines()), logger, 'info')
                #keep_logging('Make sure the number of final variant call vcf files looks correct', 'Make sure the number of final variant call vcf files looks correct', logger, 'info')
                #keep_logging('Running core_prep on these files...', 'Running core_prep on these files...', logger, 'info')
                with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
                    for file in list_of_files.splitlines():
                        depth = "grep -vE '^sample|Total' %s | awk -F'\t' '{print $3}'"
                        proc = subprocess.Popen(["grep -vE '^sample|Total' %s | awk -F'\t' '{print $3}'" % file.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_depth_of_coverage.sample_summary')], stdout=subprocess.PIPE, shell=True)
                        (out2, err2) = proc.communicate()
                        print file
                        cov_depth = int(float(out2.strip()))
                        if float(out2.strip()) > float(ConfigSectionMap(filter_criteria, Config)['dp']):
                            out_fp.write(os.path.basename(file)+'\n')
                out_fp.close()
            except:
                keep_logging('Error: The variant calling results were not found in %s. Please check if variant calling step finished properly without any errors. '
                             'This can be done by checking if all the variant call results folder contains final variant call vcf file: *.vcf_no_proximate_snp.vcf file' % core_temp_dir, 'Error: The variant calling results were not found in %s. Please check if variant calling step finished properly without any errors. '
                                                                                                                                                                                              'This can be done by checking if all the variant call results folder contains final variant call vcf file: *.vcf_no_proximate_snp.vcf file' % core_temp_dir, logger, 'exception')
                exit()





        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
        core_prep_pipeline_cmd = run_core_prep_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger, config_file)
        core_All_cmds.append(core_prep_pipeline_cmd)
        time_taken = datetime.now() - start_time_2
        keep_logging('Core Prep Logs will be recorded in file with extension log.txt in %s\n' % core_prep_logs_folder, 'Core Prep Logs will be recorded in file with extension log.txt in %s\n' % core_prep_logs_folder, logger, 'info')
        #keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')

        """ Generate Core Variants from core_prep intermediate files """
        core_logs_folder = logs_folder + "/core"
        make_sure_path_exists(core_logs_folder)
        logger = generate_logger(core_logs_folder, args.analysis_name, log_unique_time)
        #keep_logging('START: Extract core snps and generate diagnostic plots','START: Extract core snps and generate diagnostic plots', logger, 'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, core_logs_folder, log_unique_time, args.analysis_name),
             logger)
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        core_results_dir = args.output_folder + "/%s_core_results/" % log_unique_time
        make_sure_path_exists(core_results_dir)

        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)[
            'ref_name']

        """ If Phaster Summary file doesn't exist in reference genome folder """
        if not os.path.isfile("%s/summary.txt" % os.path.dirname(reference)):
            if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
                keep_logging('Preparing Functional class filters\n', 'Preparing Functional class filters\n', logger,
                             'info')
                if ConfigSectionMap("functional_filters", Config)['find_phage_region'] == "yes":
                    # Submit Phaster jobs to find ProPhage region in reference genome.
                    run_phaster(reference, os.path.dirname(reference), logger, Config)



        # Parse Phaster results file to extract phage region.
        if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
            keep_logging('Parsing Functional class filters\n', 'Parsing Functional class filters\n', logger,
                         'info')
            functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % core_temp_dir
            f1 = open(functional_class_filter_positions, 'w+')
            if ConfigSectionMap("functional_filters", Config)['find_phage_region'] == "yes":
                phage_region_positions = parse_phaster(reference, core_temp_dir, logger, Config)
                with open(phage_region_positions, 'rU') as fp:
                    for line in fp:
                        f1.write(line)
                fp.close()
            if ConfigSectionMap("functional_filters", Config)['find_repetitive_region'] == "yes":
                # Find repeat regions in reference genome
                repeat_region_positions = nucmer_repeat(reference, core_temp_dir, logger, Config)
                with open(repeat_region_positions, 'rU') as fp:
                    for line in fp:
                        f1.write(line)
                fp.close()
            if ConfigSectionMap("functional_filters", Config)['mask_region'] == "yes":
                # Mask custom region/Positions
                if ConfigSectionMap("functional_filters", Config)['mask_file']:
                    mask_file = ConfigSectionMap("functional_filters", Config)['mask_file']
                    mask_extension = os.path.splitext(mask_file)[1]
                    if mask_extension == ".bed":
                        mask_positions_file = mask_regions(mask_file, core_temp_dir, logger, Config)
                        keep_logging(
                            'Mask positions in this file %s will be filtered out' % mask_positions_file,
                            'Mask positions in this file %s will be filtered out' % mask_positions_file,
                            logger, 'info')
                    else:
                        # mask_positions_file = mask_file
                        os.system("cp %s %s/mask_positions.txt" % (mask_file, core_temp_dir))
                        mask_positions_file = "%s/mask_positions.txt" % core_temp_dir
                        keep_logging(
                            'Mask positions in this file %s will be filtered out' % mask_positions_file,
                            'Mask positions in this file %s will be filtered out' % mask_positions_file,
                            logger, 'info')
                    with open(mask_positions_file, 'rU') as fp:
                        for line in fp:
                            f1.write(line)
                    fp.close()
            f1.close()

        copy_phage_results = "cp %s/summary.txt %s/phage_region_positions.txt %s/detail.txt" % (os.path.dirname(reference), os.path.dirname(reference), os.path.dirname(reference))
        call(copy_phage_results, logger)

        core_pipeline_cmd = run_core_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger,
                          core_results_dir, config_file)
        core_All_cmds.append(core_pipeline_cmd)
        time_taken = datetime.now() - start_time_2
        keep_logging('Core Step Logs will be recorded in file with extension log.txt in %s\n' % core_logs_folder,
                     'Core Step Logs will be recorded in file with extension log.txt in %s\n' % core_logs_folder, logger, 'info')
        #keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger,'info')

        """ Generate Reports and organize core results folder """
        report_logs_folder = logs_folder + "/report"
        make_sure_path_exists(report_logs_folder)
        logger = generate_logger(report_logs_folder, args.analysis_name, log_unique_time)
        #keep_logging('START: Generating alignment and variant calling report','START: Generating alignment and variant calling report', logger, 'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, report_logs_folder, log_unique_time, args.analysis_name),
             logger)
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        # core_results_dir = args.output_folder + "/%s_core_results/" % log_unique_time
        # make_sure_path_exists(core_results_dir)
        proc = subprocess.Popen(["ls -1ad %s/*_core_results | tail -n1" % args.output_folder], stdout=subprocess.PIPE,
                                shell=True)
        (out2, err2) = proc.communicate()
        core_results_dir = out2.strip()
        #keep_logging('Moving final results and reports to %s' % core_results_dir,'Moving final results and reports to %s' % core_results_dir, logger, 'info')

        list_of_label_files = glob.glob("%s/*_label" % core_temp_dir)
        list_of_vcf_files = []
        with open("%s/vcf_filenames" % core_temp_dir, 'r') as out_fp:
            for line in out_fp:
                list_of_vcf_files.append(line)
        for i in list_of_label_files:
            if os.stat(i).st_size == 0:
                keep_logging('The file {} is empty. Please rerun core_prep step again.\n'.format(i),
                             'The file {} is empty. Please rerun core_prep step again.\n'.format(i), logger,
                             'exception')
                exit()
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)[
            'ref_name']
        run_report_analysis_cmd = run_report_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger,
                            core_results_dir, config_file)
        time_taken = datetime.now() - start_time_2
        keep_logging('Report Step Logs will be recorded in file with extension log.txt in %s\n' % report_logs_folder,
                     'Report Step Logs will be recorded in file with extension log.txt in %s\n' % report_logs_folder, logger, 'info')
        #keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger,'info')

        core_All_cmds.append(run_report_analysis_cmd)

        """ Generate Trees """
        tree_logs_folder = logs_folder + "/tree"
        make_sure_path_exists(tree_logs_folder)
        logger = generate_logger(tree_logs_folder, args.analysis_name, log_unique_time)
        #keep_logging('START: Generating FastTree and RAxML trees', 'START: Generating FastTree and RAxML trees', logger,'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, tree_logs_folder, log_unique_time, args.analysis_name),
             logger)
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        # core_results_dir = args.output_folder + "/%s_core_results/" % log_unique_time
        proc = subprocess.Popen(["ls -1ad %s/*_core_results | tail -n1" % args.output_folder], stdout=subprocess.PIPE,
                                shell=True)
        (out2, err2) = proc.communicate()
        core_results_dir = out2.strip()
        #keep_logging('Generating RaxML and fasttree trees in %s/trees' % core_results_dir,'Generating RaxML and fasttree trees in %s/trees' % core_results_dir, logger, 'info')
        print core_results_dir
        list_of_label_files = glob.glob("%s/*_label" % core_temp_dir)
        list_of_vcf_files = []
        with open("%s/vcf_filenames" % core_temp_dir, 'r') as out_fp:
            for line in out_fp:
                list_of_vcf_files.append(line)
        for i in list_of_label_files:
            if os.stat(i).st_size == 0:
                keep_logging('The file {} is empty. Please rerun core_prep step again.\n'.format(i),
                             'The file {} is empty. Please rerun core_prep step again.\n'.format(i), logger,
                             'exception')
                exit()
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)[
            'ref_name']
        run_tree_analysis_cmd = run_tree_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time,
                                                  args.cluster, logger, core_results_dir, config_file)
        time_taken = datetime.now() - start_time_2
        keep_logging('Tree Step Logs will be recorded in file with extension log.txt in %s\n' % tree_logs_folder,
                     'Tree Step Logs will be recorded in file with extension log.txt in %s\n' % tree_logs_folder, logger, 'info')
        #keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger,'info')

        core_All_cmds.append(run_tree_analysis_cmd)

        combine_job_name = core_temp_dir + "/" + log_unique_time + "_" + args.analysis_name + "_core_All.pbs"

        Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                          % (ConfigSectionMap("scheduler", Config)['email'],
                             ConfigSectionMap("scheduler", Config)['notification'],
                             ConfigSectionMap("scheduler", Config)['large_resources'],
                             ConfigSectionMap("scheduler", Config)['queue'],
                             ConfigSectionMap("scheduler", Config)['flux_account'])

        with open(combine_job_name, 'w') as out:
            job_title = "#PBS -N %s_%s_core_All" % (log_unique_time, args.analysis_name)
            out.write(job_title + '\n')
            out.write(Pbs_model_lines + '\n')
            out.write(
                "#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
            out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
            out.write("cd %s" % core_temp_dir + '\n')
            for cmds in core_All_cmds:
                out.write(cmds + '\n')
        out.close()

        keep_logging('Running: %s\n' % combine_job_name,
                     'Running: %s\n' % combine_job_name,
                     logger, 'info')

        call("qsub %s" % combine_job_name, logger)


    elif "core_prep" in args.steps:
        """ Set Up Core Prep logs folder/logger object, cluster mode and copy config files to it"""
        core_prep_logs_folder = logs_folder + "/core_prep"
        make_sure_path_exists(core_prep_logs_folder)
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        logger = generate_logger(core_prep_logs_folder, args.analysis_name, log_unique_time)
        keep_logging('START: Copying variant calling results to %s and running core_prep step.' % core_temp_dir, 'START: Copying variant calling results to %s and running core_prep step.' % core_temp_dir, logger, 'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, core_prep_logs_folder, log_unique_time, args.analysis_name), logger)
        make_sure_path_exists(core_temp_dir)
        keep_logging('Copying vcf files to %s' % core_temp_dir, 'Copying vcf files to %s' % core_temp_dir, logger, 'info')
        cp_command = "cp %s/*/*_vcf_results/*_filter2_indel_final.vcf %s/*/*_vcf_results/*_aln_mpileup_raw.vcf %s/*/*_vcf_results/*_raw.vcf_5bp_indel_removed.vcf* %s/*/*_vcf_results/*filter2_final.vcf* %s/*/*_vcf_results/*vcf_no_proximate_snp.vcf* %s/*/*_vcf_results/*array %s/*/*unmapped.bed_positions %s" % (args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, core_temp_dir)
        # cp_command = "cp %s/*/*_filter2_indel_final.vcf %s/*/*_aln_mpileup_raw.vcf %s/*/*_raw.vcf_5bp_indel_removed.vcf* %s/*/*filter2_final.vcf* %s/*/*vcf_no_proximate_snp.vcf* %s/*/*array %s/*/*unmapped.bed_positions %s" % (
        # args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder,
        # args.output_folder, args.output_folder, core_temp_dir)
        call(cp_command, logger)

        """ Decompress zipped files in core temp folder"""
        list_of_gzipped_files = glob.glob("%s/*.gz" % core_temp_dir)
        keep_logging('Decompressing gzipped files in %s' % core_temp_dir, 'Decompressing gzipped files in %s' % core_temp_dir, logger, 'info')
        gzipped_command_list = []
        for i in list_of_gzipped_files:
            gzipped_command_list.append("gzip -df %s" % i)

        """ Set the number of cores variable for Parallel-local mode """
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command_list)(i) for i in gzipped_command_list)

        """ Generate a custom vcf file list to process for core prep step. If file is not provided, it will consider all the samples in output folder"""
        """ Perform Sanity check to confirm that variant calling results are present in core temp folder """
        if args.filenames:
            list_of_files = get_filenames(args.dir, args.type, args.filenames, args.analysis_name, args.suffix)
            list_of_vcf_files = generate_custom_vcf_file_list(sorted(list_of_files), logger)
            keep_logging('Number of final variant call vcf files: %s' % len(list_of_vcf_files), 'Number of final variant call vcf files: %s' % len(list_of_vcf_files), logger, 'info')
            keep_logging('Make sure the number of final variant call vcf files looks correct', 'Make sure the number of final variant call vcf files looks correct', logger, 'info')
            keep_logging('Running core_prep on these files...', 'Running core_prep on these files...', logger, 'info')
            with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
                for file in list_of_vcf_files:
                    out_fp.write(os.path.basename(file)+'\n')
            out_fp.close()
        else:
            keep_logging('Checking if all the variant calling results exists in %s' % core_temp_dir, 'Checking if all the variant calling results exists in %s' % core_temp_dir, logger, 'info')
            #call("ls -1a %s/*.vcf_no_proximate_snp.vcf > %s/vcf_filenames" % (core_temp_dir, core_temp_dir), logger)
            try:
                list_cmd = "ls -1a %s/*.vcf_no_proximate_snp.vcf" % core_temp_dir
                list_of_files = subprocess.check_output(list_cmd, shell=True)

                keep_logging('Number of final variant call vcf files: %s' % len(list_of_files.splitlines()), 'Number of final variant call vcf files: %s' % len(list_of_files.splitlines()), logger, 'info')
                keep_logging('Make sure the number of final variant call vcf files looks correct', 'Make sure the number of final variant call vcf files looks correct', logger, 'info')
                keep_logging('Running core_prep on these files...', 'Running core_prep on these files...', logger, 'info')
                with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
                    for file in list_of_files.splitlines():
                        out_fp.write(os.path.basename(file)+'\n')
                out_fp.close()
            except:
                keep_logging('Error: The variant calling results were not found in %s. Please check if variant calling step finished properly without any errors. '
                             'This can be done by checking if all the variant call results folder contains final variant call vcf file: *.vcf_no_proximate_snp.vcf file' % core_temp_dir, 'Error: The variant calling results were not found in %s. Please check if variant calling step finished properly without any errors. '
                                                                                                                                                                                              'This can be done by checking if all the variant call results folder contains final variant call vcf file: *.vcf_no_proximate_snp.vcf file' % core_temp_dir, logger, 'exception')
                exit()



            # list_cmd = "ls -1a %s/*.vcf_no_proximate_snp.vcf" % core_temp_dir
            # list_of_files = subprocess.check_output(list_cmd, shell=True)
            # print list_of_files
            # with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
            #     for file in list_of_files.splitlines():
            #         out_fp.write(os.path.basename(file)+'\n')
            # out_fp.close()
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
        core_prep_pipeline_cmd = run_core_prep_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger, config_file)
        time_taken = datetime.now() - start_time_2
        keep_logging('Logs were recorded in file with extension log.txt in %s' % core_prep_logs_folder, 'Logs were recorded in file with extension log.txt in %s' % core_prep_logs_folder, logger, 'info')
        keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')

    elif "core" in args.steps:
        """ Generate Core Variants from core_prep intermediate files """
        core_logs_folder = logs_folder + "/core"
        make_sure_path_exists(core_logs_folder)
        logger = generate_logger(core_logs_folder, args.analysis_name, log_unique_time)
        keep_logging('START: Extract core snps and generate diagnostic plots', 'START: Extract core snps and generate diagnostic plots', logger, 'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, core_logs_folder, log_unique_time, args.analysis_name), logger)
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        core_results_dir = args.output_folder + "/%s_core_results/" % log_unique_time
        #make_sure_path_exists(core_results_dir)
        if args.filenames:
            list_of_files = get_filenames(args.dir, args.type, args.filenames, args.analysis_name, args.suffix)
            list_of_vcf_files = generate_custom_vcf_file_list(sorted(list_of_files), logger)
            list_of_label_files = []

            for i in list_of_vcf_files:
                list_of_label_files.append(i + '_positions_label')
            if len(list_of_label_files) == len(list_of_vcf_files):
                for i in list_of_label_files:
                    if os.stat("%s/%s" % (core_temp_dir, i)).st_size == 0:
                        keep_logging('The file {} is empty. Please rerun core_prep step again.\n'.format(i), 'The file {} is empty. Please rerun core_prep step again.\n'.format(i), logger, 'exception')
                        exit()
            else:
                keep_logging('Problem in core_prep results. Rerun the core_prep step\n', 'Problem in core_prep results. Rerun the core_prep step\n', logger, 'exception')
                exit()

            with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
                for file in list_of_vcf_files:
                    out_fp.write(os.path.basename(file)+'\n')
            out_fp.close()
        else:
            list_of_label_files = glob.glob("%s/*_no_proximate_snp.vcf_positions_label" % core_temp_dir)
            list_of_vcf_files = glob.glob("%s/*_filter2_final.vcf_no_proximate_snp.vcf" % core_temp_dir)
            #print sorted(list_of_vcf_files)
            if len(list_of_label_files) == len(list_of_vcf_files):
                for i in list_of_label_files:
                    if os.stat(i).st_size == 0:
                        keep_logging('The file {} is empty. Please rerun core_prep step again.\n'.format(i), 'The file {} is empty. Please rerun core_prep step again.\n'.format(i), logger, 'exception')
                        exit()
            else:
                keep_logging('Problem in core_prep results. Rerun the core_prep step\n', 'Problem in core_prep results. Rerun the core_prep step\n', logger, 'exception')
                exit()
            with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
                for file in sorted(list_of_vcf_files):
                    out_fp.write(os.path.basename(file)+'\n')
            out_fp.close()
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']

        # Parse Phaster results file to extract phage region.
        if ConfigSectionMap("functional_filters", Config)['apply_functional_filters'] == "yes":
            keep_logging('Preparing Functional class filters\n', 'Preparing Functional class filters\n', logger,
                         'info')
            functional_class_filter_positions = "%s/Functional_class_filter_positions.txt" % core_temp_dir
            f1 = open(functional_class_filter_positions, 'w+')
            if ConfigSectionMap("functional_filters", Config)['find_phage_region'] == "yes":
                phage_region_positions = parse_phaster(reference, core_temp_dir, logger, Config)
                with open(phage_region_positions, 'rU') as fp:
                    for line in fp:
                        f1.write(line)
                fp.close()
            if ConfigSectionMap("functional_filters", Config)['find_repetitive_region'] == "yes":
                # Find repeat regions in reference genome
                repeat_region_positions = nucmer_repeat(reference, core_temp_dir, logger, Config)
                with open(repeat_region_positions, 'rU') as fp:
                    for line in fp:
                        f1.write(line)
                fp.close()
            if ConfigSectionMap("functional_filters", Config)['mask_region'] == "yes":
                # Mask custom region/Positions
                if ConfigSectionMap("functional_filters", Config)['mask_file']:
                    mask_file = ConfigSectionMap("functional_filters", Config)['mask_file']
                    mask_extension = os.path.splitext(mask_file)[1]
                    if mask_extension == ".bed":
                        mask_positions_file = mask_regions(mask_file, core_temp_dir, logger, Config)
                        keep_logging(
                            'Mask positions in this file %s will be filtered out' % mask_positions_file,
                            'Mask positions in this file %s will be filtered out' % mask_positions_file,
                            logger, 'info')
                    else:
                        # mask_positions_file = mask_file
                        os.system("cp %s %s/mask_positions.txt" % (mask_file, core_temp_dir))
                        mask_positions_file = "%s/mask_positions.txt" % core_temp_dir
                        keep_logging(
                            'Mask positions in this file %s will be filtered out' % mask_positions_file,
                            'Mask positions in this file %s will be filtered out' % mask_positions_file,
                            logger, 'info')
                    with open(mask_positions_file, 'rU') as fp:
                        for line in fp:
                            f1.write(line)
                    fp.close()
            f1.close()

        core_pipeline_cmd = run_core_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger, core_results_dir, config_file)
        time_taken = datetime.now() - start_time_2
        keep_logging('Logs were recorded in file with extension log.txt in %s' % core_logs_folder, 'Logs were recorded in file with extension log.txt in %s' % core_logs_folder, logger, 'info')
        keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')

    elif "report" in args.steps:
        report_logs_folder = logs_folder + "/report"
        make_sure_path_exists(report_logs_folder)
        logger = generate_logger(report_logs_folder, args.analysis_name, log_unique_time)
        keep_logging('START: Generating alignment and variant calling report', 'START: Generating alignment and variant calling report', logger, 'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, report_logs_folder, log_unique_time, args.analysis_name), logger)
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        # core_results_dir = args.output_folder + "/%s_core_results/" % log_unique_time
        # make_sure_path_exists(core_results_dir)
        proc = subprocess.Popen(["ls -1ad %s/*_core_results | tail -n1" % args.output_folder], stdout=subprocess.PIPE, shell=True)
        (out2, err2) = proc.communicate()
        core_results_dir = out2.strip()
        keep_logging('Moving final results and reports to %s' % core_results_dir,
                     'Moving final results and reports to %s' % core_results_dir, logger, 'info')


        list_of_label_files = glob.glob("%s/*_label" % core_temp_dir)
        list_of_vcf_files = []
        with open("%s/vcf_filenames" % core_temp_dir, 'r') as out_fp:
            for line in out_fp:
                list_of_vcf_files.append(line)
        for i in list_of_label_files:
            if os.stat(i).st_size == 0:
                keep_logging('The file {} is empty. Please rerun core_prep step again.\n'.format(i), 'The file {} is empty. Please rerun core_prep step again.\n'.format(i), logger, 'exception')
                exit()
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
        run_report_analysis_cmd = run_report_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger, core_results_dir, config_file)
        time_taken = datetime.now() - start_time_2
        keep_logging('Logs were recorded in file with extension log.txt in %s' % report_logs_folder, 'Logs were recorded in file with extension log.txt in %s' % report_logs_folder, logger, 'info')
        keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')

    elif "tree" in args.steps:
        tree_logs_folder = logs_folder + "/tree"
        make_sure_path_exists(tree_logs_folder)
        logger = generate_logger(tree_logs_folder, args.analysis_name, log_unique_time)
        keep_logging('START: Generating FastTree and RAxML trees', 'START: Generating FastTree and RAxML trees', logger, 'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, tree_logs_folder, log_unique_time, args.analysis_name), logger)
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        #core_results_dir = args.output_folder + "/%s_core_results/" % log_unique_time
        proc = subprocess.Popen(["ls -1ad %s/*_core_results | tail -n1" % args.output_folder], stdout=subprocess.PIPE,
                                shell=True)
        (out2, err2) = proc.communicate()
        core_results_dir = out2.strip()
        keep_logging('Generating RaxML and fasttree trees in %s/trees' % core_results_dir,
                     'Generating RaxML and fasttree trees in %s/trees' % core_results_dir, logger, 'info')
        print core_results_dir
        list_of_label_files = glob.glob("%s/*_label" % core_temp_dir)
        list_of_vcf_files = []
        with open("%s/vcf_filenames" % core_temp_dir, 'r') as out_fp:
            for line in out_fp:
                list_of_vcf_files.append(line)
        for i in list_of_label_files:
            if os.stat(i).st_size == 0:
                keep_logging('The file {} is empty. Please rerun core_prep step again.\n'.format(i), 'The file {} is empty. Please rerun core_prep step again.\n'.format(i), logger, 'exception')
                exit()
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
        run_tree_analysis_cmd = run_tree_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger, core_results_dir, config_file)
        time_taken = datetime.now() - start_time_2
        keep_logging('Logs were recorded in file with extension log.txt in %s' % tree_logs_folder, 'Logs were recorded in file with extension log.txt in %s' % tree_logs_folder, logger, 'info')
        keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')

    else:
        logger = generate_logger(logs_folder, args.analysis_name, log_unique_time)
        keep_logging('Please provide argument -steps to run pipeline', 'Please provide argument -steps to run pipeline', logger, 'info')
        time_taken = datetime.now() - start_time_2
        keep_logging('Logs were recorded in file with extension log.txt in %s' % logs_folder, 'Logs were recorded in file with extension log.txt in %s' % logs_folder, logger, 'info')
        keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')


