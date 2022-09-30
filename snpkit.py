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
from modules.snpeff import *

""" Command Line Argument Parsing """
def parser():
    parser = argparse.ArgumentParser(description='\nSNPKIT - A workflow for Microbial Variant Calling, Recombination detection and Phylogenetic tree reconstruction.\n', formatter_class=RawTextHelpFormatter)
    input = parser.add_argument_group('INPUT')
    output = parser.add_argument_group('OUTPUT')
    resources = parser.add_argument_group('RESOURCES')
    subsample = parser.add_argument_group('SUBSAMPLE')
    advanced = parser.add_argument_group('ADVANCED')
    optional = parser.add_argument_group('Optional')
    development = parser.add_argument_group('Development Phase')
    input.add_argument('-type', action='store', dest="type", help='Type of illumina reads. Options: SE for single-end, PE for paired-end', required=True)
    input.add_argument('-readsdir', action='store', dest="dir", help='path to input sequencing reads data folder.', required=True)
    input.add_argument('-index', action='store', dest="index", help='Reference genome index name (prefix) as described in -config file.', required=True)
    input.add_argument('-steps', action='store', dest="steps", help='Run this part of snpkit. Options: All, core_All\n'
    'All: run first part of snpkit - trimming, mapping, calling variants;\n'
    'core_All: run second part of snpkit - filter variants, generate core/non-core multi-fasta alignments, SNP/Indel Matrices.')
    input.add_argument('-analysis', action='store', dest="prefix", help='prefix for output files.', required=True)
    output.add_argument('-outdir', action='store', dest="output", help='output directory path ending with output directory name.', required=True)
    
    resources.add_argument('-cluster', action='store', dest='cluster', help='run snpkit in local or cluster mode. Default: local\n')
    resources.add_argument('-scheduler', action='store', dest="scheduler",
                          help='Type of HPC job scheduler. Supports PBS, SLURM')
    optional.add_argument('-config', action='store', dest="config", help='Path to YAML config file. Default: snpkit config - snpkit/config', required=False)
    optional.add_argument('-suffix', action='store', dest="suffix", help='Custom fastq reads suffix. Default: fastq.gz. supports *.fastq, *.fastq.gz, *.fq.gz, *.fq; ', required=False)
    optional.add_argument('-filenames', action='store', dest="filenames", help='Run snpkit on a subset of files in -readsdir folder. The file should contain single-end filename per line.', required=False)
    
    optional.add_argument('-clean', action="store_true", help='Delete all intermediate files. Default: yes Options: yes, no')
    development.add_argument('-extract_unmapped', action='store', dest="extract_unmapped", help='Extract unmapped reads, assemble it and detect AMR genes using ariba')
    optional.add_argument('-gubbins', action='store', dest="gubbins", help='Run Gubbins. Options: yes,no. Default: no')
    development.add_argument('-outgroup', action='store', dest="outgroup", help='Outgroup sample name. Alpha testing version. Not recommended.')
    subsample.add_argument('-downsample', action='store', dest="downsample", help='Subsample reads to a default depth of 100X or user specified -coverage_depth')
    subsample.add_argument('-coverage_depth', action='store', dest="coverage_depth",
                          help='Downsample reads to this depth')
    subsample.add_argument('-genomesize', action='store', dest="genomesize",
                          help='Genome size to calculate raw coverage')
    
    optional.add_argument('-dryrun', action='store_true', dest="dryrun",
                          help='Perform a trial run without running any jobs.',
                          required=False)
    optional.add_argument('-mask', action='store_true', dest="mask",
                          help='Mask Gubbins detected recombinant region and run Iqtree on masked alignment')
    optional.add_argument('-clip', action="store_true",
                          help='Filter SAM file for soft and hard clipped alignments. Default: OFF')
    return parser

""" Sanity checks and maintenance methods """
def file_exists(path1):
    """Checks if the file path exists.
        Args:
            path: file path

        Output:
            True, if the file path exists. False and exists if the file path is not found.
    """
    if not os.path.isfile(path1):
        file_basename = os.path.basename(path1)
        keep_logging('The file {} does not exists. Please provide another file with full path or check the files path.\n'.format(file_basename), 'The input file {} does not exists. Please provide another file or check the files path.\n'.format(file_basename), logger, 'exception')
        exit()

def make_sure_path_exists(out_path):
    """Checks the directory path exists. If not, creates a new directory.
        Args:
            path: Path to Directory

        Output:
            True, if the directory exists or if not, a new directory is created.
    """
    try:
        os.makedirs(out_path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            keep_logging('Errors in output folder path! please change the output path or analysis name.', 'Errors in output folder path! please change the output path or analysis name', logger, 'exception')
            exit()

def get_filenames(dir, type, filenames, analysis, suffix):
    """Get a list of file with specific suffix from a directory
        Args:
            dir: directory to the files with suffix
            type: type of reads. PE/SE
            filename: list of files to fetch from directory
            analysis: unique analysis name.
            suffix: fastq suffix.

        Output:
            True, if the file path exists. False and exists if the file path is not found.
    """
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

def get_scheduler_directive(scheduler, Config):
    """ Generate Cluster Directive lines for a scheduler provided with args.scheduler"""
    # Scheduler Changes here; current changes
    if scheduler and scheduler == "SLURM":
        script_Directive = "#SBATCH"
        job_name_flag = "--job-name="
        scheduler_directives = "#SBATCH --mail-user=%s\n#SBATCH --mail-type=%s\n#SBATCH --export=ALL\n#SBATCH --partition=%s\n#SBATCH --account=%s\n#SBATCH %s\n" \
                          % (ConfigSectionMap("slurm", Config)['email'],
                             ConfigSectionMap("slurm", Config)['notification'],
                             ConfigSectionMap("slurm", Config)['partition'],
                             ConfigSectionMap("slurm", Config)['flux_account'],
                             ConfigSectionMap("slurm", Config)['resources'])
    elif scheduler and scheduler == "PBS":
        script_Directive = "#PBS"
        job_name_flag = "-N"
        scheduler_directives = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                          % (ConfigSectionMap("scheduler", Config)['email'],
                             ConfigSectionMap("scheduler", Config)['notification'],
                             ConfigSectionMap("scheduler", Config)['resources'],
                             ConfigSectionMap("scheduler", Config)['queue'],
                             ConfigSectionMap("scheduler", Config)['flux_account'])
    else:
        script_Directive = "#SBATCH"
        job_name_flag = "--job-name="
        scheduler_directives = "#SBATCH --mail-user=%s\n#SBATCH --mail-type=%s\n#SBATCH --export=ALL\n#SBATCH --partition=%s\n#SBATCH --account=%s\n#SBATCH %s\n" \
                               % (ConfigSectionMap("slurm", Config)['email'],
                                  ConfigSectionMap("slurm", Config)['notification'],
                                  ConfigSectionMap("slurm", Config)['partition'],
                                  ConfigSectionMap("slurm", Config)['flux_account'],
                                  ConfigSectionMap("slurm", Config)['resources'])
    return scheduler_directives, script_Directive, job_name_flag

""" Methods to generate jobs for various pipeline tasks """
def create_varcall_jobs(filenames_array, type, output_folder, reference, steps, config_file, logger):
    """Takes a list of files and other arguments, generate variant calling jobs.
        Args:
            filenames_array: list of fastq file names
            output_folder: output directory to save the results/jobs
            type: type of reads. PE/SE
            reference: Reference Genome name

        Output:
            List of variant calling jobs to run/submit
    """
    jobs_temp_dir = "%s/temp_jobs" % output_folder
    make_sure_path_exists(jobs_temp_dir)
    keep_logging('Generating cluster jobs in temporary directory %s' % jobs_temp_dir, 'Generating cluster jobs in temporary directory %s' % jobs_temp_dir, logger, 'exception')

    scheduler_directives, script_Directive, job_name_flag = get_scheduler_directive(args.scheduler, Config)

    for file in filenames_array:
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base or "R1.fastq.gz" in filename_base or "1_combine.fastq.gz" in filename_base or "1_sequence.fastq.gz" in filename_base or "_forward.fastq.gz" in filename_base or "R1_001.fastq.gz" in filename_base or "_1.fastq.gz" in filename_base or ".1.fastq.gz" in filename_base or "_R1.fastq.gz" in filename_base or "_L001_R1_001.fastq.gz" in filename_base:
            # Forward reads file name and get analysis name from its name
            first_file = file
            # Get the name of reverse reads files
            if "R1_001_final.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "R1_001.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "_R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
                first_part_split = filename_base.split('_R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1.fastq.gz", "R2.fastq.gz")
                first_part_split = filename_base.split('R1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "1_combine.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_combine.fastq.gz", "2_combine.fastq.gz")
                first_part_split = filename_base.split('1_combine.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "1_sequence.fastq.gz" in filename_base:
                second_part = filename_base.replace("1_sequence.fastq.gz", "2_sequence.fastq.gz")
                first_part_split = filename_base.split('1_sequence.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "_forward.fastq.gz" in filename_base:
                second_part = filename_base.replace("_forward.fastq.gz", "_reverse.fastq.gz")
                first_part_split = filename_base.split('_forward.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "R1_001.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001.fastq.gz", "R2_001.fastq.gz")
                first_part_split = filename_base.split('R1_001.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif "_1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_1.fastq.gz", "_2.fastq.gz")
                first_part_split = filename_base.split('_1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)
            elif ".1.fastq.gz" in filename_base:
                second_part = filename_base.replace(".1.fastq.gz", ".2.fastq.gz")
                first_part_split = filename_base.split('.1.fastq.gz')
                first_part = first_part_split[0].replace('_L001', '')
                first_part = re.sub("_S[0-9].*_", "", first_part)

            """ Have a standard filename preparation step"""
            # else:
            #     print "Using Standard second file naming convention"
            #     second_part = filename_base.replace("_R1_", "_R2_")
            #     first_part_split = filename_base.split('_R1.fastq.gz')
            #     first_part = first_part_split[0].replace('_L001', '')
            #     first_part = re.sub("_S.*_", "", first_part)
            second_file = args.dir + "/" + second_part
            if args.scheduler == "SLURM":
                job_name = jobs_temp_dir + "/" + first_part + ".sbat"
            else:
                job_name = jobs_temp_dir + "/" + first_part + ".pbs"
            if not steps:
                steps == "All"
            if type == "SE":
                if args.clean:
                    command = "python %s/pipeline.py -PE1 %s -o %s/%s -analysis %s -index %s -type SE -config %s -steps %s -clean" % (
                    os.path.dirname(os.path.abspath(__file__)), first_file, output_folder, first_part, first_part,
                    reference, config_file, steps)
                else:
                    command = "python %s/pipeline.py -PE1 %s -o %s/%s -analysis %s -index %s -type SE -config %s -steps %s -clean" % (os.path.dirname(os.path.abspath(__file__)), first_file, output_folder, first_part, first_part, reference, config_file, steps)
                if args.genomesize:
                    command = command + " -genome_size %s" % args.genomesize
            else:
                if args.clean:
                    command = "python %s/pipeline.py -PE1 %s -PE2 %s -o %s/%s -analysis %s -index %s -type PE -config %s -steps %s -clean" % (
                    os.path.dirname(os.path.abspath(__file__)), first_file, second_file, output_folder, first_part,
                    first_part, reference, config_file, steps)
                else:
                    command = "python %s/pipeline.py -PE1 %s -PE2 %s -o %s/%s -analysis %s -index %s -type PE -config %s -steps %s -clean" % (os.path.dirname(os.path.abspath(__file__)), first_file, second_file, output_folder, first_part, first_part, reference, config_file, steps)
                if args.genomesize:
                    command = command + " -genome_size %s" % args.genomesize

            # # Adding Downsampling support 2019-06-20
            if args.downsample == "yes":
                if args.coverage_depth:
                    depth = args.coverage_depth
                else:
                    depth = 100

                command = command + " -downsample yes -coverage_depth %s" % depth

            # Adding samclip feature July 2020
            if args.clip:
               command = command + " -clip"
            with open(job_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, first_part)
                out.write("#!/bin/sh" + '\n')
                out.write(job_title+'\n')
                out.write(scheduler_directives+'\n')
                out.write("cd %s/temp_jobs" % output_folder + '\n')
                out.write(command+'\n')
        elif "R2_001_final.fastq.gz" in filename_base or "R2.fastq.gz" in filename_base or "2_combine.fastq.gz" in filename_base or "2_sequence.fastq.gz" in filename_base or "_reverse.fastq.gz" in filename_base or "R2_001.fastq.gz" in filename_base or "_2.fastq.gz" in filename_base or ".2.fastq.gz" in filename_base or "_R2.fastq.gz" in filename_base:
            continue
        else:
            keep_logging('Error while generating cluster jobs. Make sure the fastq filenames ends with one of these suffix: R1_001_final.fastq.gz, R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, R1_001.fastq.gz, _1.fastq.gz, .1.fastq.gz, _R1.fastq.gz', 'Error while generating cluster jobs. Make sure the fastq filenames ends with one of these suffix: R1_001_final.fastq.gz, R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, R1_001.fastq.gz, _1.fastq.gz, .1.fastq.gz, _R1.fastq.gz', logger, 'exception')
            print filename_base
            exit()
    if args.scheduler == "SLURM":
        list_of_jobs = glob.glob("%s/*.sbat" % jobs_temp_dir)
    else:
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
    call("%s" % command, logger)
    done = "Command Run completed: %s" % command
    return done

def run_varcall_jobs(list_of_jobs, cluster, log_unique_time, analysis_name, output_folder, logger):
    command_list = ""
    command_list_qsub = []
    for job in list_of_jobs:
        command_list = command_list + "bash %s\n" % job
        command_list_qsub.append(job)

    job_id_array = []
    if cluster == "cluster":
        keep_logging('Running Jobs in cluster mode', 'Running Jobs in cluster mode', logger, 'info')
        for job in command_list_qsub:
            keep_logging('Submitting Job: %s' % job, 'Submitting Job: %s' % job, logger, 'info')
            if args.scheduler == "SLURM":
                #call("sbatch %s" % job, logger)
                keep_logging("sbatch %s" % job, "sbatch %s" % job, logger, 'info')
                if not args.dryrun:
                    proc = subprocess.Popen(["sbatch %s" % job], stdout=subprocess.PIPE, shell=True)
                    (out, err) = proc.communicate()
                    # job_id_array.append(out.split(' ')[3].strip())
                    # job_id_array.append('123xx')
                    # break
            elif args.scheduler == "PBS":
                #call("qsub %s" % job, logger)
                keep_logging("qsub %s" % job, "qsub %s" % job, logger, 'info')
                if not args.dryrun:
                    proc = subprocess.Popen(["qsub %s" % job], stdout=subprocess.PIPE, shell=True)
                    (out, err) = proc.communicate()

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
            if not args.dryrun:
                call("bash %s" % job, logger)
    return job_id_array

""" Pipeline individual task methods """
def run_core_prep_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger, config_file):
    file_exists(reference)

    core_prep_pipeline = "python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 1 -jobrun %s -config %s -scheduler %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, config_file, args.scheduler)

    if args.scheduler == "SLURM":
        job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".sbat"
    else:
        job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"

    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                      % (ConfigSectionMap("scheduler", Config)['email'],
                         ConfigSectionMap("scheduler", Config)['notification'],
                         ConfigSectionMap("scheduler", Config)['resources'],
                         ConfigSectionMap("scheduler", Config)['queue'],
                         ConfigSectionMap("scheduler", Config)['flux_account'])

    if args.outgroup:
        core_prep_pipeline = core_prep_pipeline + " -outgroup %s" % args.outgroup
    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        core_prep_pipeline = "python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 1 -jobrun %s -config %s -scheduler %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, config_file, args.scheduler)
    return core_prep_pipeline

def run_core_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger, core_results_dir, config_file):
    file_exists(reference)
    core_pipeline = "python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 2 -jobrun %s -results_dir %s -config %s -scheduler %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file, args.scheduler)
    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"

    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n" \
                      % (ConfigSectionMap("scheduler", Config)['email'],
                         ConfigSectionMap("scheduler", Config)['notification'],
                         ConfigSectionMap("scheduler", Config)['large_resources'],
                         ConfigSectionMap("scheduler", Config)['queue'],
                         ConfigSectionMap("scheduler", Config)['flux_account'])

    if args.outgroup:
        core_pipeline = core_pipeline + " -outgroup %s" % args.outgroup
    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        pass
    elif cluster == "parallel-cluster":
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
    return core_pipeline

def run_report_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger, core_results_dir, config_file):
    file_exists(reference)

    core_pipeline = "python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 3 -jobrun %s -results_dir %s -config %s -scheduler %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file, args.scheduler)
    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"
    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l nodes=1:ppn=4,pmem=4000mb,walltime=92:00:00\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])

    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        pass
    elif cluster == "parallel-cluster":
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
    return core_pipeline

def run_tree_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger, core_results_dir, config_file):
    core_pipeline = "python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 4 -jobrun %s -results_dir %s -config %s -scheduler %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster, core_results_dir, config_file, args.scheduler)
    if args.gubbins == "yes":
        core_pipeline = core_pipeline + " -gubbins %s" % args.gubbins
    if args.outgroup:
        core_pipeline = core_pipeline + " -outgroup %s" % args.outgroup
    if args.mask:
        core_pipeline = core_pipeline + " -mask"


    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"
    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])

    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        pass
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
    return core_pipeline

def create_fai_index(reference, ref_fai_index):
    cmd = "%s %s %s" % (ConfigSectionMap("samtools", Config)['base_cmd'], ConfigSectionMap("samtools", Config)['faiindex'], reference)
    keep_logging('', cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Samtools FAI Indexing step. Exiting.', 'Error in Samtools FAI Indexing step. Exiting.', logger, 'exception')
        sys.exit(1)


    if not os.path.isfile(ref_fai_index):
        keep_logging('The reference fai index file {} was not created properly.\n Please try to create the samtools fai index files manually. \n'.format(ref_fai_index), 'The reference fai index file {} was not created properly.\n Please try to create the samtools fai index files manually. \n'.format(ref_fai_index), logger, 'exception')

def picard_seqdict(dict_name, reference):
    #dict_name = os.path.splitext(os.path.basename(reference_filename))[0] + ".dict"
    keep_logging('', 'Creating Sequence Dictionary using Picard.', logger, 'info')
    cmd = "%s CreateSequenceDictionary R=%s O=%s/%s" % (ConfigSectionMap("picard", Config)['base_cmd'], reference, ConfigSectionMap(args.index, Config)['ref_path'], dict_name)
    keep_logging('', cmd, logger, 'debug')
    try:
        call(cmd, logger)
    except sp.CalledProcessError:
        keep_logging('Error in Picard Sequence Dictionary creation step. Exiting.', 'Error in Picard Sequence Dictionary creation step. Exiting.', logger, 'exception')
        sys.exit(1)

def create_index(reference,ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4, ref_index_suffix5):
    aligner = ConfigSectionMap("pipeline", Config)['aligner']
    keep_logging('Creating Reference Genome Index for {} aligner.'.format(aligner), 'Creating Reference Genome Index for {} aligner'.format(aligner), logger, 'info')
    if aligner == "bwa":
        cmd = "%s %s %s" % (ConfigSectionMap("bwa", Config)['base_cmd'], ConfigSectionMap("bwa", Config)['index'], reference)
        keep_logging('', cmd, logger, 'debug')
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
                keep_logging('Error in {} Indexer. Exiting.'.format(aligner), 'Error in {} Indexer. Exiting.'.format(aligner), logger, 'exception')
                sys.exit(1)
        if not os.path.isfile(ref_index_suffix1):
            keep_logging('The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), 'The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), logger, 'exception')
    elif aligner == "bowtie":
        cmd = "%s %s %s" % ( ConfigSectionMap("bowtie", Config)['build_cmd'], reference, reference)
        keep_logging(cmd, cmd, logger, 'debug')
        try:
            call(cmd, logger)
        except sp.CalledProcessError:
                keep_logging('Error in {} Indexer. Exiting.'.format(aligner), 'Error in {} Indexer. Exiting.'.format(aligner), logger, 'exception')
                sys.exit(1)
        if not os.path.isfile(ref_index_suffix1):
            keep_logging('The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), 'The {} reference index files were not created properly. Please try to create the index files again or manually.'.format(aligner), logger, 'exception')

    else:
        print "Different Aligner in config file"

def generate_index(reference):
    if not os.path.isfile(reference):
        file_basename = os.path.basename(reference)
        keep_logging(
            'The reference fasta file {} does not exists. Please provide another with full path file with full path or check the files path.\n'.format(
                file_basename),
            'The reference fasta file {} does not exists. Please provide another file or check the files path.\n'.format(
                file_basename), logger, 'exception')
        exit()
    if ConfigSectionMap("pipeline", Config)['aligner'] == "bwa":
        ref_index_suffix1 = reference + ".bwt"
        ref_index_suffix2 = reference + ".amb"
        ref_index_suffix3 = reference + ".ann"
        ref_index_suffix4 = reference + ".sa"
        ref_index_suffix5 = reference + ".pac"
    elif ConfigSectionMap("pipeline", Config)['aligner'] == "bowtie":
        ref_index_suffix1 = reference + ".1.bt2"
        ref_index_suffix2 = reference + ".2.bt2"
        ref_index_suffix3 = reference + ".3.bt2"
        ref_index_suffix4 = reference + ".4.ebwt"
        ref_index_suffix5 = reference + ".rev.1.bt2"
        ref_index_suffix6 = reference + ".rev.2.bt2"
    if not os.path.isfile(ref_index_suffix1):
        keep_logging(
            '',
            'The reference index files given below does not exists:\n {}\n {}\n {}\n {}\n {}'.format(ref_index_suffix1,
                                                                                                     ref_index_suffix2,
                                                                                                     ref_index_suffix3,
                                                                                                     ref_index_suffix4,
                                                                                                     ref_index_suffix5),
            logger, 'warning')
        create_index(reference, ref_index_suffix1, ref_index_suffix2, ref_index_suffix3, ref_index_suffix4,
                     ref_index_suffix5)
    else:
        keep_logging('Index file already exists.', 'Index file already exists.', logger, 'info')

    ref_fai_index = reference + ".fai"
    if not os.path.isfile(ref_fai_index):
        keep_logging('Creating FAI Index using Samtools.', 'The reference fai index file {} required for samtools does not exists.'.format(ref_fai_index), logger,
                     'warning')
        create_fai_index(reference, ref_fai_index)
    else:
        keep_logging('Samtools fai Index file already exists.', 'Samtools fai Index file already exists.', logger, 'info')

    dict_name = os.path.splitext(os.path.basename(reference))[0] + ".dict"
    if not os.path.isfile(ConfigSectionMap(args.index, Config)['ref_path'] + "/" + dict_name):
        keep_logging('Creating Sequence Dictionary using Picard.',
                     'The reference seq dict file {} required for GATK and PICARD does not exists.'.format(dict_name),
                     logger, 'warning')
        picard_seqdict(dict_name, reference)
    else:
        keep_logging('The reference seq dict file required for GATK and PICARD exists.',
                     'The reference seq dict file required for GATK and PICARD exists.', logger, 'info')

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

    if args.output != '':
        if not (args.output).endswith('/'):
            args.output += '/'
    make_sure_path_exists(args.output)

    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')

    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"


    logs_folder = args.output + "Logs"
    make_sure_path_exists(logs_folder)
    # logger = generate_logger(logs_folder, args.prefix, log_unique_time)

    files_to_delete = []
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)

    # Run pipeline steps
    if "core_All" not in args.steps and "core" not in args.steps and "core_prep" not in args.steps and "report" not in args.steps and "tree" not in args.steps and "2" not in args.steps:
        """ Set Up variant calling logs folder/logger object, cluster mode and copy config files to it"""
        vc_logs_folder = logs_folder + "/variant_calling"
        make_sure_path_exists(vc_logs_folder)
        logger = generate_logger(vc_logs_folder, args.prefix, log_unique_time)
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, vc_logs_folder, log_unique_time, args.prefix), logger)
        if args.cluster:
            cluster_mode = args.cluster
        else:
            cluster_mode = "local"

        """ Generate Reference Genome Index """
        # Reference Genome file name
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
        

        keep_logging('Getting Reference Genome name from config file.',
                     'Getting Reference Genome name from config file.', logger, 'info')
        keep_logging('Reference Genome: {}'.format(reference),
                     'Reference Genome: {}'.format(reference), logger, 'info')

        generate_index(reference)
        
        prepare_snpEff_db(reference, vc_logs_folder, logger, Config)

        """ Main Variant calling Methods: Generate and Run the jobs"""
        list_of_files = get_filenames(args.dir, args.type, args.filenames, args.prefix, args.suffix)
        list_of_jobs = create_varcall_jobs(list_of_files, args.type, args.output, args.index, args.steps, config_file, logger)
        job_submitted = run_varcall_jobs(list_of_jobs, cluster_mode, log_unique_time, args.prefix, args.output, logger)

        time_taken = datetime.now() - start_time_2
        keep_logging('Logs were recorded in file with extension log.txt in %s' % vc_logs_folder, 'Logs were recorded in file with extension log.txt in %s' % vc_logs_folder, logger, 'info')
        keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')
        keep_logging('End: Variant calling Pipeline', 'End: Variant calling Pipeline', logger, 'info')

    elif "core_All" in args.steps or "2" in args.steps:
        # Add Core_All step Commands to this array
        core_All_cmds = []

        """ Set Up Core Prep logs folder/logger object, cluster mode and copy config files to it"""
        core_prep_logs_folder = logs_folder + "/core_prep"
        make_sure_path_exists(core_prep_logs_folder)
        core_temp_dir = args.output + "core_temp_dir/"
        
        logger = generate_logger(core_prep_logs_folder, args.prefix, log_unique_time)
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, core_prep_logs_folder, log_unique_time, args.prefix), logger)
        make_sure_path_exists(core_temp_dir)
        
        
        # Remove Copying data files step. THis is redundant.
        # keep_logging('\nCopying vcf files to %s\n' % core_temp_dir, '\nCopying vcf files to %s\n' % core_temp_dir, logger, 'info')
        # cp_command = "cp %s/*/*_vcf_results/*_filter2_indel_final.vcf %s/*/*_vcf_results/*_aln_mpileup_raw.vcf %s/*/*_vcf_results/*_raw.vcf_5bp_indel_removed.vcf* %s/*/*_vcf_results/*filter2_final.vcf* %s/*/*_vcf_results/*vcf_no_proximate_snp.vcf* %s/*/*_vcf_results/*array %s/*/*unmapped.bed_positions %s/*/*_vcf_results/*_indel_gatk.vcf %s/*/*_stats_results/*_depth_of_coverage.sample_summary %s" % (args.output, args.output, args.output, args.output, args.output, args.output, args.output, args.output, args.output, core_temp_dir)

        # cp_command_1 = "cp %s/*/*_vcf_results/*_filter2_indel_final.vcf* %s" % (args.output, core_temp_dir)
        # cp_command_2 = "cp %s/*/*_vcf_results/*_aln_mpileup_raw.vcf* %s" % (args.output, core_temp_dir)
        # cp_command_3 = "cp %s/*/*_vcf_results/*_raw.vcf_5bp_indel_removed.vcf* %s" % (args.output, core_temp_dir)
        # cp_command_4 = "cp %s/*/*_vcf_results/*filter2_final.vcf* %s" % (args.output, core_temp_dir)
        # cp_command_5 = "cp %s/*/*_vcf_results/*array %s" % (args.output, core_temp_dir)
        # cp_command_6 = "cp %s/*/*_vcf_results/*vcf_no_proximate_snp.vcf* %s" % (args.output, core_temp_dir)
        # cp_command_7 = "cp %s/*/*unmapped.bed_positions %s" % (args.output, core_temp_dir)
        # cp_command_8 = "cp %s/*/*_vcf_results/*_indel_gatk.vcf* %s" % (args.output, core_temp_dir)
        # cp_command_9 = "cp %s/*/*_stats_results/*_depth_of_coverage.sample_summary %s" % (args.output, core_temp_dir)

        # #call(cp_command, logger)
        # call(cp_command_1, logger)
        # call(cp_command_2, logger)
        # call(cp_command_3, logger)
        # call(cp_command_4, logger)
        # call(cp_command_5, logger)
        # call(cp_command_6, logger)
        # call(cp_command_7, logger)
        # call(cp_command_8, logger)
        # call(cp_command_9, logger)

        # No need to decompress the files. THis is redundant.
        """ Decompress zipped files in core temp folder"""
        # list_of_gzipped_files = glob.glob("%s/*.gz" % core_temp_dir)
        # keep_logging('Decompressing gzipped files in %s\n' % core_temp_dir, 'Decompressing gzipped files in %s\n' % core_temp_dir, logger, 'info')
        # gzipped_command_list = []
        # for i in list_of_gzipped_files:
        #     gzipped_command_list.append("gzip -df %s" % i)

        """ Set the number of cores variable for Parallel-local mode """
        #THis is redundant.
        # num_cores = multiprocessing.cpu_count()
        # results = Parallel(n_jobs=num_cores)(delayed(run_command_list)(i) for i in gzipped_command_list)
        
        # Get the filter criteria name from config file.
        filter_criteria = ConfigSectionMap("SNP_filters", Config)['filter_criteria']

        """ Generate a custom vcf file list to process for core prep step. If file is not provided, it will consider all the samples in output folder"""
        """ Perform Sanity check to confirm that variant calling results are present in core temp folder """
        if args.filenames:
            list_of_files = get_filenames(args.dir, args.type, args.filenames, args.prefix, args.suffix)
            list_of_vcf_files = generate_custom_vcf_file_list(sorted(list_of_files), logger)
            keep_logging('\nNumber of final variant call vcf files: %s' % len(list_of_vcf_files), '\nNumber of final variant call vcf files: %s' % len(list_of_vcf_files), logger, 'info')
            empty_files = []
            with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
                for file in list_of_files:
                    with open(file.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_depth_of_coverage.sample_summary')) as fp:
                        for line in fp:
                            line = line.strip()
                            if line.startswith('Total'):
                                linesplit = line.split('\t')
                                cov_depth = int(float(linesplit[2].strip()))
                                print cov_depth

                                if float(cov_depth) < float(ConfigSectionMap(filter_criteria, Config)['dp']):
                                    keep_logging('The coverage depth for Sample %s - %s is lower than the threshold' % (
                                    os.path.basename(file), float(cov_depth.strip())),
                                                 'The coverage depth for Sample %s - %s is lower than the threshold' % (
                                                 os.path.basename(file), float(cov_depth.strip())), logger, 'info')
                    fp.close()
                    out_fp.write(os.path.basename(file.replace('R1_001.fastq.gz', '_filter2_final.vcf_no_proximate_snp.vcf')) + '\n')

                    # Check if the vcf files are empty
                    if os.stat(file).st_size == 0:
                        empty_files.append(file)
            out_fp.close()

            # Check if the vcf files are empty
            if len(empty_files) > 0:
                keep_logging(
                    'These vcf files doesnt contain any data - \n%s. please rerun variant call jobs for these samples' % empty_files,
                    'These vcf files doesnt contain any data - \n%s. please rerun variant call jobs for these samples' % empty_files,
                    logger, 'exception')
                exit()
        else:
            keep_logging('- Checking if all the variant calling results exists in %s' % args.output, '- Checking if all the variant calling results exists in %s' % args.output, logger, 'info')
            #call("ls -1a %s/*.vcf_no_proximate_snp.vcf > %s/vcf_filenames" % (core_temp_dir, core_temp_dir), logger)
            try:
                # Not checking the files in core_temp_dir. This is redundant.
                # list_cmd = "ls -1a %s/*.vcf_no_proximate_snp.vcf" % core_temp_dir
                # list_of_files = subprocess.check_output(list_cmd, shell=True)
                
                list_of_files = glob.glob("%s/*/*_vcf_results/*.vcf_no_proximate_snp.vcf" % args.output)
                list_of_GATK_depth_files = glob.glob("%s/*/*_stats_results/*_depth_of_coverage.sample_summary" % args.output)

                keep_logging('- Number of final variant call vcf files: %s' % len(list_of_files), '- Number of final variant call vcf files: %s' % len(list_of_files), logger, 'info')

                empty_files = []
                with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
                    for file in list_of_GATK_depth_files:
                        depth = "grep -vE '^sample|Total' %s | awk -F'\t' '{print $3}'"
                        proc = subprocess.Popen(["grep -vE '^sample|Total' %s | awk -F'\t' '{print $3}'" % file.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_depth_of_coverage.sample_summary')], stdout=subprocess.PIPE, shell=True)
                        (out2, err2) = proc.communicate()
                        #print file
                        cov_depth = int(float(out2.strip()))
                        out_fp.write(os.path.basename(file.replace('_depth_of_coverage.sample_summary', '_filter2_final.vcf_no_proximate_snp.vcf ')) + '\n')
                        if float(out2.strip()) < float(ConfigSectionMap(filter_criteria, Config)['dp']):
                            keep_logging('The coverage depth for Sample %s - %s is lower than the threshold' % (
                            os.path.basename(file), float(out2.strip())),
                                         'The coverage depth for Sample %s - %s is lower than the threshold' % (
                                         os.path.basename(file), float(out2.strip())), logger, 'info')
                        # Check if the vcf files are empty
                        if os.stat(file).st_size == 0:
                            empty_files.append(file)
                        # if float(out2.strip()) > float(ConfigSectionMap(filter_criteria, Config)['dp']):
                        #     out_fp.write(os.path.basename(file)+'\n')
                out_fp.close()

                # Check if the vcf files are empty
                if len(empty_files) > 0:
                    keep_logging(
                        'These vcf files doesnt contain any data - \n%s. please rerun variant call jobs for these samples'  % empty_files,
                        'These vcf files doesnt contain any data - \n%s. please rerun variant call jobs for these samples'  % empty_files,
                        logger, 'exception')

                    exit()
            except:
                keep_logging('Error: The variant calling results were not found in %s. Please check if variant calling step finished properly without any errors. '
                             'This can be done by checking if all the variant call results folder contains final variant call vcf file: *.vcf_no_proximate_snp.vcf file' % core_temp_dir, 'Error: The variant calling results were not found in %s. Please check if variant calling step finished properly without any errors. '
                                                                                                                                                                                              'This can be done by checking if all the variant call results folder contains final variant call vcf file: *.vcf_no_proximate_snp.vcf file' % core_temp_dir, logger, 'exception')
                exit()


        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
        core_prep_pipeline_cmd = run_core_prep_analysis(core_temp_dir, reference, args.prefix, log_unique_time, args.cluster, logger, config_file)
        core_All_cmds.append(core_prep_pipeline_cmd)
        time_taken = datetime.now() - start_time_2
        keep_logging('- Core Prep Logs will be recorded in %s' % core_prep_logs_folder, '- Core Prep Logs will be recorded in %s' % core_prep_logs_folder, logger, 'info')

        """ Generate Core Variants from core_prep intermediate files """
        core_logs_folder = logs_folder + "/core"
        make_sure_path_exists(core_logs_folder)
        logger = generate_logger(core_logs_folder, args.prefix, log_unique_time)
        #keep_logging('START: Extract core snps and generate diagnostic plots','START: Extract core snps and generate diagnostic plots', logger, 'info')
        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, core_logs_folder, log_unique_time, args.prefix),
             logger)
        core_temp_dir = args.output + "/core_temp_dir/"
        
        core_results_dir = args.output + "/%s_core_results/" % log_unique_time
        #make_sure_path_exists(core_results_dir)

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
            keep_logging('- Parsing Functional class filters', '- Parsing Functional class filters', logger,
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
                            '- Mask positions in this file %s will be filtered out' % mask_positions_file,
                            '- Mask positions in this file %s will be filtered out' % mask_positions_file,
                            logger, 'info')
                    else:
                        # mask_positions_file = mask_file
                        os.system("cp %s %s/mask_positions.txt" % (mask_file, core_temp_dir))
                        mask_positions_file = "%s/mask_positions.txt" % core_temp_dir
                        keep_logging(
                            '- Mask positions in this file %s will be filtered out' % mask_positions_file,
                            '- Mask positions in this file %s will be filtered out' % mask_positions_file,
                            logger, 'info')
                    with open(mask_positions_file, 'rU') as fp:
                        for line in fp:
                            f1.write(line)
                    fp.close()
            f1.close()

        copy_phage_results = "cp %s/summary.txt %s/phage_region_positions.txt %s/detail.txt" % (os.path.dirname(reference), os.path.dirname(reference), os.path.dirname(reference))
        call(copy_phage_results, logger)

        core_pipeline_cmd = run_core_analysis(core_temp_dir, reference, args.prefix, log_unique_time, args.cluster, logger,
                          core_results_dir, config_file)
        core_All_cmds.append(core_pipeline_cmd)
        time_taken = datetime.now() - start_time_2
        keep_logging('- Core Step Logs will be recorded in file with extension log.txt in %s' % core_logs_folder,
                     '- Core Step Logs will be recorded in file with extension log.txt in %s' % core_logs_folder, logger, 'info')


        """ Generate Reports and organize core results folder """
        report_logs_folder = logs_folder + "/report"
        make_sure_path_exists(report_logs_folder)
        logger = generate_logger(report_logs_folder, args.prefix, log_unique_time)

        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, report_logs_folder, log_unique_time, args.prefix),
             logger)
        core_temp_dir = args.output + "/core_temp_dir/"

        proc = subprocess.Popen(["ls -1ad %s/*_core_results | tail -n1" % args.output], stdout=subprocess.PIPE,
                                shell=True)
        (out2, err2) = proc.communicate()
        core_results_dir = out2.strip()


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
        run_report_analysis_cmd = run_report_analysis(core_temp_dir, reference, args.prefix, log_unique_time, args.cluster, logger,
                            core_results_dir, config_file)
        time_taken = datetime.now() - start_time_2
        keep_logging('- Report Step Logs will be recorded in file with extension log.txt in %s' % report_logs_folder,
                     '- Report Step Logs will be recorded in file with extension log.txt in %s' % report_logs_folder, logger, 'info')


        core_All_cmds.append(run_report_analysis_cmd)

        """ Generate Trees """
        tree_logs_folder = logs_folder + "/tree"
        make_sure_path_exists(tree_logs_folder)
        logger = generate_logger(tree_logs_folder, args.prefix, log_unique_time)

        call("cp %s %s/%s_%s_config_copy.txt" % (config_file, tree_logs_folder, log_unique_time, args.prefix),
             logger)
        core_temp_dir = args.output + "/core_temp_dir/"

        # This code snippet is not working, Commenting it for now. Test it later.
        # proc = subprocess.Popen(["ls -1ad %s/*_core_results | tail -n1" % args.output], stdout=subprocess.PIPE,
        #                         shell=True)
        # (out2, err2) = proc.communicate()
        # core_results_dir = out2.strip()

        #print core_results_dir
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
        run_tree_analysis_cmd = run_tree_analysis(core_temp_dir, reference, args.prefix, log_unique_time,
                                                  args.cluster, logger, core_results_dir, config_file)
        time_taken = datetime.now() - start_time_2
        keep_logging('- Tree Step Logs will be recorded in file with extension log.txt in %s' % tree_logs_folder,
                     '- Tree Step Logs will be recorded in file with extension log.txt in %s' % tree_logs_folder, logger, 'info')

        core_All_cmds.append(run_tree_analysis_cmd)

        if args.scheduler == "SLURM":
            combine_job_name = core_temp_dir + "/" + log_unique_time + "_" + args.prefix + "_core_All.sbat"
            keep_logging("- Submit Job: %s" % combine_job_name, "- Submit Job: %s" % combine_job_name, logger, 'info')
            scheduler_directives, script_Directive, job_name_flag = get_scheduler_directive(args.scheduler, Config)
            with open(combine_job_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(combine_job_name))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % core_temp_dir + '\n')
                for cmds in core_All_cmds:
                    out.write(cmds + '\n')
            out.close()
            if not args.dryrun:
                proc = subprocess.Popen(["sbatch %s" % combine_job_name], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()
        else:
            combine_job_name = core_temp_dir + "/" + log_unique_time + "_" + args.prefix + "_core_All.pbs"
            keep_logging("qsub %s" % combine_job_name, "qsub %s" % combine_job_name, logger, 'info')
            scheduler_directives, script_Directive, job_name_flag = get_scheduler_directive(args.scheduler, Config)
            with open(combine_job_name, 'w') as out:
                job_title = "%s %s%s" % (script_Directive, job_name_flag, os.path.basename(combine_job_name))
                out.write("#!/bin/sh" + '\n')
                out.write(job_title + '\n')
                out.write(scheduler_directives + '\n')
                out.write("cd %s/" % core_temp_dir + '\n')
                for cmds in core_All_cmds:
                    out.write(cmds + '\n')
            out.close()
            if not args.dryrun:
                proc = subprocess.Popen(["qsub %s" % combine_job_name], stdout=subprocess.PIPE, shell=True)
                (out, err) = proc.communicate()

    else:
        logger = generate_logger(logs_folder, args.prefix, log_unique_time)
        keep_logging('Please provide argument -steps to run pipeline', 'Please provide argument -steps to run pipeline', logger, 'info')
        time_taken = datetime.now() - start_time_2
        keep_logging('Logs were recorded in file with extension log.txt in %s' % logs_folder, 'Logs were recorded in file with extension log.txt in %s' % logs_folder, logger, 'info')
        keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')
