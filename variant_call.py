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


# Command Line Argument Parsing
def parser():
    parser = argparse.ArgumentParser(description='\nVariant Calling pipeline for Illumina PE/SE data.\n', formatter_class=RawTextHelpFormatter)
    required = parser.add_argument_group('Required arguments')
    optional = parser.add_argument_group('Optional arguments')
    required.add_argument('-type', action='store', dest="type", help='Type of reads: SE or PE', required=True)
    required.add_argument('-readsdir', action='store', dest="dir", help='Path to Sequencing Reads Data directory. NOTE: Provide full/absolute path.', required=True)
    required.add_argument('-outdir', action='store', dest="output_folder", help='Output Folder Path ending with output directory name to save the results. Creates a new output directory path if it doesn\'t exist. NOTE: Provide full/absolute path.', required=True)
    required.add_argument('-index', action='store', dest="index", help='Reference Index Name. Most Frequently used reference genomes index options: KPNIH1 | MRSA_USA_300 | MRSA_USA_100 | CDIFF_630 | paris Make sure the paths are properly set in config file', required=True)
    required.add_argument('-steps', action='store', dest="steps", help='Variant Calling Steps in sequential order.\n'
                                                                     '1.   All: This will run all the steps starting from cleaning the reads to variant calling;\n'
                                                                     '2.   clean,align,post-align,varcall,filter,stats : This will also run all steps starting from cleaning to variant calling. \nYou can also run part of the pipeline by giving "align,post-align,varcall,filter,stats" which will skip the cleaning part.\nThe order is required to be sequential while using this option. Also, while skipping any of the step make sure you have results already present in your output folder.\n'
                                                                     '3.   coverage_depth_stats: Run Only Depth of Coverage Stats module after cleaning and read mapping steps\n'
                                                                     '4.   core_prep: Run this step before running the core steps. This will prepare the data required for generating core SNPs\n'
                                                                     '5.   core: extract core snps and generate diagnostics plot data matrices to explore filtered snps.')
    required.add_argument('-analysis', action='store', dest="analysis_name", help='Unique analysis name that will be used as prefix to saving results and log files.', required=True)
    optional.add_argument('-config', action='store', dest="config", help='Path to Config file, Make sure to check config settings before running pipeline', required=False)
    optional.add_argument('-suffix', action='store', dest="suffix", help='Fastq reads suffix such as fastq, fastq.gz, fq.gz, fq; Default: fastq.gz', required=False)
    optional.add_argument('-filenames', action='store', dest="filenames", help='fastq filenames with one single-end filename per line. \nIf the type is set to PE, it will detect the second paired-end filename with the suffix from first filename. \nUseful for running variant calling pipeline on selected files in a reads directory or extracting core snps for selected samples in input reads directory. \nOtherwise the pipeline will consider all the samples available in reads directory.', required=False)
    optional.add_argument('-cluster', action='store', dest='cluster', help='Run variant calling pipeline in one of the four modes. Default: local. Suggested mode for core snp is cluster that will run all the steps in parallel with the available cores. Make sure to provide a large memory node for this option\nThe possible modes are: cluster/parallel-cluster/parallel-local/local\ncluster: Runs all the jobs on a single large cluster. This will mimic the local run but rather on a large compute node.\nparallel-cluster: Submit variant call jobs for each sample in parallel on compute nodes. This mode is no available for core snp extraction step.\nparallel-local: Run variant call jobs for each sample in parallel locally.\nlocal: Run variant call jobs locally.\nMake Sure to check if the [scheduler] section in config file is set up correctly for your cluster.')
    return parser

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
                print "No fastq files with suffix %s found in reads directory %s" % (suffix, dir)
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

def create_varcall_jobs(filenames_array, type, output_folder, reference, steps, config_file, logger):
    jobs_temp_dir = "%s/temp_jobs" % output_folder
    make_sure_path_exists(jobs_temp_dir)
    keep_logging('Generating cluster jobs in temporary directory %s' % jobs_temp_dir, 'Generating cluster jobs in temporary directory %s' % jobs_temp_dir, logger, 'exception')
    ##Change the email address; nodes and processor requirements accordingly

    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])

    for file in filenames_array:
        filename_base = os.path.basename(file)
        if "R1_001_final.fastq.gz" in filename_base or "R1.fastq.gz" in filename_base or "1_combine.fastq.gz" in filename_base or "1_sequence.fastq.gz" in filename_base or "_forward.fastq.gz" in filename_base or "R1_001.fastq.gz" in filename_base or "_1.fastq.gz" in filename_base or ".1.fastq.gz" in filename_base or "_R1.fastq.gz" in filename_base:
            # Forward reads file name and get analysis name from its name
            first_file = file
            # Get the name of reverse reads files
            if "R1_001_final.fastq.gz" in filename_base:
                second_part = filename_base.replace("R1_001_final.fastq.gz", "R2_001_final.fastq.gz")
                first_part_split = filename_base.split('R1_001_final.fastq.gz')
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
            elif "_R1.fastq.gz" in filename_base:
                second_part = filename_base.replace("_R1.fastq.gz", "_R2.fastq.gz")
                first_part_split = filename_base.split('_R1.fastq.gz')
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
                command = "/home/apirani/anaconda/bin/python %s/pipeline.py -PE1 %s -o %s/%s -analysis %s -index %s -type SE -config %s -steps %s" % (os.path.dirname(os.path.abspath(__file__)), first_file, output_folder, first_part, first_part, reference, config_file, steps)
            else:
                command = "/home/apirani/anaconda/bin/python %s/pipeline.py -PE1 %s -PE2 %s -o %s/%s -analysis %s -index %s -type PE -config %s -steps %s" % (os.path.dirname(os.path.abspath(__file__)), first_file, second_file, output_folder, first_part, first_part, reference, config_file, steps)

            with open(job_name, 'w') as out:
                job_title = "#PBS -N %s" % first_part
                out.write(job_title+'\n')
                out.write(Pbs_model_lines+'\n')
                #out.write(cd_command+'\n')
                out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
                out.write("echo $PBS_O_WORKDIR" + '\n')
                out.write("cd %s/temp_jobs" % output_folder + '\n')
                out.write(command+'\n')
        elif "R2_001_final.fastq.gz" in filename_base or "R2.fastq.gz" in filename_base or "2_combine.fastq.gz" in filename_base or "2_sequence.fastq.gz" in filename_base or "_reverse.fastq.gz" in filename_base or "R2_001.fastq.gz" in filename_base or "_2.fastq.gz" in filename_base or ".2.fastq.gz" in filename_base or "_R2.fastq.gz" in filename_base:
            continue
        else:
            keep_logging('Error while generating cluster jobs. Make sure the fastq filenames ends with one of these suffix: R1_001_final.fastq.gz, R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, R1_001.fastq.gz, _1.fastq.gz, .1.fastq.gz, _R1.fastq.gz', 'Error while generating cluster jobs. Make sure the fastq filenames ends with one of these suffix: R1_001_final.fastq.gz, R1.fastq.gz, 1_combine.fastq.gz, 1_sequence.fastq.gz, _forward.fastq.gz, R1_001.fastq.gz, _1.fastq.gz, .1.fastq.gz, _R1.fastq.gz', logger, 'exception')
            exit()
    list_of_jobs = glob.glob("%s/*.pbs" % jobs_temp_dir)
    return list_of_jobs


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
        keep_logging('Running Job in cluster mode', 'Running Job in cluster mode', logger, 'info')
        cluster_job_dir = output_folder + "/temp_jobs/cluster/"
        make_sure_path_exists(cluster_job_dir)
        Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])
        job_name = cluster_job_dir + log_unique_time + "_" + analysis_name + ".pbs"
        with open(job_name, 'w') as out:
            job_title = "#PBS -N %s_%s" % (log_unique_time, analysis_name)
            out.write(job_title+'\n')
            out.write(Pbs_model_lines+'\n')
            out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
            out.write("echo $PBS_O_WORKDIR" + '\n')
            out.write("cd %s/temp_jobs" % output_folder + '\n')
            out.write(command_list+'\n')
        out.close()
        keep_logging('Submitting single cluster Job: qsub %s' % job_name, 'Submitting single cluster Job: qsub %s' % job_name, logger, 'info')
        #call("qsub %s" % job_name, logger)
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
        keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')
    elif cluster == "parallel-cluster":
        keep_logging('Running Jobs in parallel-cluster mode', 'Running Jobs in parallel-cluster mode', logger, 'info')
        for job in command_list_qsub:
            keep_logging('Submitting Job: qsub %s' % job, 'Submitting Job: qsub %s' % job, logger, 'info')
            call("qsub %s" % job, logger)
        keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')
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


def run_core_prep_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger):
    file_exists(reference)
    core_prep_pipeline = "/home/apirani/anaconda/bin/python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 1 -jobrun %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster)
    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"
    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])
    with open(job_name, 'w') as out:
        job_title = "#PBS -N %s_%s_core" % (log_unique_time, analysis_name)
        out.write(job_title+'\n')
        out.write(Pbs_model_lines+'\n')
        out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
        out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
        out.write("cd %s" % core_temp_dir + '\n')
        out.write(core_prep_pipeline+'\n')
    out.close()

    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        keep_logging('Running parallel-local mode: bash %s' % job_name, 'Running parallel-local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        #call("qsub %s" % job_name, logger)
        keep_logging('Submitting single cluster Job: qsub %s' % job_name, 'Submitting single cluster Job: qsub %s' % job_name, logger, 'info')
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
    elif cluster == "parallel-cluster":
        #call("qsub %s" % job_name, logger)
        keep_logging('Submitting parallel-cluster Job: qsub %s' % job_name, 'Submitting parallel-cluster Job: qsub %s' % job_name, logger, 'info')
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
    keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')


def run_core_analysis(core_temp_dir, reference, analysis_name, log_unique_time, cluster, logger):
    file_exists(reference)
    core_pipeline = "/home/apirani/anaconda/bin/python %s/modules/variant_diagnostics/core_pipeline.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_filenames %s/vcf_filenames -reference %s -steps 2 -jobrun %s" % (os.path.dirname(os.path.abspath(__file__)), core_temp_dir, core_temp_dir, reference, cluster)
    job_name = core_temp_dir + "/" + log_unique_time + "_" + analysis_name + ".pbs"
    Pbs_model_lines = "#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n"\
                      % (ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'])
    with open(job_name, 'w') as out:
        job_title = "#PBS -N %s_%s_core" % (log_unique_time, analysis_name)
        out.write(job_title+'\n')
        out.write(Pbs_model_lines+'\n')
        out.write("#  Change to the directory you submitted from\nif [ -n \"$PBS_O_WORKDIR\" ]; then cd $PBS_O_WORKDIR; fi" + '\n')
        out.write("echo \"PBS working directory: $PBS_O_WORKDIR\"" + '\n')
        out.write("cd %s" % core_temp_dir + '\n')
        out.write(core_pipeline+'\n')
    out.close()

    if cluster == "local":
        keep_logging('Running local mode: bash %s' % job_name, 'Running local mode: bash %s' % job_name, logger, 'info')
        call("bash %s" % job_name, logger)
    elif cluster == "parallel-local":
        call("bash %s" % job_name, logger)
    elif cluster == "cluster":
        #call("qsub %s" % job_name, logger)
        keep_logging('Submitting single cluster Job: qsub %s' % job_name, 'Submitting single cluster Job: qsub %s' % job_name, logger, 'info')
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
    elif cluster == "parallel-cluster":
        #call("qsub %s" % job_name, logger)
        keep_logging('Submitting parallel-cluster Job: qsub %s' % job_name, 'Submitting parallel-cluster Job: qsub %s' % job_name, logger, 'info')
        qid = subprocess.check_output("qsub %s" % job_name, shell=True)
        print qid.split('.')[0]
    #keep_logging('You can check the job status with: qstat -u USERNAME', 'You can check the job status with: qstat -u USERNAME', logger, 'info')



# Main Method
if __name__ == '__main__':

    # Set up logging modules and config file
    start_time = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    start_time_2 = datetime.now()
    args = parser().parse_args()
    global config_file
    global log_unique_time
    if args.output_folder != '':
        args.output_folder += '/'
    make_sure_path_exists(args.output_folder)
    log_unique_time = datetime.now().strftime('%Y_%m_%d_%H_%M_%S')
    if args.config:
        config_file = args.config
    else:
        config_file = os.path.dirname(os.path.abspath(__file__)) + "/config"
    global logger
    logger = generate_logger(args.output_folder, args.analysis_name, log_unique_time)
    global Config
    global files_to_delete
    files_to_delete = []
    Config = ConfigParser.ConfigParser()
    Config.read(config_file)
    call("cp %s %s/%s_%s_config_copy.txt" % (config_file, args.output_folder, log_unique_time, args.analysis_name), logger)

    # Run pipeline steps
    if "core" not in args.steps:
        if args.cluster:
            cluster_mode = args.cluster
        else:
            cluster_mode = "local"

        keep_logging('START: Variant Calling Pipeline', 'START: Variant Calling Pipeline', logger, 'info')

        list_of_files = get_filenames(args.dir, args.type, args.filenames, args.analysis_name, args.suffix)
        list_of_jobs = create_varcall_jobs(list_of_files, args.type, args.output_folder, args.index, args.steps, config_file, logger)
        run_varcall_jobs(list_of_jobs, cluster_mode, log_unique_time, args.analysis_name, args.output_folder, logger)
        #pipeline(args, logger)
        keep_logging('End: Variant calling Pipeline', 'End: Variant calling Pipeline', logger, 'info')

    elif "core_prep" in args.steps:
        keep_logging('START: Extract core snps and generate diagnostic plots', 'START: Extract core snps and generate diagnostic plots', logger, 'info')
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        make_sure_path_exists(core_temp_dir)
        keep_logging('Copying vcf files to %s' % core_temp_dir, 'Copying vcf files to %s' % core_temp_dir, logger, 'info')
        cp_command = "cp %s/*/*_aln_mpileup_raw.vcf %s/*/*_raw.vcf_5bp_indel_removed.vcf.gz %s/*/*filter2_final.vcf.gz %s/*/*vcf_no_proximate_snp.vcf.gz %s/*/*array %s/*/*unmapped.bed_positions %s" % (args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, args.output_folder, core_temp_dir)
        call(cp_command, logger)
        list_of_gzipped_files = glob.glob("%s/*.gz" % core_temp_dir)
        keep_logging('Decompressing gzipped files in %s' % core_temp_dir, 'Decompressing gzipped files in %s' % core_temp_dir, logger, 'info')
        gzipped_command_list = []
        for i in list_of_gzipped_files:
            gzipped_command_list.append("gzip -df %s" % i)
        num_cores = multiprocessing.cpu_count()
        results = Parallel(n_jobs=num_cores)(delayed(run_command_list)(i) for i in gzipped_command_list)
        call("ls -1a %s/*.vcf_no_proximate_snp.vcf > %s/vcf_filenames" % (core_temp_dir,core_temp_dir), logger)
        list_cmd = "ls -1a %s/*.vcf_no_proximate_snp.vcf" % core_temp_dir
        list_of_files = subprocess.check_output(list_cmd, shell=True)
        with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
            for file in list_of_files.splitlines():
                out_fp.write(os.path.basename(file)+'\n')
        out_fp.close()
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
        run_core_prep_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger)

    elif "core" in args.steps:
        keep_logging('START: Extract core snps and generate diagnostic plots', 'START: Extract core snps and generate diagnostic plots', logger, 'info')
        core_temp_dir = args.output_folder + "/core_temp_dir/"
        list_of_label_files = glob.glob("%s/*_label" % core_temp_dir)
        list_of_vcf_files = glob.glob("%s/*vcf_no_proximate_snp.vcf" % core_temp_dir)
        if len(list_of_label_files) == len(list_of_vcf_files):
            for i in list_of_label_files:
                if os.stat(i).st_size == 0:
                    keep_logging('The file {} is empty. Please rerun core_prep step again.\n'.format(i), 'The file {} is empty. Please rerun core_prep step again.\n'.format(i), logger, 'exception')
                    exit()
        else:
            keep_logging('Problem in core_prep results. Rerun the core_prep step\n', 'Problem in core_prep results. Rerun the core_prep step\n', logger, 'exception')
            exit()
        list_cmd = "ls -1a %s/*.vcf_no_proximate_snp.vcf" % core_temp_dir
        list_of_files = subprocess.check_output(list_cmd, shell=True)
        with open("%s/vcf_filenames" % core_temp_dir, 'w') as out_fp:
            for file in list_of_files.splitlines():
                out_fp.write(os.path.basename(file)+'\n')
        out_fp.close()
        reference = ConfigSectionMap(args.index, Config)['ref_path'] + "/" + ConfigSectionMap(args.index, Config)['ref_name']
        run_core_analysis(core_temp_dir, reference, args.analysis_name, log_unique_time, args.cluster, logger)

    else:
        keep_logging('Please provide argument -steps to run pipeline', 'Please provide argument -steps to run pipeline', logger, 'info')

    time_taken = datetime.now() - start_time_2
    keep_logging('Logs were recorded in file with extension log.txt in output folder', 'Logs were recorded in file with extension log.txt in output folder', logger, 'info')
    keep_logging('Total Time taken: {}'.format(time_taken), 'Total Time taken: {}'.format(time_taken), logger, 'info')
