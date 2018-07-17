# import ConfigParser
# from config_settings import ConfigSectionMap
# from logging_subprocess import *
# from log_modules import *
# import os
# import csv
# import glob
# from joblib import Parallel, delayed
# import multiprocessing
#
#
# """ core_prep methods """
# def run_command(i):
#     """
#     Function to run each command and is run as a part of python Parallel mutiprocessing method.
#     :param: command
#     :return: Done status
#     """
#     call("%s" % i, logger)
#     done = "Completed: %s" % i
#     return done
#
# def create_positions_filestep(vcf_filenames, filter2_only_snp_vcf_dir, temp_dir, logger, Config):
#
#     """
#     This method gathers SNP positions from each final *_no_proximate_snp.vcf file (these are the positions that passed variant filter parameters
#     from variant calling pipeline) and write to *_no_proximate_snp.vcf_position files. Use these *_no_proximate_snp.vcf_position files to generate a list of unique_position_file
#     :param: list of final vcf filenames i.e *.vcf_no_proximate_snp.vcf . These files are the final output of variant calling step for each sample.
#     :return: unique_position_file
#     """
#
#     filter2_only_snp_position_files_array = []
#     for file in vcf_filenames:
#         with open(file, 'rU') as csv_file:
#             file_name = temp_dir + "/" + os.path.basename(file) + "_positions"
#             addpositionfilenametoarray = file_name
#             filter2_only_snp_position_files_array.append(addpositionfilenametoarray)
#             f1 = open(file_name, 'w+')
#             csv_reader = csv.reader(csv_file, delimiter='\t')
#             for row in csv_reader:
#                 position = row[0]
#                 if not position.startswith('#'):
#                     p_string = row[1] + "\n"
#                     f1.write(p_string)
#             f1.close()
#         csv_file.close()
#
#     """ Create position array containing unique positiones from positions file """
#
#     position_array = []
#     for filess in filter2_only_snp_position_files_array:
#         f = open(filess, 'r+')
#         for line in f:
#             line = line.strip()
#             position_array.append(line)
#         f.close()
#     position_array_unique = set(position_array)
#     position_array_sort = sorted(position_array_unique)
#     keep_logging('\nThe number of unique variant positions:%s' % len(position_array_sort), '\nThe number of unique variant positions:%s' % len(position_array_sort), logger, 'info')
#     unique_position_file = "%s/unique_positions_file" % filter2_only_snp_vcf_dir
#     f=open(unique_position_file, 'w+')
#     for i in position_array_sort:
#         f.write(i + "\n")
#     f.close()
#     if len(position_array_sort) == 0:
#         keep_logging('ERROR: No unique positions found. Check if vcf files are empty?', 'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
#         exit()
#     return unique_position_file
#
# def create_indel_positions_filestep(vcf_filenames, filter2_only_snp_vcf_dir, temp_dir, logger, Config):
#
#     """
#     This function gathers Indel positions from each final *_indel_final.vcf (these are the positions that passed variant filter parameters
#     from variant calling pipeline) and write to *_indel_final.vcf files. Use these *_indel_final.vcf_position files to generate a list of unique_position_file
#     :param: list of final vcf filenames i.e *_indel_final.vcf . These files are the final output of variant calling step for each sample.
#     :return: unique_indel_position_file
#     """
#
#     filter2_only_indel_position_files_array = []
#     for file in vcf_filenames:
#         indel_file = file.replace('_filter2_final.vcf_no_proximate_snp.vcf', '_filter2_indel_final.vcf')
#         with open(indel_file, 'rU') as csv_file:
#             file_name = temp_dir + "/" + os.path.basename(indel_file) + "_positions"
#             addpositionfilenametoarray = file_name
#             filter2_only_indel_position_files_array.append(addpositionfilenametoarray)
#             f1 = open(file_name, 'w+')
#             csv_reader = csv.reader(csv_file, delimiter='\t')
#             for row in csv_reader:
#                 position = row[0]
#                 if not position.startswith('#'):
#                     p_string = row[1] + "\n"
#                     f1.write(p_string)
#             f1.close()
#         csv_file.close()
#
#     """ Create position array containing unique positiones from positions file """
#     position_array = []
#     for filess in filter2_only_indel_position_files_array:
#         f = open(filess, 'r+')
#         for line in f:
#             line = line.strip()
#             position_array.append(line)
#         f.close()
#     position_array_unique = set(position_array)
#     position_array_sort = sorted(position_array_unique)
#     keep_logging('\nThe number of unique indel positions:%s' % len(position_array_sort), '\nThe number of unique indel positions:%s' % len(position_array_sort), logger, 'info')
#     unique_indel_position_file = "%s/unique_indel_positions_file" % filter2_only_snp_vcf_dir
#     f=open(unique_indel_position_file, 'w+')
#     for i in position_array_sort:
#         f.write(i + "\n")
#     f.close()
#     if len(position_array_sort) == 0:
#         keep_logging('ERROR: No unique positions found. Check if vcf files are empty?', 'ERROR: No unique positions found. Check if vcf files are empty?', logger, 'info')
#         exit()
#     return unique_indel_position_file
#
# def create_job(jobrun, vcf_filenames, unique_position_file, tmp_dir, filter2_only_snp_vcf_dir, logger, Config):
#
#     """
#     This method takes the unique_position_file and list of final *_no_proximate_snp.vcf files and generates individual jobs/script.
#     Each of these jobs/scripts will generate a *label file. These label file for each sample contains a field description for each position in unique_position_file.
#     This field description denotes if the variant position made to the final variant list in a sample and if not then a reason/filter that caused it to filtered out from final list.
#     :param jobrun:
#     :param vcf_filenames:
#     :return:
#     """
#     if jobrun == "parallel-cluster":
#         """
#         Supports only PBS clusters for now.
#         """
#         for i in vcf_filenames:
#             job_name = os.path.basename(i)
#             job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
#             job_file_name = "%s.pbs" % (i)
#             f1=open(job_file_name, 'w+')
#             f1.write(job_print_string)
#             f1.close()
#         #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
#         pbs_dir = filter2_only_snp_vcf_dir + "/*vcf.pbs"
#         pbs_scripts = glob.glob(pbs_dir)
#         for i in pbs_scripts:
#             keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
#             call("qsub %s" % i, logger)
#
#     elif jobrun == "parallel-local":
#         """
#         Generate a Command list of each job and run it in parallel on different cores available on local system
#         """
#         command_array = []
#         command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
#         f3 = open(command_file, 'w+')
#
#
#         for i in vcf_filenames:
#             job_name = os.path.basename(i)
#             job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
#             job_file_name = "%s.pbs" % (i)
#             f1=open(job_file_name, 'w+')
#             f1.write(job_print_string)
#             f1.close()
#         #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
#         pbs_dir = filter2_only_snp_vcf_dir + "/*vcf.pbs"
#         pbs_scripts = glob.glob(pbs_dir)
#         for i in pbs_scripts:
#             f3.write("bash %s\n" % i)
#         f3.close()
#         with open(command_file, 'r') as fpp:
#             for lines in fpp:
#                 lines = lines.strip()
#                 command_array.append(lines)
#         fpp.close()
#         # if args.numcores:
#         #     num_cores = int(num_cores)
#         # else:
#         num_cores = multiprocessing.cpu_count()
#         print num_cores
#         print command_array
#         results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)
#
#     elif jobrun == "cluster":
#         #command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
#         #os.system("bash %s" % command_file)
#         command_array = []
#         command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
#         f3 = open(command_file, 'w+')
#
#
#         for i in vcf_filenames:
#             job_name = os.path.basename(i)
#             job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
#             job_file_name = "%s.pbs" % (i)
#             f1=open(job_file_name, 'w+')
#             f1.write(job_print_string)
#             f1.close()
#         #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
#         pbs_dir = filter2_only_snp_vcf_dir + "/*vcf.pbs"
#         pbs_scripts = glob.glob(pbs_dir)
#         for i in pbs_scripts:
#             f3.write("bash %s\n" % i)
#         f3.close()
#         with open(command_file, 'r') as fpp:
#             for lines in fpp:
#                 lines = lines.strip()
#                 command_array.append(lines)
#         fpp.close()
#         # if args.numcores:
#         #     num_cores = int(num_cores)
#         # else:
#         num_cores = multiprocessing.cpu_count()
#         results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)
#
#     elif jobrun == "local":
#         """
#         Generate a Command list of each job and run it on local system one at a time
#         """
#
#         command_array = []
#         command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
#         f3 = open(command_file, 'w+')
#
#
#         for i in vcf_filenames:
#             job_name = os.path.basename(i)
#             job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
#             job_file_name = "%s.pbs" % (i)
#             f1=open(job_file_name, 'w+')
#             f1.write(job_print_string)
#             f1.close()
#         #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
#         pbs_dir = filter2_only_snp_vcf_dir + "/*vcf.pbs"
#         pbs_scripts = glob.glob(pbs_dir)
#
#
#         for i in pbs_scripts:
#             f3.write("bash %s\n" % i)
#         f3.close()
#         with open(command_file, 'r') as fpp:
#             for lines in fpp:
#                 lines = lines.strip()
#                 command_array.append(lines)
#         fpp.close()
#         call("bash %s" % command_file, logger)
#
# def create_indel_job(jobrun, vcf_filenames, unique_position_file, tmp_dir, filter2_only_snp_vcf_dir, logger, Config):
#
#     """
#     This method takes the unique_indel_position_file and list of final *_indel_final.vcf files and generates individual jobs/script.
#     Each of these jobs/scripts will generate a *label file. These label file for each sample contains a field description of each position in unique_indel_position_file.
#     This field description denotes if the variant position made to the final variant list in a sample and if not then a reason/filter that caused it to filtered out from final list.
#     :param jobrun:
#     :param vcf_filenames:
#     :return:
#     """
#     if jobrun == "parallel-cluster":
#         """
#         Supports only PBS clusters for now.
#         """
#         for i in vcf_filenames:
#             job_name = os.path.basename(i)
#             job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_indel_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
#             job_file_name = "%s_indel.pbs" % (i)
#             f1=open(job_file_name, 'w+')
#             f1.write(job_print_string)
#             f1.close()
#         #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
#         pbs_dir = filter2_only_snp_vcf_dir + "/*vcf_indel.pbs"
#         pbs_scripts = glob.glob(pbs_dir)
#         for i in pbs_scripts:
#             keep_logging('Running: qsub %s' % i, 'Running: qsub %s' % i, logger, 'info')
#             # os.system("qsub %s" % i)
#             call("qsub %s" % i, logger)
#
#     elif jobrun == "parallel-local" or jobrun == "cluster":
#         """
#         Generate a Command list of each job and run it in parallel on different cores available on local system
#         """
#         command_array = []
#         command_file = "%s/commands_indel_list.sh" % filter2_only_snp_vcf_dir
#         f3 = open(command_file, 'w+')
#
#
#         for i in vcf_filenames:
#             job_name = os.path.basename(i)
#             job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_indel_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, ConfigSectionMap("scheduler", Config)['email'], ConfigSectionMap("scheduler", Config)['notification'], ConfigSectionMap("scheduler", Config)['resources'], ConfigSectionMap("scheduler", Config)['queue'], ConfigSectionMap("scheduler", Config)['flux_account'], args.filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
#             job_file_name = "%s_indel.pbs" % (i)
#             f1=open(job_file_name, 'w+')
#             f1.write(job_print_string)
#             f1.close()
#         #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
#         pbs_dir = filter2_only_snp_vcf_dir + "/*vcf_indel.pbs"
#         pbs_scripts = glob.glob(pbs_dir)
#         for i in pbs_scripts:
#             f3.write("bash %s\n" % i)
#         f3.close()
#         with open(command_file, 'r') as fpp:
#             for lines in fpp:
#                 lines = lines.strip()
#                 command_array.append(lines)
#         fpp.close()
#         # if args.numcores:
#         #     num_cores = int(num_cores)
#         # else:
#         num_cores = multiprocessing.cpu_count()
#         results = Parallel(n_jobs=num_cores)(delayed(run_command)(command) for command in command_array)
#
#     # elif jobrun == "cluster":
#     #     command_file = "%s/commands_list.sh" % args.filter2_only_snp_vcf_dir
#     #     os.system("bash %s" % command_file)
#     elif jobrun == "local":
#         """
#         Generate a Command list of each job and run it on local system one at a time
#         """
#
#         command_array = []
#         command_file = "%s/commands_list.sh" % filter2_only_snp_vcf_dir
#         f3 = open(command_file, 'w+')
#
#
#         for i in vcf_filenames:
#             job_name = os.path.basename(i)
#             job_print_string = "#PBS -N %s\n#PBS -M %s\n#PBS -m %s\n#PBS -V\n#PBS -l %s\n#PBS -q %s\n#PBS -A %s\n#PBS -l qos=flux\n\n/nfs/esnitkin/bin_group/anaconda2/bin/python /nfs/esnitkin/bin_group/pipeline/Github/variant_calling_pipeline_dev/modules/variant_diagnostics/reason_job_indel_debug.py -filter2_only_snp_vcf_dir %s -filter2_only_snp_vcf_file %s -unique_position_file %s -tmp_dir %s\n" % (job_name, filter2_only_snp_vcf_dir, i, unique_position_file, tmp_dir)
#             job_file_name = "%s_indel.pbs" % (i)
#             f1=open(job_file_name, 'w+')
#             f1.write(job_print_string)
#             f1.close()
#         #os.system("mv %s/*.pbs %s/temp" % (args.filter2_only_snp_vcf_dir, args.filter2_only_snp_vcf_dir))
#         pbs_dir = filter2_only_snp_vcf_dir + "/*vcf_indel.pbs"
#         pbs_scripts = glob.glob(pbs_dir)
#
#
#         for i in pbs_scripts:
#             f3.write("bash %s\n" % i)
#         f3.close()
#         with open(command_file, 'r') as fpp:
#             for lines in fpp:
#                 lines = lines.strip()
#                 command_array.append(lines)
#         fpp.close()
#         call("bash %s" % command_file, logger)
