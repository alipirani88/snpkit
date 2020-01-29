# Create pbs script based on input

# import libraries
import argparse
import subprocess

# parse command line arguments
parser = argparse.ArgumentParser(
                description='''Create pbs script based on input arguments.''')

parser.add_argument('-c', '--commands', metavar='COMMANDS_FILE', type=str,
                    required=True,
                    help='bash file with commands to run')
parser.add_argument('-o', '--outfile', metavar='OUT_FILE', type=str,
                    required=True,
                    help='output pbs file name')
parser.add_argument('-M', '--modules', metavar='MODULES', type=str,
                    default='',
                    help='modules to load (space-delimited; default: None)')

parser.add_argument('-n', '--nodes', metavar='NUM_NODES', type=str,
                    default='1',
                    help='number of nodes (default: 1)')
parser.add_argument('-p', '--ppn', metavar='NUM_CORES', type=str,
                    default='12',
                    help='number of processers (cores) per node (default: 12)')
parser.add_argument('-m', '--mem', metavar='MEMORY', type=str,
                    default='47',
                    help='amount of memory (gb; default: 47)')
parser.add_argument('-w', '--walltime', metavar='WALL_TIME', type=str,
                    default='10:00:00:00',
                    help='''amount of time needed to run program
                    (default: 10 days -- 10:00:00:00)''')

parser.add_argument('-N', '--jobname', metavar='JOB_NAME', type=str,
                    default=None,
                    help='name of pbs job (default: output file prefix)')
parser.add_argument('-e', '--email', metavar='EMAIL', type=str,
                    default=None,
                    help='''umich email,
                    default: flux login uniqname umich email''')

parser.add_argument('-P', '--pmem', action='store_true',
                    help='''flag - use pmem (memory per processor)
                    instead of mem (total memory)''')

parser.add_argument('-a', '--acct', metavar='ACCOUNT', type=str,
                    default='esnitkin_flux',
                    help='''flux account to submit pbs script to
                    (default: esnitkin_flux)''')
parser.add_argument('-wd', '--wd', metavar='WD', type=str,
                    default='$PBS_O_WORKDIR',
                    help='''directory to submit pbs script from
                    (default: $PBS_O_WORKDIR)''')

args = parser.parse_args()

num_nodes = args.nodes
num_cores = args.ppn
ppn = args.ppn
mem = args.mem
pmem = args.pmem
walltime = args.walltime
acct = args.acct
qos = acct.split('_')[-1]
modules = args.modules
commands = args.commands
outfile = args.outfile
wd = args.wd


if args.jobname is None:
    job_name = outfile.split('.')[0]
else:
    job_name = args.jobname

if args.email is None:
    uniqname = str(subprocess.check_output('whoami').rsplit()).split('\'')[1]
    email = '{}@umich.edu'.format(uniqname)
else:
    email = args.email

if pmem:
    info = '#PBS -l nodes={}:ppn={},pmem={}gb,walltime={}'.format(
           num_nodes, num_cores, mem, walltime)
else:
    info = '#PBS -l nodes={}:ppn={},mem={}gb,walltime={}'.format(
           num_nodes, num_cores, mem, walltime)

# print pbs script to output file
file = open(outfile, 'w')
file.write('#!/bin/sh\n')
file.write('\n')
file.write('#!/bin/sh\n')
file.write('#### PBS preamble\n')
file.write('\n')
file.write('#PBS -N {}\n'.format(job_name))
file.write('\n')
file.write('# User info\n')
file.write('#PBS -M {}\n'.format(email))
file.write('#PBS -m abe\n')
file.write('\n')
file.write('# Number of cores, amount of memory, and walltime\n')
file.write('%s\n' % info)
file.write('#PBS -j oe\n')
file.write('#PBS -V\n')
file.write('\n')
file.write('#PBS -A {}\n'.format(acct))
file.write('#PBS -q {}\n'.format(qos))
file.write('#PBS -l qos=flux\n')
file.write('\n')
file.write('#### End PBS preamble\n')
file.write('\n')
file.write('# Show list of CPUs you ran on, if you\'re running under PBS\n')
file.write('if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi\n')
file.write('\n')
file.write('#  Change to the directory you submitted from\n')
file.write('#cd $PBS_O_WORKDIR\n')
file.write('cd {}\n'.format(wd))
file.write('echo {}\n'.format(wd))
file.write('\n')
file.write('# Load modules\n')
file.write('module load {}\n'.format(modules))
file.write('\n')
file.write('# Job commands\n')
file.write('bash {}\n'.format(commands))
