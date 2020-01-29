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
                    help='''account to submit script to
                    (default: esnitkin_flux)''')
parser.add_argument('-wd', '--wd', metavar='WD', type=str,
                    default=None,
                    help='''directory to submit script from
                    (default: location of script)''')
parser.add_argument('-s', '--slurm', action='store_true',
                    help='''flag - create slurm script (rather than pbs script)''')

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
slurm = args.slurm


if args.wd is None:
    if slurm:
        wd = '$SLURM_SUBMIT_DIR'
    else:
        wd = '$PBS_O_WORKDIR'

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

if slurm:
    
    if pmem:
        info = '#SBATCH --nodes={} --ntasks=1 --cpus-per-task={} --mem-per-cpu={}g --time={}'.format(
        num_nodes, num_cores, mem, walltime)
    else:
        info = '#SBATCH --nodes={} --ntasks=1 --cpus-per-task={} --mem={}g --time={}'.format(
        num_nodes, num_cores, mem, walltime)

    # print slurm script to output file
    print('#!/bin/sh',
          '# Job name',
          '#SBATCH --job-name={}'.format(job_name),
          '# User info',
          '#SBATCH --mail-user={}'.format(email),
          '#SBATCH --mail-type=BEGIN,END,NONE,FAIL,REQUEUE',
          '#SBATCH --export=ALL',
          '#SBATCH --partition=standard',
          '#SBATCH --account={}'.format(acct),
          '# Number of cores, amount of memory, and walltime',
          info,
          '#  Change to the directory you submitted from',
          'cd {}'.format(wd),
          'echo {}'.format(wd),
          '',
          '# Load modules',
          'module load {}'.format(modules),
          '',
          '# Job commands',
          'bash {}'.format(commands),
          sep='\n',
          file=open(outfile, 'w'))
          
else:
    # print pbs script to output file
    print('#!/bin/sh',
          '#### PBS preamble',
          '',
          '#PBS -N {}'.format(job_name),
          '',
          '# User info',
          '#PBS -M {}'.format(email),
          '#PBS -m abe',
          '',
          '# Number of cores, amount of memory, and walltime',
          info,
          '#PBS -j oe',
          '#PBS -V',
          '',
          '#PBS -A {}'.format(acct),
          '#PBS -q {}'.format(qos),
          '#PBS -l qos=flux',
          '',
          '#### End PBS preamble',
          '',
          '# Show list of CPUs you ran on, if you\'re running under PBS',
          'if [ -n "$PBS_NODEFILE" ]; then cat $PBS_NODEFILE; fi',
          '',
          '#  Change to the directory you submitted from',
          '#cd $PBS_O_WORKDIR',
          'cd {}'.format(wd),
          'echo {}'.format(wd),
          '',
          '# Load modules',
          'module load {}'.format(modules),
          '',
          '# Job commands',
          'bash {}'.format(commands),
          sep='\n',
          file=open(outfile, 'w'))
