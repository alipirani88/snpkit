#Creates and submits PBS scripts for performing BEAST analysis with desired xmls
#
#Takes as input:
#		1) config_file	- A 1-column file where each line is the file name (minus the .xml extension) for a beast run
# UPDATED 2017-09-21 from beast2/2.4.5 to beast2/2.4.7

#PRINT USAGE
if(@ARGV == 0)
{

	print "\n\nperl beast_PBS.pl config_file\n\n";

	exit;

}#end if


#READ IN COMMAND LINE ARGUMENTS
my $config_file = shift @ARGV;


#GO THROUGH EACH LINE OF THE CONFIG FILE, CREATING AND SUBMITTING PBS FILES
open CONFIG, $config_file;
my $count = 1;

foreach my $job (<CONFIG>)
{

	#GET CURRENT ARGUMENTS
	chomp $job;
	my ($fnameorig) = split /\t/, $job;
        my ($fname) = split /\./, $fnameorig;

	#CREATE OUTPUT DIRECTORY
	my $date = `date +%Y-%m-%d-%H-%M-%S`;
	chomp $date;

	my $wd = `pwd`;
	chomp $wd;
	$out_dir = "$wd/$date\_beast$count\_$fname";
	mkdir $out_dir;


	#CREATE AND WRITE PBS FILE
	my $rand = int(rand(10000000));
	my $pbs_file = "beast$count\_$fname\_$rand.pbs";

	open PBS, ">$pbs_file";

	print PBS "#!/bin/sh\n";
	print PBS "####  PBS preamble\n";
	print PBS "\n";
	print PBS "#PBS -N beast$count\_$fname\n";
	print PBS "\n";
	print PBS "# User info\n";
	print PBS "#PBS -M zenalapp\@umich.edu\n";
	print PBS "#PBS -m abe\n";
	print PBS "\n";
	print PBS "# Change the number of cores (ppn=1), amount of memory, and walltime:\n";
	print PBS "#PBS -l nodes=1:ppn=12,pmem=4gb,walltime=10:00:00:00\n";
	print PBS "#PBS -j oe\n";
	print PBS "#PBS -V\n";
	print PBS "\n";
	print PBS "#PBS -A esnitkin_flux\n";
	print PBS "#PBS -q flux\n";
	print PBS "#PBS -l qos=flux\n";
	print PBS "\n";
	print PBS "####  End PBS preamble\n";
	print PBS "\n";
	print PBS "#  Show list of CPUs you ran on, if you're running under PBS\n";
	print PBS "if [ -n \"\$PBS_NODEFILE\" ]; then cat \$PBS_NODEFILE; fi\n";
	print PBS "\n";
	print PBS "#  Change to the directory you submitted from\n";
	print PBS "cd $out_dir";
	print PBS "\n";
	print PBS "#Load modules\n";
	print PBS "module load beast2/2.4.7\n";
        print PBS "module load beagle\n";
        #print PBS "addonmanager -add BEAST_CLASSIC\n";
        #print PBS "addonmanager -add SCOTTI\n";
	print PBS "\n";
	print PBS "#  Put your job commands here:\n";
	print PBS "beast -beagle_SSE -instances 12 -threads 12 beast$count\_$fnameorig\n";
	print PBS "\n";

	close PBS;


	#SUBMIT JOB
	`mv $pbs_file $out_dir`;
        `cp $fnameorig $out_dir/beast$count\_$fnameorig`;
        chdir($out_dir);
        #`/nfs/esnitkin/Zena/lib/change_prefixes_beast_xml.sh beast$count\_$fnameorig`;
        #`mv beast$count\_$fname\_renamed.xml beast$count\_$fnameorig`;
	`qsub $pbs_file`;
        chdir($wd);
	$count ++;

}#end foreach
