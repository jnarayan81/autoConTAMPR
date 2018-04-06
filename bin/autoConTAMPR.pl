#!/usr/bin/env perl

=head1 NAME

autoConTAMPR.pl - autoConTAMPR script by Jitendra Narayan

=head1 SYNOPSIS

autoConTAMPR.pl --validate/-v

	autoConTAMPR.pl --validate/-v --conf/-c <configuration file>

autoConTAMPR.pl --plot/-p

	autoConTAMPR.pl --plot/-p --conf/-c <configuration file>

autoConTAMPR.pl --install/-i

	autoConTAMPR.pl --install/-i --conf/-c <configuration file>

autoConTAMPR.pl --full/-f

	autoConTAMPR.pl --full/-f --conf/-c <configuration file>

=cut

use strict;
use warnings;
use 5.010;

use Cwd;
use File::chdir;
use File::Copy;
use File::Temp qw(tempfile);
use File::Spec::Functions qw(rel2abs);
use File::Basename;
use FindBin '$Bin';
use lib::abs qw{./utils};
use File::Remove;
use File::Path qw(make_path remove_tree);
use Capture::Tiny ':all';
use Getopt::Long;
use Data::Dumper;
#use Statistics::R;
use File::Find;
use Pod::Usage;
use Spreadsheet::Read;
use autodie;
use File::Temp qw(tempfile);
use File::Copy;
use Carp; 
use Path::Tiny qw(path);
use POSIX;
use IO::File;
use Bio::SCF;
use Bio::DB::Fasta;
use File::Path qw(make_path remove_tree);
use Statistics::ChiSquare;
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;
use Term::ProgressBar 2.00;
use Bio::SeqIO;
use CPAN;
use lib "$FindBin::Bin/.";
require 'autoConTAMPR_module.pm';
#use lib 'utils/.';

#Basics mandatory autoConTAMPR variables
my (
$outfile, 	# Name for autoConTAMPR's main configuration file
$conf,
$validate,
$plot,
$full,
$install,
$help,		# If you are asking for help
);

# Default settings here for autoConTAMPR
my $current_version = "0.1";	#autoConTAMPR version
my %opts; my $nam;

print <<'WELCOME';

  --.---. .-. .---
	 / \
	/ ^ \
  --autoConTAMPR v0.1--

Citation - autoConTAMPR: Automatic Contingency Table Analysis of Minority Peak Ranks
License: Creative Commons Licence
Bug-reports and requests to: https://github.com/jnarayan81/autoConTAMPR/issues

Implemented by jitendra.narayan@unamur.be, jflot@ulb.ac.be, nicolas.debortoli@unamur.be and karine.vandoninck@unamur.be

WELCOME

printf "Operating system: %s\nPerl executable at: %s\n\n", $^O,  $^X, ;
if ($^O eq 'darwin') { print "Warning: Tested only on Linux !!\n\n"; } 

$|++;
#Get options for autoConTAMPR
GetOptions(
	\%opts,
	"conf|c=s" 	=> \$conf,
	"validate|v" 	=> \$validate,
	"plot|p" 	=> \$plot,
	"full|f" 	=> \$full,
	"install|i" 	=> \$install,
	"help|h" 	=> \$help,
);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please check manual.') unless ($validate or $plot or $install or $full);
validateHelp($current_version) if (($validate) and (!$conf));
installHelp($current_version) if (($install) and (!$conf));
fullHelp($current_version) if (($full) and (!$conf));
plotHelp($current_version) if (($plot) and (!$conf));


pod2usage(-msg => 'Please supply a valid filename.') unless ($conf && -s $conf);

# used to measure total execution time
my $start_time = time();

#Store thr opted name
if ($validate) {$nam = 'validate';} elsif ($plot) {$nam = 'plot';} elsif ($install) {$nam = 'install';} elsif ($full) {$nam = 'full';} else { print "Missing parameters\n"; exit(1);}

# Absolute path of the current working directory
my $autoConTAMPR_path = dirname(rel2abs($0)); #print " Path of the dir\n$autoConTAMPR_path --------\n";

# Parameters_ref - stores all user-defined parameters, such as file locations and program parameters
my $param_ref = read_config_files(\$conf, \$autoConTAMPR_path);  # stores the configuration in a hash reference

# Check all the parameters for their correctness
parameters_validator($param_ref); #checkin if user set the parameters right

#---------------------------------------
#Check and install the perl modules
if ($install) {
installModules($0, "$param_ref->{download}"); exit;
}

#---------------------------------------
if ($plot) { # add it later or ($full)
# Annotate the results
#Write the log files for all steps
openLog($nam);
openLogErr($nam);
my $workDir=getcwd;
my @files = glob("$workDir/utils/*.{txt,conf}");

if (!-e "$param_ref->{out_dir}/results") { print "It seems you forgot to run --validate or -v steps\n"; exit;}

if (!-e "$param_ref->{out_dir}/circos") { mkdir ("$param_ref->{out_dir}/circos") || die ("Couldn't create the directory with the circos of ambigram's analysis.\nDetails: $!\n"); }

for my $file (@files) {
    copy($file, "$param_ref->{out_dir}/circos") or die "Copy failed: $!";
}
	#Lets update the circos image location
	my @plotfiles = glob("$param_ref->{out_dir}/circos/*.conf");
	for my $fileN (@plotfiles) {
	my $filename = basename($fileN);
	if ($filename eq 'circos.image.conf') {
		my $data = read_file($fileN);
		my $location="$param_ref->{out_dir}/circos";
		$data =~ s/dir   = ./dir   = $location/g;
		write_file($fileN, $data);
		}
	}

#Plot the circos
my $circosRun = circosPlot ("$param_ref->{out_dir}/results/$param_ref->{result}.validate", "$param_ref->{table}", $param_ref);
if(! $circosRun) { print "Done with drawing the connections with circos :)\nCheck the circos folder in your out directory\n\n"; exit(1); }
else { print "It seems, circos path is not set correctly, please set the circos path in autoConTAMPR_config file and re-run!!\n - Terminated abnormally !!\n"; 
exit(1);}

} #Plot ends here

#---------------------------------------
# Delete the directory if already exist
if ((-e $param_ref->{out_dir}) and ($param_ref->{force} == 1)){ remove_tree( $param_ref->{out_dir});}

# Creating the needed directories if they don't exist
if (!-e $param_ref->{out_dir}) { mkdir ($param_ref->{out_dir}) || die ("Couldn't create the directory specified as '$param_ref->{out_dir}', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n"); }
else { die("Directory $param_ref->{out_dir} already exists.\n"); }

if (!-e "$param_ref->{out_dir}/results") { mkdir ("$param_ref->{out_dir}/results") || die ("Couldn't create the directory with the results of ambigram's analysis.\nDetails: $!\n"); }

#Copy file to the locations
copy($conf, "$param_ref->{out_dir}/project_config");

#Write the log files for all steps
open (RESULT, ">", "$param_ref->{out_dir}/results/$param_ref->{result}.$nam") || die ('Could not create result file. Please check writing permission in your current directory', "\n");
#Write the log files for all steps
open (LOG, ">", "$param_ref->{out_dir}/log.$nam") || die ('Could not create log file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
open (LOG_ERR, ">", "$param_ref->{out_dir}/log.err.$nam") || die ('Could not create log.err file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
open (SUMMARY, ">", "$param_ref->{out_dir}/results/$param_ref->{summary}.$nam") || die ('Could not create summary file. Please check writing permission in your current directory', "\n");


print "WORKING ON $nam $param_ref->{mode} mode\n";


if ($validate and ($param_ref->{mode} eq "general")) {
#my ($table, $consensus, $markerDir, $outfile, $draw, $brutal, $verbose, $logfile)=@_;
allConTAMPRLocal ("$param_ref->{table}", "$param_ref->{contaminated_seq}", "$param_ref->{data_dir}", "$param_ref->{out_dir}/results/$param_ref->{result}.$nam", "$param_ref->{draw}", "brutal", "no", "$param_ref->{out_dir}/log.$nam", $param_ref); # concated both species name in $SpsArray_ref

}

#---------------------------------------
if ($validate and ($param_ref->{mode} eq "selected")) {

#my ($table, $consensus, $markerDir, $outfile, $draw, $brutal, $verbose, $logfile)=@_;
if (!$param_ref->{table}) {print "It seems you forgot to provide me the TABLE of individual, please check readme for detail\n";}

splConTAMPR("$param_ref->{table}", "$param_ref->{contaminated_seq}", "$param_ref->{data_dir}", "$param_ref->{out_dir}/results/$param_ref->{result}.$nam", "$param_ref->{draw}", "brutal", "no", "$param_ref->{out_dir}/log.$nam", $param_ref); # concated both species name in $SpsArray_ref

}

#-------------------------------------------------------------------------------

#Lets clean the files
if ($param_ref->{clean} == 1) {  
	print "Working on to cleaning the temporary files ... \n";
	system ("rm -rf $param_ref->{out_dir}/results/plotData"); 
}

#print LOG "Analysis finished, closing autoConTAMPR.  \nGood bye :) \n" if $param_ref->{verbose};


close(LOG);
close(LOG_ERR);
close(RESULT);


############## Subs ###################

sub openLog {
my $nam =shift;
open (LOG, ">", "$param_ref->{out_dir}/log.$nam") || die ('Could not create log file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
}

sub openLogErr {
my $nam =shift;
open (LOG_ERR, ">", "$param_ref->{out_dir}/log.err.$nam") || die ('Could not create log.err file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
}

sub openResult {
my $nam =shift;
open (RESULT, ">", "$param_ref->{out_dir}/results/$param_ref->{result}.$nam") || die ('Could not create result file. Please check writing permission in your current directory', "\n");
}

__END__
