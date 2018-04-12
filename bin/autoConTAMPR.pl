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

=cut

use strict;
use warnings;

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
use Scalar::Util qw(looks_like_number);
use IO::File;
use Bio::SCF;
use Bio::DB::Fasta;
use Bio::SeqIO;
use File::Path qw(make_path remove_tree);
use Statistics::ChiSquare;
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;
use Term::ProgressBar 2.00;
use CPAN;
use lib "$FindBin::Bin/.";
require 'autoConTAMPR_module.pm';
#use lib 'utils/.';

#Basics mandatory autoConTAMPR variables
my (
$conf,
$validate,
$plot,
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

local $SIG{__WARN__} = sub {
    my $message = shift;
    logger($message);
  };

#Get options for autoConTAMPR
GetOptions(
	\%opts,
	"conf|c=s" 	=> \$conf,
	"validate|v" 	=> \$validate,
	"plot|p" 	=> \$plot,
	"install|i" 	=> \$install,
	"help|h" 	=> \$help,
);
pod2usage(-verbose => 1) if ($help);
pod2usage(-msg => 'Please check manual.') unless ($validate or $plot or $install);
validateHelp($current_version) if (($validate) and (!$conf));
installHelp($current_version) if (($install) and (!$conf));
plotHelp($current_version) if (($plot) and (!$conf));


pod2usage(-msg => 'Please supply a valid filename.') unless ($conf && -s $conf);

# used to measure total execution time
my $start_time = time();

#Store thr opted name
if ($validate) {$nam = 'validate';} elsif ($plot) {$nam = 'plot';} elsif ($install) {$nam = 'install';} else { print "Missing parameters\n"; exit(1);}

# Absolute path of the current working directory
my $autoConTAMPR_path = dirname(rel2abs($0)); #print " Path of the dir\n$autoConTAMPR_path --------\n";

# Parameters_ref - stores all user-defined parameters, such as file locations and program parameters
my $param_ref = read_config_files(\$conf, \$autoConTAMPR_path);  # stores the configuration in a hash reference

# Check all the parameters for their correctness
if ($param_ref->{check} == 1) {
parameters_validator($param_ref); #checkin if user set the parameters right
}

#---------------------------------------
#Check and install the perl modules
if ($install) {
installModules($0, "$param_ref->{download}"); exit;
}

#---------------------------------------
if ($plot) {
if (!-e "$param_ref->{out_dir}/results") { print "It seems you forgot to run --validate or -v steps\n"; exit;}
print "Generating circos files.\n";
# Annotate the results
#Write the log files for all steps
openLog($nam);
openLogErr($nam);
my $workDir=getcwd;
my @files = glob("$workDir/utils/*.{txt,conf}");

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
#check the status
if(! $circosRun) { print "Done with drawing the connections with circos :)\nCheck the circos folder in your out directory\n\n";}
else { print "It seems, circos path is not set correctly, please set the circos path in autoConTAMPR_config file and re-run!!\n - Terminated abnormally !!\n"; 
exit(1);}

#Lets plot separately
	if ($param_ref->{separate}) { 
		# Delete the directory if already exist
		if ((-e "$param_ref->{out_dir}/circos/separate") and ($param_ref->{force} == 1)){ remove_tree( "$param_ref->{out_dir}/circos/separate");}
		# Creating the needed directories if they don't exist
		if (!-e "$param_ref->{out_dir}/circos/separate") { mkdir ("$param_ref->{out_dir}/circos/separate") || die ("Couldn't create the directory specified as '$param_ref->{out_dir}/circos/separate', check if you are trying to create a subdirectory in a non-existent directory or if you don't have permission to create a directory there.\n"); }
		#Store into hash
		my $hash_ref = storeIntoHash ($param_ref->{table}, 3);
		my %hash = %$hash_ref; my $lastKey="NA";
		
		foreach my $key (keys %hash) {
			my $data = read_file("$param_ref->{out_dir}/circos/circos.image.conf");
			my $location="$key.png";
			my $n = "autoConTMPRcircos.png"; my $m = "$lastKey.png";
			$data =~ s/$m|$n/$location/g;
			write_file("$param_ref->{out_dir}/circos/circos.image.conf", $data);
			extractSpsLine("$param_ref->{out_dir}/results/$param_ref->{result}.validate","$param_ref->{out_dir}/circos/result_$key.tmp", $key);
			circosPlot ("$param_ref->{out_dir}/circos/result_$key.tmp", "$param_ref->{table}", $param_ref);
			move("$param_ref->{out_dir}/circos/$key.svg","$param_ref->{out_dir}/circos/separate/"); 
			move("$param_ref->{out_dir}/circos/$key.png","$param_ref->{out_dir}/circos/separate/"); 
			$lastKey = $key;
		unlink "$param_ref->{out_dir}/circos/result_$key.tmp" or warn "Can't unlink $key $!\n"; 
		}
	}
exit(1); #stop the plot script here
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
my ($WARN, $LOG, $LOG_ERR, $RESULT, $SUMMARY);
open ($WARN, ">", "$param_ref->{out_dir}/warn.$nam") || die ('Could not create log file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
open ($LOG, ">", "$param_ref->{out_dir}/log.$nam") || die ('Could not create log file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
open ($LOG_ERR, ">", "$param_ref->{out_dir}/log.err.$nam") || die ('Could not create log.err file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");

print "$nam the chromatogram files using autoConTAMPR $param_ref->{mode} mode\n";


if ($validate and ($param_ref->{mode} eq "general")) {
#my ($table, $consensus, $markerDir, $outfile, $draw, $brutal, $verbose, $logfile)=@_;
allConTAMPRLocal ("$param_ref->{table}", "$param_ref->{consensus}", "$param_ref->{data_dir}", "$param_ref->{out_dir}/results/$param_ref->{result}.$nam", "$param_ref->{draw}", "brutal", "no", $LOG, $param_ref, $WARN, $LOG_ERR); # concated both species name in $SpsArray_ref

}

#---------------------------------------
if ($validate and ($param_ref->{mode} eq "selected")) {
print "Thanks for opting selected mode, we are working on it. Try later ...\n"; exit;
#my ($table, $consensus, $markerDir, $outfile, $draw, $brutal, $verbose, $logfile)=@_;
if (!$param_ref->{table}) {print "It seems you forgot to provide me the TABLE of individual, please check readme for detail\n";}

splConTAMPR("$param_ref->{table}", "$param_ref->{consensus}", "$param_ref->{data_dir}", "$param_ref->{out_dir}/results/$param_ref->{result}.$nam", "$param_ref->{draw}", "brutal", "no", $LOG, $param_ref, $WARN, $LOG_ERR); # concated both species name in $SpsArray_ref
}
#-------------------------------------------------------------------------------
#Lets clean the files
if ($param_ref->{clean} == 1) {  
	print $LOG "Working on to cleaning the temporary files ... \n";
	#system ("rm -rf $param_ref->{out_dir}/results/plotData"); 
}

print $LOG "Analysis finished, closing autoConTAMPR.\nGood bye :)\n" if $param_ref->{verbose};

close ($WARN);
close($LOG);
close($LOG_ERR);


############## Subs ###################
#log them
sub logger {
    my ($msg) = @_;
        print $WARN "$msg";
  }

sub openLog {
my $nam =shift;
open ($LOG, ">", "$param_ref->{out_dir}/log.$nam") || die ('Could not create log file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
}

sub openLogErr {
my $nam =shift;
open ($LOG_ERR, ">", "$param_ref->{out_dir}/log.err.$nam") || die ('Could not create log.err file in ', $param_ref->{out_dir}, '. Please check writing permission in your current directory', "\n");
}

sub openResult {
my $nam =shift;
open ($RESULT, ">", "$param_ref->{out_dir}/results/$param_ref->{result}.$nam") || die ('Could not create result file. Please check writing permission in your current directory', "\n");
}

__END__
#Convert to lowercase
awk '{ if ($0 !~ />/) {print toupper($0)} else {print $0} }' name.fasta

tr A-Z a-z < input 


