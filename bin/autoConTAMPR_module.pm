use strict;
use warnings;

use re qw(eval);
use vars qw($matchStart);

#Check perl modules
sub installModules {
my ($script, $download) = @_;
die "Provide script file name on the command line\n" unless defined $script;

my @modules = qw(
Cwd;
File::chdir;
File::Copy;
File::Temp;
File::Spec::Functions;
File::Basename;
FindBin;
lib::abs;
File::Remove;
File::Path;
Capture::Tiny;
Getopt::Long;
Data::Dumper;
File::Find;
Pod::Usage;
Spreadsheet::Read;
autodie;
File::Temp;
File::Copy;
Carp; 
Path::Tiny;
POSIX;
Scalar::Util;
IO::File;
Bio::SCF;
Bio::DB::Fasta;
File::Path;
Statistics::ChiSquare;
List::Util;
Statistics::Distributions;
Term::ProgressBar 2.00;
Bio::SeqIO;
CPAN;
);

print "Checking mandatory modules for autoConTAMPR\n\n";
checkModules(\@modules, $download); ## Provide the name of all modules

=pod	
until ( do $script ) {
    my $expire = $@;
    if ( my ($file) = $expire =~ /^Can't locate (.+?) in/ ) {
        my $module = $file;
        $module =~ s/\.(\w+)$//;
        $module = join('::', split '/', $module);
        print "Attempting to install '$module' via cpan\n";
        system(cpan => $module);
        last unless prompt(y => 'Try Again?', '', 'n');
    }
    else {
        die $expire;
    }
}
=cut

}

#Check the modules and if not install try to install it
sub checkModules {
	my ($reqMod_ref, $download)=@_;
	my @reqMod=@$reqMod_ref;
	for(@reqMod) {
    	eval "use $_";
    		if ($@) {
        		warn "Not found -> $_ module\t Need to install $_ \n" if $@;
			if ($download eq '1') { print "Trying to download using CPAN\n"; CPAN::Shell->install($_); }
    		} else {
        		print "Found : $_\n";
    		}
	}
}


#All individual against all
sub allConTAMPRLocal {
my ($table, $consensus, $markerDir, $outfile, $draw, $brutal, $verbose, $logfh, $param, $wfh, $logerr)=@_;

#my $logfh = IO::File->new("$logfile",'w') if $verbose;
my $ofh = IO::File->new("$outfile",'w') if $verbose;
my $seperator='*' x 100;

#Store into hash
my $hash_ref = storeIntoHash ($table, 1);
my %hash = %$hash_ref;

my $ss = ReadData ("$table", attr => 1)->[1]; #Reading sheet 1 
my $book = Spreadsheet::Read->new ("$table");
my $sheet = $book->sheet (1);     # OO
my $progress = Term::ProgressBar->new($ss->{maxrow});
foreach my $row (1 .. $ss->{maxrow}) {
	next if $row == 1; #avoid row
	my $speciesName; my $rVal;
	print $logfh "$seperator\nChecking on row number $row\n$seperator\n" if $verbose;
    	foreach my $col (1 .. $ss->{maxcol}) {
        	my $cell = cr2cell ($col, $row);
        	#printf "%s %-3s %d  ", $cell, $ss->{$cell}, $ss->{attr}[$col][$row]{merged};
		if ($col == 1) { $speciesName=$ss->{$cell}; next;}		
		#print "$cell\t";
		if ($col == 2) { $rVal=$ss->{$cell};}
		next if $col > 2;
		my @col = $sheet->column($col);  #foreach (@col) { print "$_\n";} #exit; # Column 2 is indivisual name
		shift @col;
		if (index($ss->{$cell}, "($speciesName)") != -1) { print $logfh "'$ss->{$cell}' contains '$speciesName'\n" if $verbose; }

		foreach my $cVal (@col) {
			#foreach my $cVal(@col) {
			#next if $col == 2; #ignore as this is indivisual name
			#next if !$ss->{$cell}; #Ignore if empty	
			#print $logfh "'$ss->{$cell}' NOT contains species '$speciesName\t$cell'\n" if $verbose;
			#my ($markerName, $mCol, $mRow)=findMarker($ss, $speciesName, $ss->{$cell}, $col);
			#print $logfh "$mCol,$mRow\t$ss->{$cell}\t$markerName --->\n" if $verbose;
			#my $indCol=1; #column with ind name !!!!!!!!!
			#my ($indRName)=findName($ss, $row, $indCol);
			#my ($indCName)=findName($ss, $mRow, $indCol);

			my ($indRName)="_"."$rVal"."_";
			my ($indCName)=$cVal;

			#print "$indRName\t$indCName\n" if $verbose;
			my ($fileArray_ref) = returnFilePath($markerDir,$indRName, $wfh);
			#my $workDir=getcwd(); 
			my @fArray= @$fileArray_ref;
			my $lHand=''; my $hHand=''; my $direction='NA';
			
			#Ignore if not scf file found
			next if (!@fArray) and ($param->{ignore} == 1);
			
			#Lets work on both scf files
			my $tSec=1; my $tThi=1; my $tFou=1;
			foreach my $f (sort {$a cmp $b} @fArray) {
				if ($f =~ m/^L/)  { print $logfh "$f\n" if $verbose; $direction='right';}
				if ($f =~ m/^H/) {print $logfh "$f\n" if $verbose; $direction='left';}
				print $logfh "$f\t$markerDir\n" if $verbose;
				my $cSeq = extractSeq($consensus,$indCName, $direction);
				my ($chroId, $chroSeq)=extractChroSeq ("$markerDir/$f", "$direction");
				my @scfLen= split '\_', $chroId;
				print $logfh "Chromatogram:$chroId\n$chroSeq\n\n" if $verbose; 
				print $logfh "Consensus:$indCName\n$cSeq\n\nAlignment\n" if $verbose; 
				my ($cAln, $chroAln, $sti, $stj) = alnSW ($cSeq, $chroSeq, $verbose, $logfh, $scfLen[-1], $chroId, $indCName);
				my ($sec, $thi, $fou) = checkIntensityLocal ($cAln, $chroAln, "$markerDir/$f", $sti, $stj, $verbose, $logfh, $indCName); #$sti is scf and $stj is consensus
				$tSec=$tSec+$sec; $tThi=$tThi+$thi; $tFou=$tFou+$fou; 
			}
			my @allInt = ($tSec, $tThi, $tFou);
			#check the randomness
			my $normal= (sum( @allInt)/3);
			my $chi = chisquare( @allInt);
			my $chiVal = chi_quared( observed => [ @allInt ], expected => [ ($normal) x 3 ]); 
			#And that prints out 0.018360, the probability of that set of freq occurring by chance.
			my $chisprob=Statistics::Distributions::fprob(2,3,$chi);
			my $chichi='NA';
			$tSec= $tSec-1; $tThi= $tThi-1; $tFou= $tFou-1; # Just to count 1 less (because cnt start from 1)
			my $spsNameConnection=$hash{$ss->{$cell}};
			$indRName =~ s/\_//g;
			print $ofh "$speciesName\t$indRName\t$indCName\t$spsNameConnection\t$tSec\t$tThi\t$tFou\t$chiVal\n";

			if ($brutal eq "yes") { print $logfh "\nExiting with ONE run --- \n\n" if $verbose; exit; }	
		}
	
	}
	$progress->update($row)
     	#print "\n";
}
close $logfh;
close $ofh;
#Reformat the result data
reformatOut($outfile, "$param->{out_dir}/results/plotData");
my $workDir = getcwd();
	if ($draw) {
	print $logfh "\nYou have opted to draw final table- Drawing ...\n";
	my $rDrawer = system "Rscript $workDir/utils/conPlot.R $param->{out_dir}/results/plotData $param->{out_dir}/results/seeStat.pdf";
		if(! $rDrawer) { print $logfh "Done with drawing final table with R :)\n\n"; exit(1); }
		else { print "\nIt seems, R path is not set correctly, please set the R path in autoConTAMPR_config file and re-run!!\n - Terminated final R drawing abnormally !!\n"; exit; }
	}
}

#Read the file
sub read_file {
my ($filename) = @_;
open my $in, '<:encoding(UTF-8)', $filename or die "Could not open '$filename' for reading $!";
local $/ = undef;
my $all = <$in>;
close $in;
return $all;
}


#store into hash
sub storeIntoHash {
my ($table, $col)=@_;
my %hash;
open(my $fh, '<:encoding(UTF-8)', $table) or die "Could not open file '$table' $!";
while (my $row = <$fh>) {
  chomp $row;
  my @tmp = split '\t', $row;
  if ($col == 1) {
  $hash{$tmp[1]}=$tmp[0];
  }
  elsif($col == 3) {
  $hash{$tmp[0]}=$tmp[3];
  }
}
return \%hash;
}

#Write to the file
sub write_file {
my ($filename, $content) = @_;
open my $out, '>:encoding(UTF-8)', $filename or die "Could not open '$filename' for writing $!";;
print $out $content;
close $out;
return;
}


#All individual against all --not in use just for testing
sub allConTAMPRGlobal {
my ($table, $consensus, $markerDir, $outfile, $plot, $brutal, $verbose, $logfile)=@_;

my $logfh = IO::File->new("$logfile",'w') if $verbose;
my $ofh = IO::File->new("$outfile",'w') if $verbose;
my $seperator='*' x 100;

my $ss = ReadData ("$table", attr => 1)->[1]; #Reading sheet 1 
my $book = Spreadsheet::Read->new ("$table");
my $sheet = $book->sheet (1);     # OO
my $progress = Term::ProgressBar->new($ss->{maxrow});
foreach my $row (1 .. $ss->{maxrow}) {
	next if $row == 1; #avoid row
	my $speciesName; my $rVal;
	print $logfh "$seperator\nChecking on row number $row\n$seperator\n" if $verbose;
    	foreach my $col (1 .. $ss->{maxcol}) {
        	my $cell = cr2cell ($col, $row);
        	#printf "%s %-3s %d  ", $cell, $ss->{$cell}, $ss->{attr}[$col][$row]{merged};
		if ($col == 1) { $speciesName=$ss->{$cell}; next;}		
		#print "$cell\t";
		if ($col == 2) { $rVal=$ss->{$cell}; print "$rVal ->>\n";}
		next if $col > 2;
		my @col = $sheet->column($col); foreach (@col) { print "$_ <<-\n";} #exit; # Column 2 is indivisual name
		shift @col;
		if (index($ss->{$cell}, "($speciesName)") != -1) { print $logfh "'$ss->{$cell}' contains '$speciesName'\n" if $verbose; }
		foreach my $cVal (@col) {
			my ($indRName)="_"."$rVal"."_";
			my ($indCName)=$cVal;
			print "$indRName\t$indCName ++ \n" if $verbose;
			my ($fileArray_ref) = returnFilePath($markerDir,$indRName);
			my @fArray= @$fileArray_ref;
			my $lHand=''; my $hHand=''; my $direction='NA';
			#Lets work on both scf files
			my $tSec=1; my $tThi=1; my $tFou=1;
			foreach my $f (sort {$a cmp $b} @fArray) {
				if ($f =~ m/^L/)  { print $logfh "$f\n" if $verbose; $direction='right';}
				if ($f =~ m/^H/) {print $logfh "$f\n" if $verbose; $direction='left';}
				print $logfh "$f\t$markerDir/$f ++++\n" if $verbose;
				my $cSeq = extractSeq($consensus,$indCName, $direction);
				my ($chroId, $chroSeq)=extractChroSeq ("$markerDir/$f", "$direction");
				print $logfh "Chromatogram:$chroId\n$chroSeq\n\n" if $verbose; 
				print $logfh "Consensus:$indCName\n$cSeq\n\nAlignment\n" if $verbose; 
				my ($cAln, $chroAln) = alnNW ($cSeq, $chroSeq, $verbose, $logfh);
				my ($sec, $thi, $fou) = checkIntensity ($cAln, $chroAln, "$markerDir/$f", $verbose, $logfh);
				$tSec=$tSec+$sec; $tThi=$tThi+$thi; $tFou=$tFou+$fou; 
			}
			my @allInt = ($tSec, $tThi, $tFou);
			#check the randomness
			my $normal= (sum( @allInt)/3);
			my $chi = chisquare( @allInt);
			my $chiVal = chi_quared( observed => [ @allInt ], expected => [ ($normal) x 3 ]); 
			#And that prints out 0.018360, the probability of that set of freq occurring by chance.
			my $chisprob=Statistics::Distributions::fprob(2,3,$chi);
			my $chichi='NA';
			$tSec= $tSec-1; $tThi= $tThi-1; $tFou= $tFou-1; # Just to count 1 less (because cnt start from 1)
			print $ofh "$speciesName\t$indRName\t$indCName\t$ss->{$cell}\t$tSec\t$tThi\t$tFou\t$chiVal\t$chi\t$chisprob\t$chichi\n";

			if ($brutal eq "yes") { print $logfh "\nExiting with ONE run --- \n\n" if $verbose; exit; }	
		}
	
	}
	$progress->update($row)
     	#print "\n";
}
close $logfh;
close $ofh;

reformatOut($outfile, 'plotData');
my $workDir = getcwd();
if ($plot) { system "Rscript $workDir/utils/conPlot.R plotData";}

}

# Read config files in the form element = value #comment --->
sub read_config_files {
  my $project_config_file = shift;
  my $autoConTAMPR_path = shift;
  my %param;
  #There is two config files, one for general settings and other for third party software
  open(my $user_config, "<", "$$project_config_file") || die ("Couldn't open the project configuration file: $!\n");
  open(my $autoConTAMPR_config, "<", "$$autoConTAMPR_path/../config_files/autoConTAMPR_config") || die ("The configuration file 'autoConTAMPR_config' couldn't be read, please check if the file exists and if its permissions are properly set.\n");

# INPUT FILES --->
  $param{data_dir} = read_config('data_dir', $param{data_dir}, $user_config);

# PROJECT NAME --->
  $param{out_dir} = read_config('out_dir', $param{out_dir}, $user_config);

# PROJECT CONFIGURATION --->
  $param{consensus} = read_config('consensus', '', $user_config);
  $param{verbose} = read_config('verbose', '', $user_config);
  $param{force} = read_config('force', '', $user_config);
  $param{check} = read_config('check', '', $user_config);
  $param{draw} = read_config('draw', '', $user_config);
  $param{clean} = read_config('clean', '', $user_config);
  $param{separate} = read_config('separate', '', $user_config);
  $param{ignore} = read_config('ignore', '', $user_config);

  $param{download} = read_config('download', '', $user_config);
  $param{table} = read_config('table', '', $user_config);

  $param{threshold} = read_config('threshold', '', $user_config);
  $param{len} = read_config('len', '', $user_config);

#GENERAL SETTINGS
  $param{mode} = read_config('mode', '', $user_config);

  close($user_config);

# PATH TO EXTERNAL PROGRAMS --->

  $param{circos} = read_config('circos', $param{reads_dir}, $autoConTAMPR_config);
  $param{rscript} = read_config('rscript', $param{reads_dir}, $autoConTAMPR_config);

# OUTPUT NAMES --->
  $param{result} = read_config('result', '', $autoConTAMPR_config);
  $param{summary} = read_config('summary', '', $autoConTAMPR_config);
  $param{warnings} = read_config('warnings', '', $autoConTAMPR_config);

  close($autoConTAMPR_config);
  return \%param;
}

############################################################################################""
sub read_config { # file format element = value
  my ($parameter, $reads_dir, $config_file) = @_;

  seek($config_file, 0, 0);              # added to allow any order of parameters in the config files, preventing unfriendly error messages if the user changes the order
  while (my $line = <$config_file>){
    if ($line =~ /^\s*$parameter\s*=/) {    # the string to be searched in the file
      chomp ($line);
      $line =~ s/^\s*$parameter\s*=\s*//;   # removing what comes before the user input
      $line =~ s/#.*$//;                    # removing what comes after the user input (commentaries)
      $line =~ s/\s*$//;                    # removing what comes after the user input (space caracteres)
      $line =~ s/\$data_dir/$reads_dir/;     # allows the use of "$reads_dir" in the config file as a reference to the said parameter
      if ($line eq 'undef' || $line eq '') { return; }
      else { return $line; }
    }
  }
  return;
}

############################################################################################""
# function to identify errors in the configuration files and direct the user to the needed adjustments
sub parameters_validator { #check for all parameters,
  my $param = shift;

  my $config_path = getcwd();
  $config_path =~ s/\/\w+$/\/config_files/;

# INPUT FILES --->
  if (!defined $param->{data_dir}) { die ("No path to the nucleotide files was specified in your project's configuration file, please fill the parameter 'data_dir'.\n"); }
  if (!-d $param->{data_dir}) { die ("The path to your project's nucleotide files isn't a valid directory, please check if the path in 'data_dir' is correct: $param->{data_dir}\n"); }
  if (!-r $param->{data_dir}) { die ("You don't have permission to read in your project's nucleotide directory, please redefine your permissions.\n"); }


# PROJECT CONFIGURATION --->

  if (!defined $param->{verbose}) { $param->{verbose} = 1; } # default value
  if (!defined $param->{force}) { $param->{force} = 1; } # default value
  if (!defined $param->{check}) { $param->{force} = 1; } # default value
  if (!defined $param->{draw}) { $param->{force} = 1; } # default value
  if (!defined $param->{clean}) { $param->{force} = 1; } # default value
  if (!defined $param->{separate}) { $param->{force} = 1; } # default value
  if (!defined $param->{ignore}) { $param->{force} = 1; } # default value
  if (!defined $param->{download}) { $param->{force} = 1; } # default value
  if (!defined $param->{threshold}) { $param->{force} = 0.05; } # default value
  if (!defined $param->{len}) { $param->{force} = 600; } # default value
  if (!defined $param->{force}) { $param->{force} = 1; } # default value

  if (! looks_like_number($param->{verbose})) { print "$param->{verbose} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{force})) { print "$param->{force} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{check})) { print "$param->{check} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{draw})) { print "$param->{draw} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{clean})) { print "$param->{clean} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{separate})) { print "$param->{separate} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{ignore})) { print "$param->{ignore} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{download})) { print "$param->{download} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{threshold})) { print "$param->{threshold} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{len})) { print "$param->{len} should be number in your config file\n"; exit;}
  if (! looks_like_number($param->{force})) { print "$param->{force} should be number in your config file\n"; exit;}

  if (!defined $param->{out_dir}) {die "Project directory not configured. Please set out_dir element in configuration file\n";}

}

# now the script loads all nucleotide sequence files to a hash structure,
# checks their validity and translates them to protein sequence
sub parse_genome {
  my ($param) = shift;
  opendir (my $nt_files_dir, $param->{data_dir}) || die ("Path to asseembly fasta files not found: $!\n");
  my (%sequence_data);
  print LOG ('Parsing overlapsd genome/contigs/scaffolds files', "\n") if $param->{verbose};
  my $id_numeric_component = 1;  # Generate unique IDs with each sequence later on
  while (my $file = readdir ($nt_files_dir)) {
    if (($file eq '.') || ($file eq '..') || ($file =~ /^\./) || ($file =~ /~$/)) { next; }  # Prevents from reading hidden or backup files
    my $file_content = new Bio::SeqIO(-format => 'fasta',-file => "$param->{data_dir}/$file");
    print LOG ('Reading file ', $file, "\n") if $param->{verbose};
    while (my $gene_info = $file_content->next_seq()) {
      my $sequence = $gene_info->seq();
      my $len = length ($sequence);
      my $accession_number = $gene_info->display_id;
      $sequence_data{$accession_number}{status} = "OK"; #everybody starts fine
      $sequence_data{$accession_number}{problem_desc} = "-"; #everybody starts fine
      if ($sequence_data{$accession_number}{status} eq "OK") { # Add check points here <<<<<<
        $sequence_data{$accession_number}{nuc_seq} = $sequence;
	$sequence_data{$accession_number}{len} = $len;
      }
    }
  }
  print LOG ('Done', "\n") if $param->{verbose};
  closedir ($nt_files_dir);
  return (\%sequence_data);
}

sub parse_gene_id {
  my @aux = split (/\(/, $_[0]);
  my $specie = $aux[1];
  $specie =~ s/\)//g;
  return ($aux[0], $specie);  # aux[0] has the if o the gene
}

sub mean {
  my @tmp = @{$_[0]};
  my $soma = 0;
  foreach my $value(@tmp) {
    $soma = $soma + $value;
  }
  my $mean = ($soma/($#tmp+1));
  return $mean;
}

############################################################################################""
# Prints the time taken by the tasks of the group before the codeml runs in a file. Used for summary
sub print_task_time {
  my ($ortholog_dir, $task_time, $name) = @_;
  my $f_tree_time = time() - $$task_time;
  open(my $fh_time_write, ">", "$$ortholog_dir/$name");
  print $fh_time_write ('Time taken by task: ', $f_tree_time, "\n");
#  print "$$fh_time_write\n";
  close ($fh_time_write);
  return;
}

############################################################################################""
sub write_summary {
  my ($param, $clusters_ref, $start_time) = @_;
  my ($f_tree, $model1, $model2, $model7, $model8) = (0,0,0,0,0);


  # printing time spent by the program to run
  my $total_time = time() - $$start_time;
  my ($hours, $minutes, $seconds) = divide_time(\$total_time);

  print SUMMARY ("Time spent: $hours:$minutes:$seconds ($total_time seconds)\n");
  foreach my $ortholog_group (keys %{$clusters_ref}) {
    if (-s "$param->{out_dir}/intermediate_files/$ortholog_group/time_id_rec") {
      open (my $fh_time_read, "<", "$param->{out_dir}/intermediate_files/$ortholog_group/time_id_rec");
      my $line = <$fh_time_read>;
      $f_tree += $1 if ($line =~ /:\s(\d+)/);
      close($fh_time_read);
    }
  }

  my $sequential_time = $f_tree+$model1+$model2+$model7+$model8;

  ($hours, $minutes, $seconds) = divide_time(\$sequential_time);
  print SUMMARY ("Total time (sequential run): $hours:$minutes:$seconds ($sequential_time seconds)\n");

  ($hours, $minutes, $seconds) = divide_time(\$f_tree);
  print SUMMARY (" - total time on building phylogenetic trees: $hours:$minutes:$seconds ($f_tree seconds)\n");

  return;
}

############################################################################################""
#Time check
sub divide_time {
  my $total_time = shift;

  my $hours = POSIX::floor( $$total_time / 3600 );
  my $minutes = POSIX::floor(($$total_time % 3600) / 60);
  if ($minutes < 10) { $minutes = '0' . $minutes; }
  my $seconds = $$total_time % 60;
  if ($seconds < 10) { $seconds = '0' . $seconds; }

  return ($hours, $minutes, $seconds);
}

############################################################################################""!!
#Function to move the files with wild
sub moveFiles {
    my ( $source_ref, $arc_dir ) = @_;
    my @old_files = @$source_ref;
    foreach my $old_file (@old_files)
         {
    #my ($short_file_name) = $old_file =~ m~/(.*?\.dat)$~;
    #my $new_file = $arc_dir . $short_file_name;
    move($old_file, $arc_dir) or die "Could not move $old_file to $arc_dir: $!\n";
   }
}

=pod

Comment here

=cut



########################################################################################"
#Store fasta file to hash
sub fastafile2hash {
    my $fastafile = shift @_;
    my %sequences;
    my $fh = &read_fh($fastafile);
    my $seqid;
    while (my $line = <$fh>) {
        if ($line =~ /^>(\S+)(.*)/) {
            $seqid = $1;
            $sequences{$seqid}{desc} = $2;
        }
        else {
            chomp $line;
            $sequences{$seqid}{seq}     .= $line;
            $sequences{$seqid}{len}     += length $line;
            $sequences{$seqid}{gc}      += ($line =~ tr/gcGC/gcGC/);
            $line =~ s/[^atgc]/N/ig;
            $sequences{$seqid}{nonatgc} += ($line =~ tr/N/N/);
        }
    }
    close $fh;
    return \%sequences;
}


########################################################################################"
#Open and Read a file
sub read_fh {
    my $filename = shift @_;
    my $filehandle;
    if ($filename =~ /gz$/) {
        open $filehandle, "gunzip -dc $filename |" or die $!;
    }
    else {
        open $filehandle, "<$filename" or die $!;
    }
    return $filehandle;
}


sub extractSeq_new {
my ($name, $st, $end, $db)=@_;
my $seq = $db->seq($name, $st => $end);
return $seq;
}


########################################################################################"
#Process the GCAT sequence
sub processGCAT {
    my $sequence = shift;
    my @letters = split(//, $sequence);
    my $gccount = 0; my $totalcount = 0; my $gccontent = 0;
    my $acount = 0; my $tcount = 0; my $gcount = 0; my $ccount = 0; my $atcontent =0;
    foreach my $i (@letters) {
	if (lc($i) =~ /[a-z]/) { $totalcount++;}
	if (lc($i) eq "g" || lc($i) eq "c") { $gccount++; }
	if (lc($i) eq "a") { $acount++;}
	if (lc($i) eq "t") { $tcount++;}
	if (lc($i) eq "g") { $gcount++;}
	if (lc($i) eq "c") { $ccount++;}
    }
    if ($totalcount > 0) {
	$gccontent = (100 * $gccount) / $totalcount;
    }
    else {
	$gccontent = 0;
    }
    my $others=($totalcount-($acount+$tcount+$gcount+$ccount));
    return ($gccontent,$others,$totalcount,$gcount,$ccount,$acount,$tcount);

}


########################################################################################"
#Print the hash values
sub print_hash_final {
    my ($href,$fhandler)  = @_;
    while( my( $key, $val ) = each %{$href} ) {
        print $fhandler "$key\n";
	#print $fhandler "$key\t=>$val\n";
    }
}

########################################################################################"
#Mean of GCAT
use List::Util qw(sum);
sub meanGCAT { return @_ ? sum(@_) / @_ : 0 }
#sub meanGCAT { return sum(@_)/@_; }


########################################################################################"
#  Sigmoid function
sub sigmoidFun {
my $val = shift;
   my ($h) = @_;
   return 1.0 / ( 1.0 + exp(-$val) );
}


########################################################################################"

sub sumArray {
    return( defined $_[0] ? $_[0] + sumArray(@_[1..$#_]) : 0 );
}

########################################################################################"

sub spacer { return (" " x 20000);}


########################################################################################"

=pod
sub round {

    my ($nr,$decimals) = @_;
    return (-1)*(int(abs($nr)*(10**$decimals) +.5 ) / (10**$decimals)) if $nr<0;
    return int( $nr*(10**$decimals) +.5 ) / (10**$decimals);

}
=cut


sub checkFreq {
my ($str, $letter, $mutation, $param) = @_;
my @all; my $finalScore=0; my $longString=0;
while ($str =~ /($letter+)/g) {
push @all, length($1);
	@all = grep { $_ != 1 } @all; #Keep all except 1 (i.e delete 1)

	if (@all) {
	my $max = (sort { $b <=> $a } @all)[0];
	if ($max > $param->{consider}) { $longString=$max/$param->{consider}}
	}

	my $scoreGen = scalar @all/4; #4 is reads
	if ($scoreGen == 1) { $mutation = 0;}
	my $scorePen = $scoreGen + ($mutation * $param->{mutation}) + $longString;
	my $scoreFin = $scorePen;
	$finalScore=sigmoidFun ($scoreFin);
}
$str=join ',', @all;
return "$str\t$finalScore";
}

########################################################################################"
sub get_genome_sequence {
   my ($fpath) = @_;
   open GENOME, "<$fpath" or die "Can't open $fpath: $!\n";
   $_ = <GENOME>; # discard first line
   my $sequence = "";
   while (<GENOME>) {
      chomp($_);                  # Remove line break
      s/\r//;                     # And carriage return
      $_ = uc $_;                 # Make the line upper-case;
      $sequence = $sequence . $_; # Append the line to the sequence
   }
   return $sequence;
}

########################################################################################"
#Get reverse complement
sub rcomplement {
        my $dna = shift;
	# reverse the DNA sequence
	my $revcomp = reverse($dna);
	# complement the reversed DNA sequence
        $revcomp =~ tr/ACGTacgt/TGCAtgca/;
        return $revcomp;
}
#print the lines

########################################################################################"
use Carp qw/croak/;
sub file_write {
        my $file = shift;
        open IO, ">$file" or croak "Cannot open $file for output: $!\n";
        print IO @_;
        close IO;
}

########################################################################################"
#Open and write a file
sub write_fh {
    my $filename = shift @_;
    my $filehandle;
    open $filehandle, ">$filename" or die $!;
    return $filehandle;
}


########################################################################################"
sub pLine {
my $msg = shift;
print LOG "$msg" x 80 . "\n";
#print ($msg x 20);
}


########################################################################################"

sub sorter {
	$a->[1] cmp $b->[1] ||
     $a->[3] <=> $b->[3]
  || $b->[4] <=> $a->[4]
}

########################################################################################"

sub detectStat {
my ($file, $outfile, $score, $sequence_data_ref)=@_;
  my @terms;
  my %genome=%{$sequence_data_ref};

  my $fh = &read_fh($file);
  my $out =&write_fh($outfile);
  while (<$fh>) { chomp; push @terms, [split /\t/]; }
  my $Ids_ref=extractIds(\@terms);
  foreach my $id (@$Ids_ref) {
  my ($all, %pal, %str);
  my $len = $genome{$id}{len};
  my $palN=0; my $palY=0; my $strP=0; my $strR=0;

	for my $item (@terms) {
    		if ($item->[1] eq $id) {
			if ($item->[9] >= $score) { # 0 score means all the count
			$pal{$item->[10]}++;
			$str{$item->[7]}++;
			}
		$all++;
    		}
	}
  if ($pal{0}) { $palN=$pal{0};}
  if ($pal{1}) { $palY=$pal{1};}
  if ($str{P}) { $strP=$str{P};}
  if ($str{R}) { $strR=$str{R};}
  print $out "$id\t$all\t$palN\t$palY\t$strP\t$strR\n";
}

}

########################################################################################"

sub analyticStat {
my ($file, $outfile, $score, $sequence_data_ref)=@_;
  my @terms;
  my %genome=%{$sequence_data_ref};

  my $fh = &read_fh($file);
  my $out =&write_fh($outfile);
  while (<$fh>) { chomp; push @terms, [split /\t/]; }
  my $Ids_ref=extractIds(\@terms);
  my $feature_ref=extractFeature(\@terms);

  foreach my $id (@$Ids_ref) {
  my %feature;
	for my $item (@terms) {
	if ($item->[11] eq $id) { if ($item->[9] >= $score) { $feature{$item->[13]}++; } } #$item->[9] >= 0 for score filter
	}
  print $out "$id\t";
  my $len = $genome{$id}{len};

  foreach my $f(@$feature_ref) {
	if ($feature{$f}) { print $out "$feature{$f}/$len\t"; }
	else { print $out "0\t";} }
  print $out "\n";
  }

}


########################################################################################"
sub extractIds {
my ($terms_ref)=@_;
my @allIds;
for my $item (@$terms_ref) {
     if ($item->[1]) { push @allIds , $item->[1];}
}
my @allIds_uniq=uniq(@allIds);
return \@allIds_uniq;
}

########################################################################################"
# Checks if a provided two coordinates overlaps or not it return 1 if overlaps
sub checkCorOverlaps {
my ($x1, $x2, $y1, $y2)=@_;
return $x1 <= $y2 && $y1 <= $x2;
}

sub checkCov {
my ($file, $name, $cor1, $cor2)=@_;
my @covData; my $sum=0; my $avCov=0;
open(my $fh, '<:encoding(UTF-8)', $file) or die "Could not open file '$file' $!";
while (my $row = <$fh>) {
  	chomp $row;
	my @data = split /\t/, $row;
	next if $name ne $data[0];
	my $res=checkCorOverlaps ($data[1], $data[2], $cor1, $cor2);
	if ($res) {push @covData, $data[3]}

}
$sum += $_ for @covData;
$avCov=$sum/scalar(@covData);
return $avCov;
}

sub extractSeq_genome {
my ($file, $chr, $st, $ed)=@_;
my $db = Bio::DB::Fasta->new($file);
my $seq = $db->seq($chr, $st => $ed);
return $seq;
}

sub checkPal {
my ($seq, $param)=@_;
my $pp = qr/(?: (\w) (?1) \g{-1} | \w? )/ix;
    while ($seq =~ /(?=($pp))/g) {
        return "$-[0] - $1" if length($1) > $param->{palsize};
    }
}


##Recontruct the breakpoints
sub reconstructTar {
my ($finalEBA, $allHSB, $threshold, $length, $spsFile, $SpsNumber, $refName, $param)=@_;

my $version=0.1;

if (-d "$param->{out_dir}/output_$refName") {
deldir("$param->{out_dir}/output_$refName"); # or deldir($ARGV[0]) to make it commandline
} else { mkdir "$param->{out_dir}/output_$refName"; }

if (-f "Reconstruction_$refName.stats") { unlink "Reconstruction_$refName.stats";}
my $InFile=$finalEBA; #final_classify.eba file
print "$InFile\n";

#Store species
my @SpsArray;
open SPSFILE, "$spsFile" or die $!;
if (-f "Reconstruction_$refName.stats") { unlink "Reconstruction_$refName.stats";}
while (<SPSFILE>) {
	my $SpsLine=$_; chomp $SpsLine; @SpsArray=split /,/, lc($SpsLine);  $SpsNumber = scalar (@SpsArray); } ## It read the species names from sps.txt file ... need to improve !!!
	my $SpsArrayTabed=join("\t", @SpsArray);

close SPSFILE or die "could not close file: $!\n";


foreach my $spsName(@SpsArray) {
my $outFile1="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.tar1";
my $outfile2="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.tar2";
my $outfile3="$param->{out_dir}/output_$refName/$spsName"."_brk_$refName.tar3";

open (OUTFILE1, ">$outFile1") or die "$0: open $outFile1: $!";
open (OUTSTAT, ">>$param->{out_dir}/output_$refName/Reconstruction_$refName.stats") or die "$0: open $param->{out_dir}/output_$refName/Reconstruction_$refName.stats: $!";

	open INFILE,  $InFile or die "$0: open $InFile: $!";
	my @array; my @index; my @nameArray; my $in; my $countReal; my @done; my $countBrk; my $countGap; my $total;
	while (<INFILE>) {
		chomp;
		my $line=trim($_);
		my @tmp = split /\t/, lc ($line);
		if ($. == 1) {
			@nameArray = split /\t/, $line;
			(@index)= grep { $nameArray[$_] eq "$spsName" } 0..$#nameArray;
			$in=$index[0];
			next;
		}
		next if $_ =~ /^\s*#/;
		next if !$tmp[$in];

		my @val = split /\,/, $tmp[$in];

		# In case more than one breakpoints ( separated with comma)
		foreach my $l(@val) {
		my $line=trim($l);

		my $seenIn=isInList ($l, @done);
		if (!$seenIn) { $countReal++; }
		if (index($l, "breakpoints") != -1) { $countBrk++; } elsif (index($l, "gap") != -1) { $countGap++; }
		$total++;

		print OUTFILE1 "$l\t$tmp[$SpsNumber+1]\t$tmp[$SpsNumber+7]\t$tmp[$SpsNumber+9]\t$tmp[$SpsNumber+10]\t$tmp[$SpsNumber+11]\t$tmp[$SpsNumber+12]\t$tmp[$SpsNumber+13]\n";

		push @done, $l;
		#push (@array,$_);
		#72412203--72417432=Breakpoints+0.925	Breakpoints	chicken:6.22711733452034e-07	237.8835341	0.05	20	1	18
		}
	}
	if (-z "Reconstruction_$refName.stats") { print OUTSTAT "spsName\tcountReal\tcountBrk\tcountGap\ttotal\n";}
	print OUTSTAT "$spsName\t$countReal\t$countBrk\t$countGap\t$total\n";
	undef @done;
	close INFILE or die "could not close $InFile file: $!\n";
close OUTSTAT or die "could not close OUTSTAT file: $!\n";
close OUTFILE1 or die "could not close $outFile1 file: $!\n";

reconstructBrk($outFile1, $outfile2, $spsName, $outfile3, $threshold, $allHSB, $param);
undef @index;
}

sub isInList {
	my $needle = shift;
	my @haystack = @_;
	foreach my $hay (@haystack) {
		if ( $needle eq $hay ) {
			return 1;
		}
	}
	return 0;
}


sub deldir {
  my $dirtodel = pop;
  my $sep = '\\'; #change this line to "/" on linux.
  opendir(DIR, $dirtodel);
  my @files = readdir(DIR);
  closedir(DIR);

  @files = grep { !/^\.{1,2}/ } @files;
  @files = map { $_ = "$dirtodel$sep$_"} @files;

  @files = map { (-d $_)?deldir($_):unlink($_) } @files;

  rmdir($dirtodel);
}

}


sub uniq {
my %seen;
return grep { !$seen{$_}++ } @_;
}

# Checks if a provided two coordinates overlaps or not
#my $OverRes = EBALib::CommonSubs::checkCorOverlaps ($val1[0],$val1[1],$val2[0],$val2[1]);
#if ($OverRes) { $j1++;} else {$j2++; }
sub checkOverlaps {
my ($x1, $x2, $y1, $y2)=@_;
return $x1 <= $y2 && $y1 <= $x2;
}

#Store the classification file
sub trim($)
{
my $string = shift;
$string =~ s/^[\t\s]+//;
$string =~ s/[\t\s]+$//;
$string =~ s/[\r\n]+$//; ## remove odd or bad newline ...
return $string;
}

sub reformatTable {
my ( $infile1, $finalOut) = @_;

open(my $frout, ">$finalOut") or die "Could not open file '$finalOut' $!"; #countSTAT_<SPS> file

open(my $infh, '<:encoding(UTF-8)', $infile1) or die "Could not open file '$infile1' $!"; #countSTAT_<SPS> file
while (<$infh>) {
	chomp $_; my $line=$_;
	$line=trim($line);
	next if $line =~ /^\s*$/;
	my @tmpLine = split /\t/, $line;
	for (my $val1=0; $val1<=scalar(@tmpLine); $val1++) {
		if ($val1 == 1) {
		print $frout "$tmpLine[0]\tBreak\t$tmpLine[1]\n";
		}
		elsif ($val1 == 2) {
		print $frout "$tmpLine[0]\tDetected\t$tmpLine[3]\n";
		}
		elsif ($val1 == 3) {
		print $frout "$tmpLine[0]\tReal\t$tmpLine[4]\n";
		}
		else {

		}
	}
    }
close $infh;

}


#subroutines here ------

sub Who {
print "autoConTAMPR by Jitendra, Nico, JF and Karine\n";
}

#To recerse complement
sub reverse_complement_IUPAC {
my $dna = shift;
my $revcomp = reverse($dna);
$revcomp =~ tr/ABCDGHMNRSTUVWXYabcdghmnrstuvwxy/TVGHCDKNYSAABWXRtvghcdknysaabwxr/;
return $revcomp;
}

#Reformat the outfile to plot
sub reformatOut { 
my ($inFile, $outFile)=@_;
my $infh = IO::File->new("$inFile",'r');
my $oufh = IO::File->new("$outFile",'w');
if (-z $outFile) { print $oufh "name\tcnt\tfreq\n";}
    while (<$infh>) {
	chomp;
	my @tmpLine = split '\t', $_;
	my @foo = splice( @tmpLine, 4, 3 );
	my $cnt=0; my $freqName='NA';
		foreach my $tLine (@foo) {
			if ($cnt==0) { $freqName='Second'; } elsif ($cnt == 1) { $freqName='Third'; } else {$freqName='Fourth';}
			print $oufh "$tmpLine[1]:$tmpLine[2]:$tmpLine[3]\t$tLine\t$freqName\n";
			$cnt++;
		}
   }
close $infh;

}

#Return name of both files
#File should be named in this format:H_D14-HCO1_*.scf direction_indivisualName-markerName_*.scf
sub returnFilePath {
my ($dir, $mName, $wfh)=@_;
my @fArray; my $wDir;
    opendir(DIR, $dir) or die  die "Can't open $dir: $!";
    while (my $file = readdir(DIR)) {
        # Use a regular expression to ignore files beginning with a period
        next if ($file =~ m/^\./);
	if (index(lc ($file), lc ($mName)) != -1) {
	push @fArray, $file;
	}
    }
    closedir(DIR);
if (scalar(@fArray) > 2) { print "Something looks wrong, have ($#fArray+1) AND I expect 2 file\n"; exit;} 
elsif (!@fArray) { print $wfh "No chromatogram file found for $mName, Check filename carefully\n";}
return (\@fArray);
}


#Check the intensity of the peak
sub checkIntensity {
my ($conAln, $chromoAln, $chromoFile, $verbose, $logfh) = @_;
#both aln will be of same size
my @conAln = (split //, $conAln); 
my @chromoAln = (split //, $chromoAln);
my $gCnt=0; my $cgCnt=0; my $second=0; my $third=0; my $four=0;
for (my $i=0; $i<=$#chromoAln; $i++){
	next if $conAln[$i] eq $chromoAln[$i]; #if aligned properly
	#Now check only mismatches
	if ($chromoAln[$i] eq "*") {$gCnt++; next;} #Next if '-' in alignment file !!!! why such cases need to resolve 
	if ($conAln[$i] eq "*") {$cgCnt++; next;}
	my $iNEW=$i-$gCnt; #To ignore the '-' count for index
	print $logfh "$conAln[$i] <=> $chromoAln[$i] --- $i => $iNEW <<<\n" if $verbose;
	my ($gb,$ab,$tb,$cb)=intExtractor($chromoFile, $iNEW, $chromoAln[$i]);
	
	my @allFreq=("$gb","$ab","$tb","$cb");
	my %baseHash; my @allScr;
	foreach my $base (@allFreq) { 
		my @tmpBase = split '\:', $base; 
		$baseHash{$tmpBase[0]}=$tmpBase[1];
		push @allScr, $tmpBase[1];
	}
	my @allScore = sort {$b<=>$a} @allScr;
	my( $index ) = grep { $allScore[$_] eq $baseHash{$conAln[$i]}} 0..$#allScore; #Index will be one lower ... begin with 0
	if ($index == 1) {$second++;} elsif ($index == 2) {$third++;} elsif ($index == 3) {$four++;}
	print $logfh "$baseHash{$conAln[$i]}\t$index\n" if $verbose;
	print $logfh "$gb,$ab,$tb,$cb -->>\n" if $verbose;
	undef %baseHash; undef @allScr;
}
print $logfh "$second\t$third\t$four\n" if $verbose;
#print "$second\t$third\t$four\n" if $verbose;
return ($second, $third, $four);
}

#Intensity of the base index
sub intExtractor {
my ($scfFile, $index,$iTop)=@_;
my $scf = Bio::SCF->new("$scfFile");
my $sample_index = $scf->index($index);
my ($g,$a,$t,$c) = map { $scf->sample($_,$sample_index) } qw(G A T C);
#print "G:$g,A:$a,T:$t,C:$c\n";
return ("g:$g","a:$a","t:$t","c:$c");
}

#Read scf file and conver it in fasta
sub extractChroSeq {
my ($inFile, $direction) = @_;
my %scf;
tie %scf, 'Bio::SCF', "$inFile";
my $scf = Bio::SCF->new("$inFile");
my $chroSeq=''; my $chroId='NA';
#my $write_fh = IO::File->new("$outFile",'w');
#$read_fh = IO::File->new("/tmp/msg",'r');
#for bases of each position
for (my $i=0; $i < $scf{bases_length}; $i++){
	#print $write_fh ">$inFile._$scf{bases_length}\n" if $i == 0;
#print $write_fh $scf{bases}[$i];
	$chroId= ">$inFile._$scf{bases_length}" if $i == 0;
	$chroSeq=lc ("$chroSeq"."$scf{bases}[$i]");
}
#print $write_fh "\n";
return ($chroId, $chroSeq);
}


#To find marker in other species
sub findMarker {
my ($worksheet, $speciesName, $marker, $yesCol)=@_;
foreach my $row (1 .. $worksheet->{maxrow}) {
	next if $row == 1; #avoid row contain header
	next if $speciesName eq $worksheet->{"A$row"}; #Not check in same species -- as it is know to be already contaminated
    	foreach my $col (1 .. $worksheet->{maxcol}) {
        	my $cell = cr2cell ($col, $row);
		#next if $noCell ne $cell; # Check only on the same marker cell.
		my @am = split '\s', $marker;
		my @arrayMark = split '\-', $am[0];
		foreach my $mark (@arrayMark) {
			if (index($worksheet->{$cell}, $mark) != -1) {
				next if $yesCol != $col; #remove this line if want to seqrch all cross the table/sheet
   				#print "'$worksheet->{$cell}' contains '$speciesName'\n";
				#return the first one identified -- might need to store all and return				
				return ($worksheet->{$cell}, $col, $row);
			}
		}
	}
}
return ('Not Found',0,0);
}

#find the marker name
sub findName {
my ($worksheet, $mRow, $nCol)=@_;
my @row = Spreadsheet::Read::row($worksheet, $mRow);
if ($row[$nCol]) { return $row[$nCol];} else {return 'NA';}
}


#extract the consensus sequence
sub extractSeq {
my ($mFasta, $mInd, $direction)=@_;
my $mIndLC = lc($mInd);
my $db = Bio::DB::Fasta->new( $mFasta );
#my $write_fh = IO::File->new("tmpConsensus.fa",'w');
my $sequence = $db->seq($mIndLC);
if  (!defined( $sequence )) {
      print "Sequence $mIndLC not found. STRANGE !!!\n Seems your consensus file have some missing sequences\n"; exit;
    }
#print ">$mInd\n", "$sequence\n";
#print $write_fh ">$mId\n", "$sequence\n";
my $fSeq='NA';
if ($direction eq 'left') { $fSeq = reverse_complement_IUPAC($sequence); } else { $fSeq = $sequence; }
return ($mIndLC, $fSeq);
}


# Needleman-Wunsch Algorithm (global alignment)
sub alnNW {
my ($seq1, $seq2, $verbose, $logfh) = @_;
# scoring scheme used for autoConTAMPR
my $MATCH    =  1; # +1 for letters that match
my $MISMATCH = -2; # -1 for letters that mismatch
my $GAP      = -3; # -1 for any gap

# initialization the matrix
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
    $matrix[0][$j]{score}   = $GAP * $j;
    $matrix[0][$j]{pointer} = "left";
}
for (my $i = 1; $i <= length($seq2); $i++) {
    $matrix[$i][0]{score}   = $GAP * $i;
    $matrix[$i][0]{pointer} = "up";
}

# lets fill the matrix
for(my $i = 1; $i <= length($seq2); $i++) {
    for(my $j = 1; $j <= length($seq1); $j++) {
	# dScore / DiagonalScore 
	# lScore / leftScore 
	# uScore / UpScore
        my ($dScore, $lScore, $uScore); 

        # calculate the match score
        my $letter1 = substr($seq1, $j-1, 1);
        my $letter2 = substr($seq2, $i-1, 1);                            
        if ($letter1 eq $letter2) {
            $dScore = $matrix[$i-1][$j-1]{score} + $MATCH;
        }
        else {
            $dScore = $matrix[$i-1][$j-1]{score} + $MISMATCH;
        }

        # calculate the gap scores
        $uScore   = $matrix[$i-1][$j]{score} + $GAP;
        $lScore = $matrix[$i][$j-1]{score} + $GAP;

        # choose the best score for alignment
        if ($dScore >= $uScore) {
            if ($dScore >= $lScore) {
                $matrix[$i][$j]{score}   = $dScore;
                $matrix[$i][$j]{pointer} = "diagonal";
            }
        else {
                $matrix[$i][$j]{score}   = $lScore;
                $matrix[$i][$j]{pointer} = "left";
            }
        } else {
            if ($uScore >= $lScore) {
                $matrix[$i][$j]{score}   = $uScore;
                $matrix[$i][$j]{pointer} = "up";
            }
            else {
                $matrix[$i][$j]{score}   = $lScore;
                $matrix[$i][$j]{pointer} = "left";
            }
        }
    }
}

# trace-back the matrix 

my $align1 = "";
my $align2 = "";

# start at last cell of matrix
my $j = length($seq1);
my $i = length($seq2);

while (1) {
    last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix

    if ($matrix[$i][$j]{pointer} eq "diagonal") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= substr($seq2, $i-1, 1);
        $i--;
        $j--;
    }
    elsif ($matrix[$i][$j]{pointer} eq "left") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= "*";
        $j--;
    }
    elsif ($matrix[$i][$j]{pointer} eq "up") {
        $align1 .= "*";
        $align2 .= substr($seq2, $i-1, 1);
        $i--;
    }  
}

$align1 = reverse $align1;
$align2 = reverse $align2;
print $logfh "$align1\n" if $verbose; #Consensus
print $logfh "$align2\n" if $verbose; #chromoSeq

return ($align1, $align2);
}



#ContTAMPR on specific species 
sub splConTAMPR {
my ($table, $consensus, $markerDir, $outfile, $plot, $brutal, $verbose, $logfile)=@_;
my $logfh = IO::File->new("$logfile",'w') if $verbose;
my $ofh = IO::File->new("$outfile",'w') if $verbose;
my $seperator='*' x 100;

my $ss = ReadData ("$table", attr => 1)->[1]; #Reading sheet 1 
my $progress = Term::ProgressBar->new($ss->{maxrow});
foreach my $row (1 .. $ss->{maxrow}) {
	next if $row == 1; #avoid row
	my $speciesName;
	print $logfh "$seperator\nChecking on row number $row\n$seperator\n" if $verbose;
    	foreach my $col (1 .. $ss->{maxcol}) {
        	my $cell = cr2cell ($col, $row);
        	#printf "%s %-3s %d  ", $cell, $ss->{$cell}, $ss->{attr}[$col][$row]{merged};
		if ($col == 1) { $speciesName=$ss->{$cell}; next;}		
		#print "$cell\t";
		if (index($ss->{$cell}, "($speciesName)") != -1) {
   			print $logfh "'$ss->{$cell}' contains '$speciesName'\n" if $verbose;
		}
		else {
			next if $col == 2; #ignore as this is indivisual name
			next if !$ss->{$cell}; #Ignore if empty	
			print $logfh "'$ss->{$cell}' NOT contains species '$speciesName\t$cell'\n" if $verbose;
			my ($markerName, $mCol, $mRow)=findMarker($ss, $speciesName, $ss->{$cell}, $col);
			print $logfh "$mCol,$mRow\t$ss->{$cell}\t$markerName --->\n" if $verbose;
			my $indCol=1; #column with ind name !!!!!!!!!
			my ($indRName)=findName($ss, $row, $indCol);
			my ($indCName)=findName($ss, $mRow, $indCol);
			print $logfh "$indRName\t$indCName\n" if $verbose;
			my ($fileArray_ref) = returnFilePath($markerDir,$indRName);
			my @fArray= @$fileArray_ref;
			my $lHand=''; my $hHand=''; my $direction='NA';

			#Lets work on both scf files
			my $tSec=1; my $tThi=1; my $tFou=1;
			foreach my $f (sort {$a cmp $b} @fArray) {
				if ($f =~ m/^L/)  { print $logfh "$f\n" if $verbose; $direction='left';}
				if ($f =~ m/^H/) {print $logfh "$f\n" if $verbose; $direction='right';}
				print $logfh "$f\t$markerDir/$f\n" if $verbose;
				my $cSeq = extractSeq($consensus,$indCName, $direction);
				my ($chroId, $chroSeq)=extractChroSeq ("$markerDir/$f", "$direction");
				print $logfh "Chromatogram:$chroId\n$chroSeq\n\n" if $verbose; 
				print $logfh "Consensus:$indCName\n$cSeq\n\nAlignment\n" if $verbose; 
				my ($cAln, $chroAln) = alnNW ($cSeq, $chroSeq, $verbose, $logfh);
				my ($sec, $thi, $fou) = checkIntensity ($cAln, $chroAln, "$markerDir/$f", $verbose, $logfh);
				$tSec=$tSec+$sec; $tThi=$tThi+$thi; $tFou=$tFou+$fou; 
			}
			my @allInt = ($tSec, $tThi, $tFou);
			#check the randomness
			my $normal= (sum( @allInt)/3);
			my $chi = chisquare( @allInt);
			my $chiVal = chi_quared( observed => [ @allInt ], expected => [ ($normal) x 3 ]); 
			#And that prints out 0.018360, the probability of that set of freq occurring by chance.
			my $chisprob=Statistics::Distributions::fprob(2,3,$chi);
			my $chichi='NA';
			$tSec= $tSec-1; $tThi= $tThi-1; $tFou= $tFou-1; # Just to count 1 less (because cnt start from 1)
			print $ofh "$speciesName\t$indRName\t$indCName\t$ss->{$cell}\t$tSec\t$tThi\t$tFou\t$chiVal\t$chi\t$chisprob\t$chichi\n";

			if ($brutal eq "yes") { print $logfh "\nExiting with ONE run --- \n\n" if $verbose; exit; }
		}
	}
	$progress->update($row)
     	#print "\n";
}
close $logfh;
close $ofh;

reformatOut($outfile, 'plotData');
my $workDir = getcwd;
if ($plot) { system "Rscript $workDir/utils/conPlot.R plotData";}


}

#Plot the contamination graph
sub circosPlot {
my ($ConTAMPR_out, $table, $param, $log)=@_;

open(OUTseq, ">$param->{out_dir}/circos/circos.sequences.txt") or die("Cannot create circos.sequences.txt\n");   
open(OUTsegdup, ">$param->{out_dir}/circos/circos.segdup.txt") or die("Cannot create circos.segdup.txt\n");     

#Store into hash
my $hash_ref = storeIntoHash ($table, 3);
my %hash = %$hash_ref;

my $ss = ReadData ("$table", attr => 1)->[1]; #Reading sheet 1 
my $book = Spreadsheet::Read->new ("$table");
my $sheet = $book->sheet (1);     # OO
my @col = $sheet->column(2); #foreach (@col) { print "$_\n";} #exit; # Column 2 is indivisual name
shift @col;

my $seqCnt; my %seqHash; 
foreach my $cVal(@col) {
	$seqCnt++;
	my $seq="seq"."$seqCnt";
	my $markerLen=$param->{len}; #Legth of the expected size of the markers
	my $tmpRow=$seqCnt+1;# Depends upon unsorted array ... dont sort
	my $cell = cr2cell (1,$tmpRow); 
        my $spsName=trim ($ss->{$cell});
	my $blkColor='NA';
	$blkColor=$hash{$spsName};

	#chr - seq1132 contig_1131 0 18652167 red_a4
	my $seqLine="chr - $seq $cVal 0 $markerLen $blkColor";	
		print OUTseq "$seqLine\n";
		$seqHash{$cVal}=$seq;
	}
close OUTseq;

#print "Generating circos segdup file\n";

foreach my $cVal(@col) {
	open(INConTAMPR, $ConTAMPR_out) or die("Cannot open $ConTAMPR_out\n -Terminated !! \n Either you forgot to run --validation step OR Folder name is changed\n");
	my $markerLen=$param->{len};
	while (<INConTAMPR>) {
		next if /^\s*$/;
		my @splitline = split(/\t/, $_);
		#my $blkCnt++;
		$splitline[1] =~ s/\_//g; #just for nex _A24_ underscore correction
		next if $splitline[7] > $param->{threshold}; # This threshold to ignore non-interesting links in circos plot (0.01)
		if ($splitline[1] eq $cVal) {
			my $blkName="block_000"."$.";
			#block_0000 seq1132 1 18652167 color=chr13_a2
			#print OUTsegdup "$splitline[1] eq $cVal\n";
			my $clink=$seqHash{$splitline[2]};
			my $seq=$seqHash{$cVal};
			my $num=int(rand(20));
			my $clr="chr"."$num"."_a2";
			next if $splitline[7] >= 0.05;
			print OUTsegdup "$blkName $seq 1 $markerLen color=$clr\n";
			print OUTsegdup "$blkName $clink 1 $markerLen color=$clr\n";
			#print OUTcrules "</rule>\n";
		}
	} 
close INConTAMPR;	
}
close OUTsegdup;

system ("$param->{circos}/circos -silent -conf $param->{out_dir}/circos/circos.conf");
}


## Subroutine to open & read logfile 
sub extractSpsLine {
my ($infile, $outfile, $sps)=@_;  
    open (IN, $infile); 
    open (OUT, ">", $outfile) or die "$!"; 
    while (<IN>) {
	chomp $_;
	my @tmpLine = split '\t', $_;
        if ($tmpLine[0] eq $sps) {  
            print OUT "$_\n"; 
        } 
    } 
    close OUT;
} 


# Smith-Waterman  Algorithm LOCAL
sub alnSW {
my ($seq1, $seq2, $verbose, $logfh, $scfLen, $chroId, $indCName) = @_;

# scoring scheme
my $MATCH    =  1; # +1 for letters that match
my $MISMATCH = -1; # -1 for letters that mismatch
my $GAP      = -10; # -1 for any gap

# initialization
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
    $matrix[0][$j]{score}   = 0;
    $matrix[0][$j]{pointer} = "none";
}
for (my $i = 1; $i <= length($seq2); $i++) {
    $matrix[$i][0]{score}   = 0;
    $matrix[$i][0]{pointer} = "none";
}

# fill
my $max_i     = 0;
my $max_j     = 0;
my $max_score = 0;

for(my $i = 1; $i <= length($seq2); $i++) {
    for(my $j = 1; $j <= length($seq1); $j++) {
        my ($diagonal_score, $left_score, $up_score);
        
        # calculate match score
        my $letter1 = substr($seq1, $j-1, 1);
        my $letter2 = substr($seq2, $i-1, 1);       
        if ($letter1 eq $letter2) {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
        }
        else {
            $diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
        }
        
        # calculate gap scores
        $up_score   = $matrix[$i-1][$j]{score} + $GAP;
        $left_score = $matrix[$i][$j-1]{score} + $GAP;
        
        if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
            $matrix[$i][$j]{score}   = 0;
            $matrix[$i][$j]{pointer} = "none";
            next; # terminate this iteration of the loop
        }
        
        # choose best score
        if ($diagonal_score >= $up_score) {
            if ($diagonal_score >= $left_score) {
                $matrix[$i][$j]{score}   = $diagonal_score;
                $matrix[$i][$j]{pointer} = "diagonal";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        } else {
            if ($up_score >= $left_score) {
                $matrix[$i][$j]{score}   = $up_score;
                $matrix[$i][$j]{pointer} = "up";
            }
            else {
                $matrix[$i][$j]{score}   = $left_score;
                $matrix[$i][$j]{pointer} = "left";
            }
        }
        
        # set maximum score
        if ($matrix[$i][$j]{score} > $max_score) {
            $max_i     = $i;
            $max_j     = $j;
            $max_score = $matrix[$i][$j]{score};
        }
    }
}

# trace-back

my $align1 = "";
my $align2 = "";

my $j = $max_j;
my $i = $max_i;

my $starti=0; my $startj=0; 
while (1) {
    last if $matrix[$i][$j]{pointer} eq "none";
    
    if ($matrix[$i][$j]{pointer} eq "diagonal") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= substr($seq2, $i-1, 1);
        $i--; $j--;
$starti=$i; $startj=$j;
    }
    elsif ($matrix[$i][$j]{pointer} eq "left") {
        $align1 .= substr($seq1, $j-1, 1);
        $align2 .= "*";
        $j--;
$starti=$i; $startj=$j;
    }
    elsif ($matrix[$i][$j]{pointer} eq "up") {
        $align1 .= "*";
        $align2 .= substr($seq2, $i-1, 1);
        $i--;
$starti=$i; $startj=$j;
    }   
}

#print "$starti\t$startj\n";
$align1 = reverse $align1;
$align2 = reverse $align2;
print $logfh "$align1\n" if $verbose;
print $logfh "$align2\n" if $verbose;

my $alnLen=length ($align1); 
#print "$scfLen $alnLen\n";
if ($scfLen < ($alnLen/2)) { print "The length of alignment is too small between $chroId (length: $scfLen), $indCName (alignment length:$alnLen) !! \nPossible orientation error. Flip the scf file strand and try again\n"; exit(1);}

return ($align1, $align2, $starti, $startj);
}


#Check the intensity of the peak
sub checkIntensityLocal {
my ($conAln, $chromoAln, $chromoFile, $sti, $stj, $verbose, $logfh, $cName) = @_;
#both aln will be of same size
my @conAln = (split //, $conAln); 
my @chromoAln = (split //, $chromoAln);
my $gCnt=0; my $cgCnt=0; my $second=0; my $third=0; my $four=0;
print $logfh "\nReal start for i(scf) => $sti, for j(consensus) => $stj\n\n" if $verbose;
print $logfh "\n#/MISMATCH location and INTENSITY\n\n" if $verbose;
#start from nex index $sti
for (my $i=0; $i<=$#chromoAln; $i++){
	next if $conAln[$i] eq $chromoAln[$i]; #if aligned properly
	#Now check only mismatches
	if ($chromoAln[$i] eq "*") {$gCnt++; next;} #Next if '*' in alignment file !!!! why such cases need to resolve 
	if ($conAln[$i] eq "*") {$cgCnt++; next;}
	my $iNEW=(($i+$sti)-$gCnt); #real index begin at ... ignore * $gCnt
	print $logfh "$conAln[$i] <=> $chromoAln[$i] --- $i => $iNEW <<<\n" if $verbose;
	my ($gb,$ab,$tb,$cb)=intExtractor($chromoFile, $iNEW, $chromoAln[$i]);
	
	my @allFreq=("$gb","$ab","$tb","$cb");
	my %baseHash; my @allScr;
	foreach my $base (@allFreq) { 
		my @tmpBase = split '\:', $base; 
		$baseHash{$tmpBase[0]}=$tmpBase[1];
		push @allScr, $tmpBase[1];
	}
	my @allScore = sort {$b<=>$a} @allScr;
	my( $index ) = grep { $allScore[$_] eq $baseHash{$conAln[$i]}} 0..$#allScore; #Index will be one lower ... begin with 0
	if ($index == 1) {$second++;} elsif ($index == 2) {$third++;} elsif ($index == 3) {$four++;}
	print $logfh "$baseHash{$conAln[$i]}\t$index\n" if $verbose;
	print $logfh "$gb,$ab,$tb,$cb -->>\n" if $verbose;
	undef %baseHash; undef @allScr;
}
print $logfh "\nFinal count:\n second=> $second\n third=> $third\n fourth=> $four\n" if $verbose;
print $logfh "\nGAP in scf $gCnt, GAP in consensus $cgCnt : Between $chromoFile, $cName\n\n" if $verbose;
print $logfh "\nEND/#\n\n" if $verbose;
#print "$second\t$third\t$four\n" if $verbose;
return ($second, $third, $four);
}


sub validateHelp {
  my $ver = $_[0];
  print "\n  autoConTAMPR --validate $ver \n\n";
  print "    Usage: autoConTAMPR.pl --validate/-v --conf/-c <configuration file>\n\n";
  print "    To check the contamination in the consensus markers fasta sequences\n\n";
  print "    The path to a valid autoConTAMPR configuration file. This file contains all parameters needed to execute  autoConTAMPR.\n\n";
print "Point to NOTE:
1.Consensus fasta file should be in lower case : tr A-Z a-z < input.fa
2.Name of the all the markers file should begin with L or H
3.The constamination table have special format (see the default table in sample data)\n";
exit(1);
}

sub installHelp {
  my $ver = $_[0];
  print "\n  autoConTAMPR --full $ver \n\n";
  print "    Usage: autoConTAMPR.pl --install/-i --conf/-c <configuration file>\n\n";
  print "    To run complete in one GO \n\n";
  print "    The path to a valid autoConTAMPR configuration file. This file contains all parameters needed to execute  autoConTAMPR.\n\n";

exit(1);
}

sub plotHelp {
  my $ver = $_[0];
  print "\n  autoConTAMPR --plot $ver \n\n";
  print "    Usage: autoConTAMPR.pl --plot/-p --conf/-c <configuration file>\n\n";
  print "    To plot the autoConTAMPR results \n\n";
  print "    The path to a valid autoConTAMPR configuration file. This file contains all parameters needed to execute  autoConTAMPR.\n\n";

exit(1);
}



#Chi-squared test calculation
sub chi_quared {
    my %arg_for = @_;
    my $observed = $arg_for{observed} 
      // croak q(Argument "observed" required);
    my $expected = $arg_for{expected}
      // croak q(Argument "expected" required);

    @$observed == @$expected or croak q(Input arrays must have same length);

    my $chiSquared = sum map { 
        ( $observed->[$_] - $expected->[$_] )**2 
        / 
        $expected->[$_]
    } 0 .. $#$observed;

    my $dof_freedom = @$observed - 1;
    my $prob        = chisqrprob(
        $dof_freedom,
        $chiSquared
    );

    return $prob;
}
1;

__END__
