#!/usr/bin/perl

use strict;
use warnings;
use Spreadsheet::Read;
use autodie;
use File::Temp qw(tempfile);
use File::Copy;
use Cwd qw();
use Carp; 
use IO::File;
use Bio::SCF;
use Bio::DB::Fasta;
use Getopt::Long;
use Cwd;
use File::Remove 'remove';
use File::Path qw(make_path remove_tree);
use Fcntl qw( LOCK_EX LOCK_NB );
use Cwd 'abs_path';
use Statistics::ChiSquare;
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;
use Term::ProgressBar 2.00;
use Bio::Perl;
use 5.012;

print  <<'WELCOME';  
	  /\	
	 /^^\ 
      autoConTAMPR v0.1 [March 1, 2018]  - [+]
Citation - automated ConTAMPR
License: Creative Commons Licence
Bug-reports and requests to: jitendra.narayanATunamur.be, nicolas.debortoliATunamur.be and karine.vandoninck@unamur.be
-------------------------------------------------------------------
USAGE: perl perl autoConTAMPR.pl -t table.csv -m markers -o outTable -p -c COI_full_LC.fas --brutal no --verbose

WELCOME

my (
	$table,
	$consensus,
	$markerDir,
	#$infile,
	$outfile,
	$plot,
	$brutal,
	$logfile,
);

# Default option setting for palindromer tool
my $VERSION='v0.1';
my $verbose=0; 		# Verbose set to 0;
my %options = ();
$brutal='no';
$logfile='autoConTAMPR.log';

GetOptions(
	\%options,
    	'table|t=s'    	=> \$table,        	## Table
    	'consensus|c=s' => \$consensus,        	## consensus marker file
    	#'infile|f=s'    => \$infile,        	## Infile
    	'markerDir|m=s' => \$markerDir,        	## Infile
    	'outfile|o=s' 	=> \$outfile,           ## Outfile
	'plot|p' 	=> \$plot, 		## plot dotplot with R "yes or no"
	'brutal|b=s'	=> \$brutal,
    	'help|?|h!'     => sub { printUsage($VERSION) },
   	'who|w!'     	=> sub { Who($VERSION) },
	'verbose' 	=> \$verbose,
    	'logfile=s' 	=> \$logfile,		## logfile
	
) or die printUsage($VERSION);

if ((!$table) or (!$consensus) or (!$outfile) or (!$markerDir) ) { 
print "ERROR: You might forgot to provide right flags\n";
printUsage($VERSION); exit; }

my $logfh = IO::File->new("$logfile",'w') if $verbose;
my $ofh = IO::File->new("$outfile",'w') if $verbose;
my $seperator='*' x 100;

#Set up the path of all required tools - stored in the same folder
#$ENV{'PATH'} = "/bin:/usr/bin:/usr/bin/env:$ENV{PWD}/lastz-distrib-1.03.73/src:$ENV{PWD}/trf";

#my $workbook  =  ReadData ("$table", sep => "\t");
#my $cell  = $workbook->[1]{B3};         # content of field A3 of sheet 1
#my @row = Spreadsheet::Read::row ($workbook->[1], 3);
#print "$cell\n";
#foreach (@row) { print "$_\n";}

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
			my $workDir=getcwd(); 
			my @fArray= @$fileArray_ref;
			my $lHand=''; my $hHand=''; my $direction='NA';

			#Lets work on both scf files
			my $tSec=1; my $tThi=1; my $tFou=1;
			foreach my $f (sort {$a cmp $b} @fArray) {
				if ($f =~ m/^L/)  { print $logfh "$f\n" if $verbose; $direction='left';}
				if ($f =~ m/^H/) {print $logfh "$f\n" if $verbose; $direction='right';}
				print $logfh "$f\t$workDir\n" if $verbose;
				my $cSeq = extractSeq($consensus,$indCName, $direction);
				my ($chroId, $chroSeq)=extractChroSeq ("$workDir/$markerDir/$f", "$direction");
				print $logfh "Chromatogram:$chroId\n$chroSeq\n\n" if $verbose; 
				print $logfh "Consensus:$indCName\n$cSeq\n\nAlignment\n" if $verbose; 
				my ($cAln, $chroAln) = alnNW ($cSeq, $chroSeq);
				my ($sec, $thi, $fou) = checkIntensity ($cAln, $chroAln, "$workDir/$markerDir/$f");
				$tSec=$tSec+$sec; $tThi=$tThi+$thi; $tFou=$tFou+$fou; 
			}
			my @allInt = ($tSec, $tThi, $tFou);
			#check the randomness
			my $normal= (sum( @allInt)/3);
			my $chi = chisquare( @allInt);
			my $chiVal = chi_quared( observed => [ @allInt ], expected => [ ($normal) x 3 ]); 
			#And that prints out 0.018360, the probability of that set of freq occurring by chance.
			my $chisprob=Statistics::Distributions::fprob(3,3,$chi);
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

if ($plot) { system "Rscript conPlot.R plotData";}


#subroutines here ------
#Help section
sub printUsage {
  my $ver = $_[0];
  print "\n autoConTAMPR $ver\n\n";

  print "Usage: $0 perl autoConTAMPR.pl -t table.csv -m markers -o outTable -p -c COI_full_LC.fas --brutal no --verbose \n\n";
  print	"Options:\n";
  print "	--table|-t	table with the contamination info (specially formated)\n";
  print "	--consensus|--c	consensus markers fasta sequence\n";
  print "	--infile|--f	infile\n";
  print "	--markerDir|-m	marker files in a directory (special name format)\n";
  print "	--outfile|-o	outfile for results\n";
  print "	--plot|-p	plot the result\n";
  print "	--brutal|-b	brutally delete all intermediate file\n";
  print "	--help|-h|?	print the help\n";
  print "	--who|-w	print the developer name\n";
  print "	--verbose	print the log file\n";
  print "	--logfile	name of the logfile\n\n";

print "Point to NOTE:\n";
print "Consensus fasta file should be in lower case : tr A-Z a-z < input.fa\n";
print "Name of the all the markers file should begin with L or H\n";
print "The constamination table have special format (see the default table in sample data)\n\n";

exit;
}

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
			if ($cnt==0) { $freqName='Second'; } elsif ($cnt == 1) { $freqName='third'; } else {$freqName='Fourth';}
			print $oufh "$tmpLine[1]:$tmpLine[2]:$tmpLine[3]\t$tLine\t$freqName\n";
			$cnt++;
		}
   }
close $infh;

}

#Return name of both files
#File should be named in this format:H_D14-HCO1_*.scf direction_indivisualName-markerName_*.scf
sub returnFilePath {
my ($dir, $mName)=@_;
my @fArray; my $wDir;
    opendir(DIR, $dir) or die  die "Can't open $dir: $!";
    while (my $file = readdir(DIR)) {
        # Use a regular expression to ignore files beginning with a period
        next if ($file =~ m/^\./);
	if (index($file, $mName) != -1) {
	push @fArray, $file;
	}
    }
    closedir(DIR);
if (scalar(@fArray) > 2) { print "Something looks wrong, have ($#fArray+1) AND I expect 2 file\n"; exit;} 
elsif (!@fArray) { print "No chromatogram file found for $mName, Check filename carefully\n";}
return (\@fArray);
}


#Check the intensity of the peak
sub checkIntensity {
my ($conAln, $chromoAln, $chromoFile) = @_;
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
      print "Sequence $mIndLC not found. STRANGE !!!\n"; exit;
    }
#print ">$mInd\n", "$sequence\n";
#print $write_fh ">$mId\n", "$sequence\n";
my $fSeq='NA';
if ($direction eq 'left') { $fSeq = reverse_complement_IUPAC($sequence); } else { $fSeq = $sequence; }
return ($mIndLC, $fSeq);
}


# Needleman-Wunsch Algorithm (global alignment)
sub alnNW {
my ($seq1, $seq2) = @_;
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
