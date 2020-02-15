# autoConTAMPR v0.1
```
   --.---. .-. .---
	 / \
	/ ^ \
  --autoConTAMPR v0.1--
```
# autoConTAMPR
autoConTAMPR: Automatic Contingency Table Analysis of Minority Peak Ranks

Implemented by ~~jitendra.narayan@unamur.be~~, jnarayan81@gmail.com, jflot@ulb.ac.be, nicolas.debortoli@unamur.be and karine.vandoninck@unamur.be

## LICENSE

This file is part of autoConTAMPR.

autoConTAMPR is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

autoConTAMPR is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with autoConTAMPR, in a file called COPYING. If not, see <http://www.gnu.org/licenses/>.

## INSTALLATION

See INSTALLATION. This is a standard setuptools setup, if you already know the procedure.

USING PERL 5.x

autoConTAMPR does not support Perl 5.x, and no plans exist to provide such support. For a strong biological analysis package for perl 5, see perl https://www.perl.org/ and bioperl: http://bioperl.org/

Perl modules

You also need to have Perl, BioPerl, and some other modules installed in your
machine:

```
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

```
You can check if you have them installed in your machine with 
'perl -M<module> -e 1'. It will return an error message if it isn't installed.

E.g. "perl -MBio::SeqIO -e 1"

To install these modules, you can do through the CPAN or manually downloading
(http://search.cpan.org/) and compiling them. To use CPAN, you can do by 
writing:

> perl -MCPAN -e 'install <module>'

To install manually, search for the most recent version of these modules and,
for each, download and type the following (should work most of the time, except
for modules like BioPerl, of course):

> tar -zxvf <module.tar.gz>
> perl Makefile.PL
> make
> make test
> make install

For more detail, you could visit: http://bioinformaticsonline.com/blog/view/710/how-to-install-perl-modules-manually-using-cpan-command-and-other-quick-ways

## DOCUMENTATION

For documentation, see autoConTAMPR_manual_v0.1.pdf in the base project folder.

## Quick run

perl autoConTAMPR.pl -t table.csv -m newMarkers -o outTable -p -c COI_full_LC.fas --brutal no --verbose

This documentation may be out of date depending on whether or not the developers did their job and re-generated the documentation before the release. If you suspect that the documentation is out of date, or if you are using code from the repository (and not from a release), you can re-generate the documentation or contact the authors.

## RELEASE HISTORY

0.1.0 - 1 March 2018

## OUTPUT FORMAT

 autoConTAMPR outfile columns:
 
    * SPS         Name of the species
    * NAME	  Name of the marker
    * CON         Name of the consesus marker
    * SPSNAME     Name of the species it is contaminated with
    * SEC         Count of the second best peaks
    * THI         Count of the third best peaks
    * FOR         Count of the fourth best peaks
    * CHI         Chi square values
    * DETAIL      Detail for undersatanding the chi value

 autoConTAMPR logfile :

    * It contain the detail alignment log file for all markers
 
 ConTAMPR plotData column
 
    * NAME      Name of the markers
    * CNT       Total count of the peak
    * FREQ      Peak frequency name
    
 
## ANNOTATION DESCRIPTIONS
Coming Soon

## EXTRA
 Coming Soon

## FAQ

Who can I report bugs to or ask questions?
Please report your issues to ticketing system.

## CONTRIBUTION

Feel free to clone this repository and use it under the licensing terms.

Additionally, as the project is on github, you may submit patches, ticket requests, edit the wiki, send pull requests - anything you like and have the permissions to do. I will enjoy any and all contributions, I'm sure! :)

As always, you can contact the authors at jitendra.narayan@unamur.be, jflot@ulb.ac.be, nicolas.debortoli@unamur.be and karine.vandoninck@unamur.be.

## CITATION

Reply to Cross-Contamination Explains "Inter and Intraspecific Horizontal Genetic Transfers" between Asexual Bdelloid Rotifers (Wilson, Nowell & Barraclough 2018) Jean-François Flot, Nicolas Debortoli, Bernard Hallet, Jitendra Narayan, Karine Van Doninck doi: https://doi.org/10.1101/368209

Cross-Contamination Explains “Inter and Intraspecific Horizontal Genetic Transfers” between Asexual Bdelloid Rotifers , Christopher G. Wilson Reuben W. Nowell Timothy G. Barraclough https://www.cell.com/current-biology/abstract/S0960-9822(18)30706-1
