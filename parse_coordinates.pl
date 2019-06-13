#!/usr/bin/env perl

####################################################################################################
#                                                                                                  #
# parse_coordinates.pl                                                                             #
#                                                                                                  #
# Author: Katharina Hoff                                                                           # 
#                                                                                                  #
# Contact: katharina.hoff@uni-greifswald.de                                                        #
#                                                                                                  #
# This script is under the GNU General Public License v3.0                                         #
####################################################################################################

use strict;
use warnings;
use Getopt::Long;
use File::Basename;


my $usage = <<'ENDUSAGE';

parse_coordinates.pl     Parse coordinates from custom XML files of electron microscopy images

SYNOPSIS

parse_coordinates.pl --in_stem=stem_name --out=matrix.txt

INPUT FILE OPTIONS

--in_stem=stem_name                 Name stem of input files for a single 3D object,
                                    e.g. amaP_immunogold. Script assumes that name_stem
                                    is followed by increasing integers starting by 1,
                                    script will stop searching for further files if there
                                    is a gap between the integers.
--cfg=file                          Configuration file that contains translation from "Contour names"
                                    to cell structures (e.g. "red" may be "nucleus"). By default,
                                    the provided file translate.cfg in the same directory as 
                                    parse_coordinates.pl is used. The file must be tabulator
                                    separated and contains always contour label (possibly as
                                    pattern) in the first column, and final translation in the second
                                    column.
--out=matrix.txt                    Matrix with 3D coordinates and labels for import into R
--version                           Print version number of script
--help                              Print this help message

DESCRIPTION

Examples:

parse_coordinates.pl --in_stem=amaP_immunogold.--out=matrix.txt

ENDUSAGE

my $version = "1.0";
my $vers;
my $help;
my $in_file_stem;
my $out_file;
my $cfg;

if ( @ARGV == 0 ) {
    print "$usage\n";
    exit(0);
}

GetOptions(
    'in_stem=s'     => \$in_file_stem,
    'cfg=s'         => \$cfg,
    'out=s'         => \$out_file,
    'version!'      => \$vers,
    'help!'         => \$help);

if($help){
	print $usage;
	exit(0);
}

if($vers){
	print "version ".$version."\n";
	exit(0);
}

if(not(defined($in_file_stem))){
	print "ERROR: in file " . __FILE__ ." at line "
	      . __LINE__ ."\nNo input provided!\n";
    print $usage;
    exit(1);
}

if(not(defined($cfg))){
	my $dirname = dirname(__FILE__);
	$cfg = $dirname."/translate.cfg";
}

my $out_exists = 0;

if(not(defined($out_file))){
	print "ERROR: in file " . __FILE__ ." at line "
	      . __LINE__ ."\nNo output file name provided!\n";
    print $usage;
    exit(1);
}elsif(-e $out_file){
	print "WARNING: output file $out_file already exists. Will append results to the already existing file...\n";
	$out_exists = 1;
}

# determine cell name 
my @t = split(/\//, $in_file_stem);
my $nt = @t;
my $cell = $t[$nt-2]."_".$t[$nt-1];

# read Contour translation table
my %transl;
open(CFG, "<", $cfg) or die ("ERROR: in file " . __FILE__ ." at line "
                             . __LINE__ ."\nFailed to open file $cfg for reading!\n");
while(<CFG>){
	chomp;
	my @t = split(/\t/);
	$transl{$t[0]} = $t[1];
}
close(CFG) or die ("ERROR: in file " . __FILE__ ." at line "
                             . __LINE__ ."\nFailed to close file $cfg!\n");

# initialize variables for parsing XML files
my $c = 0;
my %data; # stores coordinate data
my $z;
my $screen = 0;
my $type; # membrane, gold, etc.
my $gold_id = 0; # to summarize single gold particles
my $contour_label;


my @x_coef;
my @y_coef;

# loop over input files
while($c < 200 && not( -e $in_file_stem.$c)){
	$c = $c + 1;
}
while(-e $in_file_stem.$c){
	print "Scanning file $c: " . $in_file_stem.$c ."\n";
	open(INFILE, "<", $in_file_stem.$c) or die ("ERROR: in file " . __FILE__ ." at line "
                                     . __LINE__ ."\nFailed to open file $in_file_stem".$c." for reading!\n");
	# screen XML format of single file for important information

	while(<INFILE>){
		if(m/<Section index="\d+" thickness="(\d+\.\d+)"/){
			$z = $1*$c;
		}
		if(m/xcoef=" (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?)"/){
			@x_coef = ($1, $4, $7, $10, $13, $16);
			# convert float, somehow perl otherwise thinks they are strings?!
			for(my $i=0; $i<6; $i++){
				if($x_coef[$i] =~ m/(-?\d+\.\d+)e(-?\d+)/){
					$x_coef[$i] = $1 * 10 ** $2;
				}
			}
			#foreach(@x_coef){print $_."\n";}
		}elsif(m/ycoef=" (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?) (-?\d+(\.\d+e?-?(\d+)?)?)"/){
			@y_coef = ($1, $4, $7, $10, $13, $16);
			for(my $i=0; $i<6; $i++){
				if($y_coef[$i] =~ m/(-?\d+\.\d+)e(-?\d+)/){
					$y_coef[$i] = $1 * 10 ** $2;
				}
			}
			#foreach(@y_coef){print $_."\n";}

		}
		if(m/<Contour name="([^""]+)/){
			$contour_label = $1;
			my $type_set = 0;
			while (my ($key, $value) = each (%transl)){
				if($contour_label =~ m/$key/){
					$type = $value;
					$type_set = 1;
					$screen = 1;
				}
			}
			if($type_set == 0){
				print("Found and ignored contour label $contour_label\n");
				$screen = 0;
			}
		}
		if(m/\/>/){
			$screen = 0;
		}
		if( ($screen == 1) && (m/(-?\d+\.\d+)\s+(-?\d+\.\d+),/)){
			my $x = $1;
			my $y = $2;
			my $a0 = $x_coef[0];
			my $a1 = $x_coef[1];
			my $a2 = $x_coef[2];
			my $b0 = $y_coef[0];
			my $b1 = $y_coef[1];
			my $b2 = $y_coef[2];
			my $trans_x =  (-$a0 * $b2 + $a2 * $b0 - $a2 * $y + $b2 * $x) / ($a1 * $b2 - $a2 * $b1);
			my $trans_y = ($b1 * ($a0 - $x) + $a1 * ($y - $b0)) / ($a1 * $b2 - $a2 * $b1);
			if($type =~ m /gold/){
				push(@{$data{"gold"}{$gold_id}{"x"}}, $trans_x);
				push(@{$data{"gold"}{$gold_id}{"y"}}, $trans_y);
				$data{"gold"}{$gold_id}{"z"} = $z;
				$data{"gold"}{$gold_id}{"cell"} = $cell;
			}else{
				push(@{$data{$type}}, $trans_x . "\t". $trans_y ."\t$z\t$type\t$cell\n");

			}
		}
	}
	close(INFILE ) or die ("ERROR: in file " . __FILE__ ." at line "
                                     . __LINE__ ."\nFailed to close file $in_file_stem".$c."!\n");
	$c++;
}

# write important content in tabular separated table output file

open(OUTFILE, ">>", $out_file) or die ("ERROR: in file " . __FILE__ ." at line "
                                     . __LINE__ ."\nFailed to open file $out_file for writing!\n");

if($out_exists == 0){
	print OUTFILE "x\ty\tz\ttype\tcell\n";
}

my @easy_types = ("cell_wall", "membrane", "nucleus_equivalent");
foreach my $et (@easy_types){
	if(exists $data{$et}){
		foreach my $l (@{$data{$et}}){
			print OUTFILE $l;
		}
	}
}

if(exists $data{"gold"}){
	while (my ($key, $value) = each (%{$data{"gold"}})){
		my $sumX = 0;
		my $sumY = 0;
		my $n = @{$$value{"x"}}."\n"; # must be the same for x and y implicitely
		# average x and y from the "stamp" coordinates for each gold particle
		for (my $i=0; $i < $n; $i++) {
       		$sumX += ${$$value{"x"}}[$i];
        	$sumY += ${$$value{"y"}}[$i];
    	}
    	print OUTFILE ($sumX/$n)."\t".($sumY/$n)."\t".$$value{"z"}."\tgold\t".$$value{"cell"}."\n";
    }
}

close(OUTFILE) or die ("ERROR: in file " . __FILE__ ." at line "
                                     . __LINE__ ."\nFailed to close file $out_file"."!\n");