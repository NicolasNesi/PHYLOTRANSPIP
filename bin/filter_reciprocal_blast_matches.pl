#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Fasta;
use Getopt::Long;

my $i = '';
my $r = '';
my $d = '';
my $s = '';
my $o = '';
GetOptions(
			'i=s' => \$i,
			'r=s' => \$r,
			'd=s' => \$d,
			's=s' => \$s,
			'o=s' => \$o
); 
		
# Open reciprocal blast file to read
my $reciprocal_blast_file = $r;
open(reciprocal_blast_file_data, '<', $reciprocal_blast_file) 
		or die "Failed to open input: $!";

my $reciprocal_blast_table={};
my $refseqID;
LINE: while (<reciprocal_blast_file_data>) {
		chomp;
			next LINE if ( /#/);
			my @data = split(/\s+/, $_);
     		my $refseqID = "$data[0]";
     		my $reciprocal_blast_hit = "$data[1]";
			${$reciprocal_blast_table}{$refseqID} = $reciprocal_blast_hit 
						unless exists ${$reciprocal_blast_table}{$refseqID};
		}

close(reciprocal_blast_file_data);

# Open infile to read
my $infile = $i;
open(blast_file_data, '<', $infile) 
		or die "Failed to open input: $!";
		
# Open outfile (outfile.out) to write
my $outfile = $o;
 open(outfile_data, '>', $outfile . '.fa')
 		or die "Failed to open input: $!";


 open(outfile_data_log, '>', $outfile . '.log')
 		or die "Failed to open input: $!";

# create database from directory of fasta files
my $sequence_data = $d;
my $db = Bio::DB::Fasta->new("$sequence_data");

my $species = $s;

LINE: while (<blast_file_data>) {
		chomp;
			next LINE if ( /#/);
			my @data = split(/\s+/, $_);
     		my $candidate_seqID = "$data[0]";
     		my $blast_hit = "$data[1]";
			if ( exists ${$reciprocal_blast_table}{$blast_hit} && 
				${$reciprocal_blast_table}{$blast_hit} eq $candidate_seqID ) {
						
						print outfile_data_log "$blast_hit matches $candidate_seqID in $species\n";
						
						my $subseq_start = $data[6];
						my $subseq_end = $data[7];
						my $subseq = $db->seq("$candidate_seqID:$subseq_start,$subseq_end");
						
						print outfile_data ">$species\|$blast_hit\|$candidate_seqID:$subseq_start-$subseq_end\n";
						print outfile_data "$subseq\n";

						}
		}

# Close infile and outfile
close(blast_file_data);
close(outfile_data);