#!/usr/bin/perl -w

use strict;
use warnings;
use Bio::SeqIO;
use Getopt::Long;

my $i = '';
my $o = '';
GetOptions(
			'i=s' => \$i,
			'o=s' => \$o
); 

# Open infile to read
my $infile = $i;
open(infile_data, '<', $infile) 
		or die "Failed to open input: $!";
		
# Open outfile (outfile.out) to write
 my $outfile = $o;
  open(outfile_data, '>', $outfile)
  		or die "Failed to open input: $!";


my $seqin = Bio::SeqIO->new( -format => 'fasta', -file => $infile);
#my $seqout = Bio::SeqIO->new( -format => 'fasta', -file => ">$outfile" );


while( (my $seq = $seqin->next_seq()) ) {
	my $taxon = ">" . $seq->id;
#	print outfile_data "$taxon\t";
	my $taxon_id = $seq->id;
	
 	my $pseq0 = $seq->translate(-frame => 0, -terminator => 'x', -unknown => 'X');
 	
 	my $taxon_aa0 = $pseq0->seq;
 	if ( $taxon_aa0 !~ m/x\w/ ) {
 		print outfile_data $taxon . "\n";
 		print outfile_data $seq->seq . "\n";
 		} 
 	
 }
		
# Close infile and outfile
close(outfile_data);
close(infile_data);
