#!/usr/bin/perl

=pod

=head1 NAME

progressive_blast_amplicon.pl

=head1 SYNOPSIS

progressive_blast_amplicon.pl -d database.fa -p primers.csv -i 1 -f 2 -r 3 -n 1 -l 500 -o report.csv -m matches.csv -u unmatched.csv 

=head1 DESCRIPTION

This script locates the genomic region that is most likely amplfied by a pair of primers
using progressively less stringent BLAST searches and checks for distance and correct orientation of the primers.

=head1 LICENSE

  Same as Perl.

=head1 AUTHORS

  Angela M. Baldo 
  Jeremy D. Edwards <jde22@cornell.edu>

=cut

use strict;
use warnings;
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::DB::Fasta;
use Bio::Tools::Run::StandAloneBlast;

our ($opt_d, $opt_p, $opt_i, $opt_f, $opt_r, $opt_n, $opt_l, $opt_o, $opt_m, $opt_u, $opt_h);
getopts('d:p:i:f:r:n:o:m:u:h');
if ($opt_h){
    help();
    exit;
}
if ((!$opt_i) || (!$opt_p) || (!$opt_i) || (!$opt_f) || (!$opt_r) || (!$opt_n) || (!$opt_l) || (!$opt_o) || (!$opt_m) || (!$opt_u)) {
    print STDERR "\nMissing required options\n\n\n";
    help();
}

# User options
my $database = $opt_d;    # FASTA file formatted for BLAST with formatdb
my $index = $opt_p;       # spreadsheet with the primers in it (marker[tab]primer1[tab]primer2)
my $ID_col = $opt_i;      # column where the marker name appears
my $primer1_col = $opt_f; # column where the primer1 appears
my $primer2_col = $opt_r; # column where the primer2 appears
my $max_matches = $opt_n; # maximum number of matches desired for a marker (usually 1)
my $max_length = $opt_l;  # maximum amplicon length
my $report = $opt_r;	  # report dump for looking at results
my $matches = $opt_m;	  # tsv of matches
my $unmatched = $opt_u;   # dump for unmatched primers

# Built in
my $spacer = "NNNNNNNNNNNNNNNNNNNN";       #desired spacer between the primers
my @evalues=(0.01, 0.1, 1, 10, 100, 1000); #array of e-values to try progressively

# Data storage
my %markers;        # hash of markers containing query sequences
my %positive_match; # list of markers for which there was at least one positive match
my %num_matches;    # list of markers and how many positive matches there were
my %matches;        # hash of markers and their output strings

# open output files
open(MATCHES,"> $matches")||die("can't open $matches: $!\n");
open(UNMATCHED,"> $unmatched")||die("can't open $unmatched: $!\n");
open(REPORT,"> $report")||die("can't open $report: $!\n");

# open input file
open(INDEX,$index)||die("can't open $index: $!\n");

# print column headers in output files
print MATCHES "Num_Hits\tMarker\tChrom\tSize\tRegion\tEval pos\tPrimer start pos\tPrimer end pos\tEval neg\tPrimer start neg\tPrimer end neg\tQuery string\n";
print UNMATCHED "Num_Hits\tMarker\tChrom\tSize\tRegion\tEval pos\tPrimer start pos\tPrimer end pos\tEval neg\tPrimer start neg\tPrimer end neg\tQuery string\n";
print REPORT "Eval_threshold\tMarker\tChrom\tNum HSPs\t\tQuery Range\tHit Strand\tEvalue\tHit Range\n";

my $linenum=0;
while (<INDEX>) {
    chomp $_;
    if ($linenum == 0) { # skip the header
	$linenum++;
	next;
    }
    my $temp = $_;
    $temp =~ s/"//g;    # strip the quotes
    $temp =~ s///g;   # strip those weird ^M characters that show up as line terminators

    #This is the order in which columns are expected in the input file:
    my @line = split(/\t/,$temp);
    my $marker_name = $line[$ID_col];
    my $primer1 = $line[$primer1_col];
    my $primer2 = $line[$primer2_col];

    # get rid of whitespace in any of these
    if (($marker_name ne "") && ($primer1 ne "") && ($primer2 ne "")) {
	$marker_name =~ s/\s//g;
	$primer1 =~ s/\s//g;
	$primer2 =~ s/\s//g;
    }

    # skip row if any of the fields are empty or invalid
    if (!(($marker_name ne "") && ($primer1 ne "") && ($primer2 ne "") && ($primer1 =~ /\A[acgtACGT]*\z/i) && ($primer1 =~ /\A[acgtACGT]*\z/i))) {
	print "Skipping amplicon $marker_name line $linenum\n";
	next;
    }

    # store concatenated string of primers
    my $string = $primer1.$spacer.$primer2;
    $markers{$marker_name}=$string;

    my $i = 0; # increment counter for evalues
    while ((not exists $positive_match{$marker_name}) && ($i <= $#evalues)) {

	# create sequence object from primer string
	my $query = Bio::Seq->new(-display_id => $marker_name, -seq => $string, -alphabet => 'dna');

	# setup BLAST
	my $factory = Bio::Tools::Run::StandAloneBlast->new('program'  => 'blastn',
							    'database' => $database,
							    'F' => 'F',
							    _READMETHOD => "Blast",
							    'e' => $evalues[$i]
	    );
	# Run BLAST
	my $blast_report = $factory->blastall($query);
	my $result = $blast_report->next_result; #only one BLAST search so only one result

	my $strandsum = 0;
	my @hsp_locations = ();
	my $match_flag = 0; # is this used?

	my %forwardhash;	# HSPs on strand 1
	my %reversehash;	# HSPs on strand -1

	while (my $hit = $result->next_hit()) { # loop through all hits

	    print REPORT "$evalues[$i]\t".$query->display_id()."\t",
	    $hit->accession,"\t",
	    $hit->num_hsps,"\t";
	    undef @hsp_locations;

	    $strandsum = 0;
	    foreach my $item ($hit->hsps) {
		if ($item->strand('hit') eq 1) {
		    $forwardhash{$item->start('hit')}=$item;
		} elsif ($item->strand('hit') eq -1) {
		    $reversehash{$item->end('hit')}=$item;
		}

		$strandsum = $strandsum + $item->strand('hit');
		push (@hsp_locations,$item->start('hit'),$item->end('hit'));
		print REPORT "\t",$item->start('query'),"-",$item->end('query'),
		"\t",$item->strand('hit'),
		"\t",$item->evalue,
		"\t",$item->start('hit'),"-",$item->end('hit'),"\t",
	    }		# while parsing hsps
	    print REPORT "\n";

	    foreach my $fwd (sort {$a<=>$b} keys %forwardhash) {
		foreach my $rev (sort {$b<=>$a} keys %reversehash) {
		    if (($fwd <= $rev) && ($rev-$fwd <= $max_length)) {

			$num_matches{$query->display_id()}++;
			my $IDandnum_matches = $query->display_id()."\t".$num_matches{$query->display_id()};
			my $size = $rev-$fwd;
			$matches{$IDandnum_matches} =
			    $query->display_id()."\t".		   # Marker
			    $hit->accession."\t".		   # Chrom
			    $size."\t".			   # size
			    "$fwd"."-"."$rev"."\t".		   # region
			    $forwardhash{$fwd}->evalue."\t". # left primer eval
			    $forwardhash{$fwd}->start('query')."\t". # left primer start position
			    $forwardhash{$fwd}->end('query')."\t". # left primer size
			    $reversehash{$rev}->evalue."\t". # right primer eval
			    $reversehash{$rev}->start('query')."\t". # right primer start position
			    $reversehash{$rev}->end('query'). # right primer size
			    "\t$string"; # query string

			$positive_match{$query->display_id()}=$hit->accession;
			$match_flag = 1;  # is this used?
		    }
		}
	    }

	    undef %forwardhash;
	    undef %reversehash;

	}			# while parsing hits
	$i++;
    }		 # while incrementing $i until we get a positive match
    $linenum++;
}
close (INDEX);
close (REPORT);

foreach my $marker (sort keys %num_matches) {
  if ($num_matches{$marker} <= $max_matches) {

    foreach my $match (sort keys %matches) {
      if ($match =~ /$marker\t/) {
	print MATCHES "$num_matches{$marker}\t$matches{$match}\n";
      }
    }
  } else {
    foreach my $match (sort keys %matches) {
      if ($match =~ /$marker\t/) {
	print UNMATCHED "$num_matches{$marker}\t$matches{$match}\n";
      }
    }
  }
}

foreach my $marker (sort keys %markers) {
  if (not exists $positive_match{$marker}) {
    print  UNMATCHED "0\t$marker\t\t\t\t\t\t\t\t\t\t$markers{$marker}\n";
  }
}
close(MATCHES);
close(UNMATCHED);




sub help {
  print STDERR <<EOF;
  $0:

    Description:

      This script locates the genomic region that is most likely amplfied by a pair of primers
      using progressively less stringent BLAST searches and checks for distance and correct orientation of the primers.

    Usage:

      progressive_blast_amplicon.pl -d database.fa -p primers.csv -i 1 -f 2 -r 3 -n 1 -l 500 -o report.csv -m matches.csv -u unmatched.csv 

    Flags:

      -d <database_file>                 FASTA database file. Must be formatted with formatdb. (mandatory)

      -p <primer_file>                   Primer file.  Comma separated. (mandatory)

      -i <id_column_number>              Column number of amplicon names/IDs in the primer file. (mandatory)

      -f <forward_primer_column_number>  Column number of forward primer in the primer file (mandatory)

      -r <reverse_primer_column_number>  Column number of reverse primer in the primer file (mandatory)

      -n <max_matches>                   Maximum number of matches (mandatory)

      -l <max_length>                    Maximum amplicon length in bp (mandatory)

      -o <report_output_file>            Output file name for the report (mandatory)

      -m <matches_output_file>           Output file name for the matched amplicons (mandatory)

      -u <unmatched_output_file>         Output file name for the unmatched amplicons (mandatory)

      -h <help>                          Help

EOF
exit (1);
}
