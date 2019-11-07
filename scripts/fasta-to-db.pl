#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper;
use Getopt::Long;
use File::Basename qw/basename/;
use DBD::SQLite;
use Bio::SeqIO;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(outfile=s help)) or die $!;
  die usage() if($$settings{help});
  die "ERROR: need --outfile or -o" if(!$$settings{outfile});
  die "ERROR: need locus.fasta" if(!@ARGV);

  my $dbh = createDb($$settings{outfile}, $settings);
  addToDatabase($dbh, \@ARGV, $settings);

  return 0;
}

sub createDb{
  my($filename, $settings) = @_;

  my $dbh = DBI->connect("dbi:SQLite:dbname=$filename","","");
  $dbh->do(qq(
    CREATE TABLE allele(
    locus  varchar(20), /* longest was 17 bytes */
    allele INT,
    seq    TEXT
    );
  ));

  $dbh->do(qq(
    CREATE INDEX locus on allele (locus);
  ));
  $dbh->do(qq(
    CREATE INDEX locus_allele on allele (locus,allele);
  ));
  $dbh->do(qq(
    CREATE INDEX allele_field on allele (allele);
  ));
  
  return $dbh;
}

sub addToDatabase{
  my($dbh, $alleleFiles, $settings) = @_;
  
  my $insertCounter = 0;
  my $prepareStr = "INSERT INTO allele (locus, allele, seq) VALUES ";
  my $valuesStr  = "";
  my @values;
  for my $filename(@$alleleFiles){
    logmsg "  Reading $filename ...";
    my $in = Bio::SeqIO->new(-file=>$filename);
    while(my $seq = $in->next_seq){
      if($seq->id !~ /(.+)[-_](\d+)/){
        logmsg "WARNING: not parsing ".$seq->id;
        next;
      }
      my($locus, $allele) = ($1,$2);
      $valuesStr .= '(?,?,?),';
      push(@values, $locus, $allele, $seq->seq);

      $insertCounter++;
      if($insertCounter % 100 == 0){
        $valuesStr =~ s/,\s*$//; # remove last comma
        my $sth = $dbh->prepare("$prepareStr $valuesStr");
        $sth->execute(@values);
        # reset
        $valuesStr = "";
        @values = ();
        # Progress
        if($insertCounter % 100000 == 0){
          logmsg "    Executed $insertCounter inserts";
          #logmsg "DEBUG"; last;
        }
      }
    }
    $in->close;
  }

  # One last insert
  if($valuesStr){
    $valuesStr =~ s/,\s*$//; # remove last comma
    my $sth = $dbh->prepare("$prepareStr $valuesStr");
    $sth->execute(@values);
  }
  logmsg "$insertCounter inserts were applied to the database";
  return $insertCounter;
}

sub usage{
  "$0: turns a multifasta file of MLST alleles into a database
  Usage: $0 -o loci.sqlite loci.fasta [loci2.fasta...]
  loci.fasta      A fasta file whose deflines are in the format of 
                  >locus-allele or >locus_allele
  "
}
