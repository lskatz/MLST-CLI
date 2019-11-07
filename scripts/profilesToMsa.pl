#!/usr/bin/env perl 

use warnings;
use strict;
use Data::Dumper qw/Dumper/;
use Getopt::Long qw/GetOptions/;
use File::Basename qw/basename dirname/;
use File::Temp qw/tempdir tempfile/;
use File::Slurp qw/read_file/;
use DBD::SQLite;

local $0 = basename $0;
sub logmsg{local $0=basename $0; print STDERR "$0: @_\n";}
exit(main());

sub main{
  my $settings={};
  GetOptions($settings,qw(help tempdir=s profiles=s db|database=s)) or die $!;
  die usage() if($$settings{help});
  $$settings{profiles} //= die "ERROR: need --profiles";
  $$settings{db} //= die "ERROR: need --database";
  $$settings{tempdir} //= tempdir("$0.XXXXXX", TMPDIR=>1, CLEANUP=>1);

  logmsg "Reading MLST profiles";
  my $genomeAlleles = readProfile($$settings{profiles}, $settings);
  logmsg "Reading database...";
  my $dbh = openDatabase($$settings{db}, $settings);

  # Turn those alleles and profiles into an MSA.
  # Do some local alignment within loci.
  my $msa = makeMsa($genomeAlleles, $dbh, $settings);

  # Print to alignment
  while(my($sample,$sequence) = each(%$msa)){
    $sequence=~s/\s+//g;
    print ">$sample\n$sequence\n";
  }

  return 0;
}

# Make an alignment from alleles.
# Return a hash of genome=>gappedSeq
sub makeMsa{
  my($alleles, $dbh, $settings) = @_;

  my %pseudoAlignmentSeq; # genome => sequence

  # Get an index of allele names in the whole set.
  my %alleleCount;
  my @genomeName;
  while(my($genome, $locusHash) = each(%$alleles)){
    push(@genomeName,$genome);
    while(my($locus,$allele) = each(%$locusHash)){
      $alleleCount{$locus}++;
    }
  }
  my @locusName = sort{$b cmp $a} keys(%alleleCount);

  for my $locus(@locusName){
    my %allele; # genome => allele-sequence
    my $alleleCounter++;
    for my $genome(@genomeName){
      my $alleleId = $$alleles{$genome}{$locus} || 0;

      # If the allele is a whole number, then query for the sequence.
      if($alleleId){
        my $sth = $dbh->prepare(qq(
          SELECT seq 
          FROM allele 
          WHERE locus=?
            AND allele=?
          LIMIT 1 /* should not be necessary, but just in case */
        ));
        my $res = $sth->execute($locus, $alleleId);
        my ($sequence) = $sth->fetchrow_array();
        # If the sequence is blank, then toss this allele.
        if(!defined($sequence) || $sequence !~ /[atcgATCG]/){
          next;
        }
        $allele{$alleleId} = $sequence;
      }
    }

    # Multiple sequence alignment with muscle
    my $gappedSeq = alignLoci(\%allele, $settings);
    # Length of each allele in the alignment
    my $firstSeq = (values(%$gappedSeq))[0];
    my $gappedLength = length($firstSeq);
    my $allGaps = '-' x $gappedLength;
    
    # Append to the whole alignment
    for my $genome(@genomeName){
      # Get the allele for this genome
      my $alleleId = $$alleles{$genome}{$locus};
      # What the gapped sequence is, or if it's all gaps
      my $sequence = $$gappedSeq{$alleleId} || $allGaps;
      # Append.
      $pseudoAlignmentSeq{$genome} .= $sequence;
    }
    #die Dumper \%pseudoAlignmentSeq if(length($pseudoAlignmentSeq{PNUSAS080689}) > 1000000);
    print STDERR ".";
  }
  print STDERR "\n";

  return \%pseudoAlignmentSeq;
}

sub alignLoci{
  my($seqs, $settings) = @_;
  my $numSeqs = scalar(keys(%$seqs));

  my %gappedSeq;  # allele =>  seq

  # Shortcut: do not align if we have fewer than two entries.
  if($numSeqs < 2){
    # If we have only one sequence, then set the gapped hash
    # of sequences to simply the first sequence.
    if($numSeqs == 1){
      %gappedSeq = %$seqs;
    }
    # If we have no sequences then return the null allele
    else {
      $gappedSeq{0} = '-';
    }
    # Return what is essentially ungapped/unaligned sequence.
    return \%gappedSeq;
  }

  # Make input file
  my($fh, $filename) = tempfile("locus_XXXXXX", SUFFIX=>".fasta", DIR=>$$settings{tempdir});
  while(my($alleleId, $sequence) = each(%$seqs)){
    print $fh ">$alleleId\n$sequence\n";
  }
  close($fh);

  # Alignment.
  # semi-fast options for muscle, output into fasta format
  system("muscle -in $filename -out $filename.aln -maxiters 5 -diags > $filename.muscle.log 2>&1");
  die "ERROR with muscle:\n".`grep ">" $filename`."\n".`cat $filename.muscle.log` if $?;

  my $alignment = read_file("$filename.aln");
  while($alignment =~ />(\S+).*?\n([^>]+)/sg){
    my $seqId = $1;
    my $seqSeq= $2;
    $gappedSeq{$seqId} .= $seqSeq;
  }

  return \%gappedSeq;
}

# Reads the tab-separated-values profile
# Returns hash of genome=>{locus=>allele}
sub readProfile{
  my($profile, $settings)=@_;

  my %genomeAllele;

  open(my $fh, $profile) or die "ERROR: could not read $profile: $!";
  my $header = <$fh>;
  chomp($header);
  my (undef, @header) = split(/\t/, $header);
  while(<$fh>){
    chomp;
    my ($genome, @allele) = split(/\t/, $_);
    # make a hash of locus => allele
    my %alleleHash;
    @alleleHash{@header} = @allele;

    $genomeAllele{$genome} = \%alleleHash;
  }
  close $fh;

  return \%genomeAllele;
}

sub openDatabase{
  my($filename, $settings) = @_;
  
  my $dbh = DBI->connect("dbi:SQLite:dbname=$filename","","");

  return $dbh;
}

sub usage{
  "$0: Outputs a multiple sequence alignment, given a multifasta file of alleles and a tsv of alleles
  Usage: $0 --profiles file.tsv [options] --database something.sqlite
  --profiles      A tab-separated file of profiles, where the first column is
                  the genome name and the next columns are loci. Missing loci
                  will be represented by gaps (-).
                  The first line is a header. The first column's header is
                  'Name' and the subsequent headers are locus names.
  --database      A database file created by fasta-to-db.pl
  "
}
