# MLST CLI

Make a database of wgMLST alleles and manipulate it

# Synopsis

## Download

Get an MLST scheme, e.g., Salmonella cgMLST v2 from Enterobase.  Uncompress the fasta files to a subfolder.

## Make the database 
    fasta-to-db.pl --outfile Salmonella.cgMLSTv2.sqlite Salmonella.cgMLSTv2.fasta

## Make a pseudoalignment

1. You need a spreadsheet of profiles whose columns are Name, ST, allele1, allele2,...
2. run `profilesToMsa.pl`
