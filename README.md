# MLST CLI

Make a database of wgMLST alleles and manipulate it

# Synopsis

## Download

Get an MLST scheme, e.g., Salmonella cgMLST v2 from Enterobase.  Uncompress the fasta files to a subfolder.

### Example from cgmlst.org

```
mkdir Salmonella_enterica.cgmlst
cd Salmonella_enterica.cgmlst
wget https://www.cgmlst.org/ncs/schema/4792159/alleles/ -O Salmonella_enterica.cgmlst.zip
# This zip uncompresses to $CWD
unzip Salmonella_enterica.cgmlst.zip
# Also fix the deflines to be in the standard format of locus_allele
for i in *.fasta; do 
  echo $i; 
  locus="$(basename $i .fasta)_"; 
  sed -i.bak "s/>/>$locus/" $i; 
done;

grep -m 2 ">" *.fasta
# The output of the grep command should show something like
#SC0831.fasta:>SC0831_1
#SC0831.fasta:>SC0831_2
#SEN0401.fasta:>SEN0401_1
#SEN0401.fasta:>SEN0401_2
#SNSL254_A1382.fasta:>SNSL254_A1382_1
#SNSL254_A1382.fasta:>SNSL254_A1382_2
#SPAB_04503.fasta:>SPAB_04503_1
#SPAB_04503.fasta:>SPAB_04503_2
#STM0834.fasta:>STM0834_1
#STM0834.fasta:>STM0834_2
# If it does not look correct, then undo with this loop!
# for i in *.fasta; mv -v $i.bak $i; done;

# If that looks correct then delete backup files
rm -v *.bak
cd -
```

## Make the database 

```
perl scripts/fasta-to-db.pl -o Salmonella_enterica.cgmlst.sqlite Salmonella_enterica.cgmlst/*.fasta
```

## Make a pseudoalignment

1. You need a spreadsheet of profiles whose columns are Name, allele1, allele2,...
2. run `profilesToMsa.pl`
