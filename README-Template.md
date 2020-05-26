# Functional codon classification for synonymous genetic variants in breast cancer

## Description

This is a pipeline for annotation of changes in codon sequences and ΔRSCU values into the INFO field of variant calling files (.vcf)

Programming languages used:
* Python (v3.8.2)
* Perl
* Bash (v5.0)

Available software used:
* CodonW (v1.4.4)
* VEP (v100)

### Prerequisites

VEP requires:
* gcc, g++ and make
* Perl version 5.0 or above
* Perl packages: Archive::Zip, DBD::mysql, DBI

## Codon usage analysis using CodonW
### CodonW installation

In order to install CodonW you need to have administrative access
```
sudo apt-get install -y codonw
```
###Obtaining open reading frames

To obtain coding sequences, UCSC Table Browser was used in the series set of steps and named "GencodeSeq"

```
https://genome.ucsc.edu/ → Tools → Table browser → clade:Mammal, genome:Human,  assembly:Dec.2013(GRCh38/hg38), group:Genes and Gene Predictions, track:All GENCODE V27, region:genome, filter(mark cdsStartStat:cmpl and cdsEndStat:cmpl), outputformat:sequence, file type returned:gzip compressed → get output → submit → Sequence Formatting Options: All upper case → get sequence
```
In order to make sequence headers only their "ENST" IDs, use a custom script named only_enst.py

```
#! /usr/bin/python3

import sys

with open(sys.argv[1], 'r') as fin:
    for line in fin:
        if line.startswith(">"):
            line = line.split(" ")
            line = line[0].split("_")
            header = line[2].rstrip()
            print('>{}'.format(header))
        else:
            sequence = line.rstrip()
            print(sequence)
```
Usage example
```
./only_enst.py GencodeSeq > GencodeSeq_enst
```
To make a single line fasta file from an interleaved one, run multi_to_single.py
```
#! /usr/bin/python3
import sys
with open(sys.argv[1], 'r') as fin, open(sys.argv[2], 'w') as fout:
    block = []
    for line in fin:
        if line.startswith('>'):
            if block:
                fout.write(''.join(block) + '\n')
                block = []
            fout.write(line)
        else:
            block.append(line.strip())

    if block:
        fout.write(''.join(block) + '\n')
```
Usage example
```
./multi_to_single.py GencodeSeq_enst GencodeSeq_single
```
Rename the file for more convenient usage
```
mv GencodeSeq_single GenSeq
```
In order to have only one representative sequence per gene, match longest transcripts with gene IDs
```
use warnings;
open (ID, "<$ARGV[0]");
while ($line = <ID>) {
    chomp $line;
    if ($line =~ /^(ENST\S+)\t.+(ENSG\S+)$/) {
	$transcript2gene{$1} = $2;
	$maxLengths{$2}[0] = "";
	$maxLengths{$2}[1] = 0;
    }
}
close ID;

open (CDS, "<$ARGV[1]");
$line = <CDS>;
if ($line =~ /^\>\S+(ENST\w+\.\d+)\s/) {
    $transcript = $1;
}
$seq = "";
while ($line = <CDS>) {
    if ($line =~ /^\>\S+(ENST\w+\.\d+)\s/) {
	$len = length($seq);
	$seq = "";
	if (exists $transcript2gene{$transcript} && $len > $maxLengths{$transcript2gene{$transcript}}[1]) {
	    $maxLengths{$transcript2gene{$transcript}}[0] = $transcript;
	    $maxLengths{$transcript2gene{$transcript}}[1] = $len;
	}
	$transcript = $1;
    } else {
	chomp $line;
	$seq .= $line;
    }
}
$len = length($seq);
if (exists $transcript2gene{$transcript} && $len > $maxLengths{$transcript2gene{$transcript}}[1]) {
    $maxLengths{$transcript2gene{$transcript}}[0] = $transcript;
    $maxLengths{$transcript2gene{$transcript}}[1] = $len;
}
close CDS;
```

Extract longest transcripts from the whole sequence set
```
awk 'NR==FNR{n[">"$0];next} f{print f ORS $0;f=""} $0 in n{f=$0}' GENCODE_V27_IDs_longest_CDS.txt GenSeq > GenSeq_longest
```

Check basic sequence integrity using CodonW and output warnings into a seperate file (better to run this command in the background if the sequence file is large in size)
```
nohup codonw GenSeq_longest -nomenu 2> warnings.txt
```

Extract IDs from the warning list
```
cat warnings.txt | grep "Warning" | cut -d " " -f 4 | cut -c2- > ids_to_remove.txt
```
Remove sequences that did not pass the integrity check
```
wk 'BEGIN{while((getline<"ids_to_remove.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' GenSeq_longest > GenSeq_filtered
```

### Analysis of codon usage

To get RSCU values of high and low bias data sets, run correspondence analysis with CodonW (better to run in the background)
```
nohup codonw GenSeq_filtered -coa_cu -nomenu -silent
```
Calculate the overall RSCU values for the whole data set
```
nohup codonw -cutot GenSeq_filtered -nomenu
```


## Annotation of codon sequence changes using VEP
### VEP installation

Download
```
curl -L -O https://github.com/Ensembl/ensembl-vep/archive/release/100.zip
unzip 100.zip
cd ensembl-vep-release-100/
```
Install (use the current release of human genome assembly "homo_sapiens_GRCh38")
```
perl INSTALL.pl
```
### Annotation
Annotate changes in codon sequences into a seperate file
```
vep -i input.vcf -o output_ann.txt -cache
```
Filter for synonymous variants using an integrated filtering option in VEP
```
filter_vep -i output_ann.txt -o codons_filtered -filter "Consequence is synonymous_variant or splice_region_variant&synonymous_variant"
```

## Annotation of codon changes into .vcf

Add changes in the codon sequences into the INFO field of the .vcf file using a custom script codon_annotation.py
```
! /usr/bin/python3
'''
Description: the first part of the script makes a dictionary of chromosomal positions as keys
and codon changes as values. The second part uses this dictionary to annotate synonymous variants when matched
by chromosomal position in the dictionary

'''

import sys
import re
with open(sys.argv[1], 'r') as syn_ann,open(sys.argv[2], 'r') as vcf, open(sys.argv[3], 'w') as fout:
    codons_dictionary= {}
    for line in syn_ann:
        if line[0] != "#":
            line = line.split("\t")
            chrom_pos = line[1]
            chrom_pos = chrom_pos.split(":")
            chrom_pos = chrom_pos[1]
            codons = line[11]
            codons_dictionary[chrom_pos]=codons
    for line in vcf:
        if line.startswith("##"):
            print(line.strip("\n"),file=fout)
        elif line.startswith("#CHROM"):
            print(('##INFO=<ID=CODCH,Number=.,Type=String,Description="Changes in the codon sequence of synonymous variants: ref/alt">'),file=fout)
            print(line.strip("\n"),file=fout)
        elif line[0] != "#" and re.search('synonymous_variant', line):
            line = line.strip("\n")
            line = line.split("\t")
            chrom = line[0]
            pos = line[1]
            id = line[2]
            ref = line[3]
            alt = line[4]
            qual = line[5]
            filter = line[6]
            info = line[7]
            format = line[8]
            end = line[9]
            if pos not in codons_dictionary:
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t'.format(chrom, pos, id, ref, alt, qual, filter, info, format, end), file=fout)
            if pos in codons_dictionary:
                print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}|CODCH={}\t{}\t{}\t'.format(chrom, pos, id, ref, alt, qual, filter, info, codons_dictionary[pos], format, end), file=fout)
        else:
            print(line.strip("\n"),file=fout)
```
Usage example
```
./codon_annotation.py codons_filtered.txt input.vcf output.vcf
```
## Annotation of ΔRSCU values into the .vcf file

In order to add the difference of relative synonymous codon usage values between reference and alternative codons of synonymous variants, run a custom script rscu_annotation.py

```
#! /usr/bin/python3

rscu_dictionary = dict(UUU = 0.92, UUC = 1.08, UUA = 0.46, UUG = 0.77,
                    CUU = 0.80, CUC = 1.16, CUA = 0.43, CUG = 2.38,
                    AUU = 1.09, AUC = 1.39, AUA = 0.52, AUG = 1.00,
                    GUU = 0.73, GUC = 0.95, GUA = 0.47, GUG = 1.85,
                    UCU = 1.12, UCC = 1.30, UCA = 0.92, UCG = 0.33,
                    CCU = 1.15, CCC = 1.30, CCA = 1.09, CCG = 0.46,
                    ACU = 1.00, ACC = 1.41, ACA = 1.15, ACG = 0.45,
                    GCU = 1.05, GCC = 1.61, GCA = 0.91, GCG = 0.44,
                    UAU = 0.89, UAC = 1.11, UAA = 0.84, UAG = 0.66,
                    CAU = 0.84, CAC = 1.16, CAA = 0.53, CAG = 1.47,
                    AAU = 0.95, AAC = 1.05, AAA = 0.88, AAG = 1.12,
                    GAU = 0.93, GAC = 1.07, GAA = 0.85, GAG = 1.15,
                    UGU = 0.92, UGC = 1.08, UGA = 1.50, UGG = 1.00,
                    CGU = 0.48, CGC = 1.11, CGA = 0.65, CGG = 1.23,
                    AGU = 0.90, AGC = 1.44, AGA = 1.26, AGG = 1.27,
                    GGU = 0.64, GGC = 1.35, GGA = 1.00, GGG = 1.01)

import sys
import re
with open(sys.argv[1], 'r') as vcf,open(sys.argv[2], 'w') as fout:
    for line in vcf:
        if line.startswith("##"):
            print(line.strip("\n"),file=fout)
        elif line.startswith("#CHROM"):
            print(('##INFO=<ID=deltaRSCU,Number=1,Type=Float,Description="Difference in relative synonymous codon usage values between alternative and reference codons"'),file=fout)
            print(line.strip("\n"),file=fout)
        elif line[0] != "#" and re.search("CODCH",line):
            line = line.strip("\n")
            line = line.split("\t")
            chrom = line[0]
            pos = line[1]
            id = line[2]
            ref = line[3]
            alt = line[4]
            qual = line[5]
            filter = line[6]
            info = line[7]
            format = line[8]
            end = line[9]
            codons = info.split('CODCH')[1]
            ref_codon = codons[1:4].upper().replace('T','U')
            alt_codon = codons[5:8].upper().replace('T','U')
            rscu_value = round((rscu_dictionary[alt_codon]-rscu_dictionary[ref_codon]),2)
            print('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}|deltaRSCU={}\t{}\t{}\t'.format(chrom, pos, id, ref, alt, qual, filter, info, rscu_value, format, end), file=fout)
        else:
            print(line.strip("\n"),file=fout)

```
Usage example
```
./rscu_annotation.py input.vcf output.vcf
```

## Creator and contact

Monika Kurgonaite – monika.kurg@gmail.com

## Acknowledgments

I wish to thank Helena Persson, whose code was used for sequence filtering in codon usage analysis step
