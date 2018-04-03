# ClassifyBySeq

Description: This program is classifying barcode+primer sequence to different sample.

## usage: python ClassifyBySeq.py [options]

optional arguments:
```
-h/--help show this help message and exit

-f1/--fqfile1 FQFILE1 Input fastq format file R1.

-f2/--fqfile2 FQFILE2 Input fastq format file R2.

-d/--outdir OUTDIR Output file director

-o/--outfile OUTFILE Output tab format file, which record search results.

-w/--window WINDOW [optional] Search Windowsize. Default=10

-m/--mismatch MISMATCH [optional] Tolerance of mismatches. Default=1

-c/--configure CONFIGURE Inquiry Sequences file, splited by tab.

-v/--version Show program's version number and exit
```

## conf-file format, splited by tab
```
A501P5  TGAACCTTAACGAATGGTATCAACTCAGAG
A502P5  TGCTAAGTAACHAATGGTATCAACTCAGAG
N701P3G  TAAGGCGAATGCAAGACCGATGGGATTTTGGTGG
N702P3G  CGTACTAATGCAAGACCGATGGGATTTTGGTGG
...
```
## Change Log
```
ClassifyBySeq.V1.py is calculate identity
ClassifyBySeq.V2.py is calculate mismatch number
```

Created by WangChengrui, 29th Marth,2018
