# ClassifyBySeq
Classify barcode+primer sequence to different sample

usage: my - program [options] usage

my - description

optional arguments:

  -h, --help            show this help message and exit
  
  -f1 FQFILE1, --fqfile1 FQFILE1
                        Input fastq format file R1.
                        
  -f2 FQFILE2, --fqfile2 FQFILE2
                        Input fastq format file R2.
                        
  -o OUTFILE, --outfile OUTFILE
                        Output tab format file, which record search results.
                        
  -w WINDOW, --window WINDOW
                        [optional] Search Windowsize. Default=10
                        
  -m MISMATCH, --mismatch MISMATCH
                        [optional] Tolerance of mismatches. Default=1
                        
  -c CONFIGURE, --configure CONFIGURE
                        Inquiry Sequences file, splited by tab. Like
                        
  -v, --version         show program's version number and exit
  

Created by WangChengrui, 29th Marth,2018
