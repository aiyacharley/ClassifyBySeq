#!/zzh_gpfs/apps/python/bin/python
from __future__ import division
import os,sys,csv
from Bio import SeqIO
import gzip
import argparse
import Levenshtein

def open_fastq_file(filename):
	if filename.endswith('gz'):
		fastq_handle = gzip.open(filename,'rU')
	else:
		fastq_handle = open(filename,'rU')
	return fastq_handle
def smooth_step(seq1,seq2):
	mylist,myseq = [],[]
	num = len(seq2)-len(seq1)+1
	for i in range(num):
		misnum = Levenshtein.hamming(seq1,seq2[i:i+len(seq1)])
		myseq.append(seq2[i:])
		mylist.append(misnum)
	if min(mylist)<=5:
		return min(mylist),myseq[mylist.index(min(mylist))][:len(seq1)]
	else:
		return 100,"NNNNNNNNNN"
def main():
	fq1 = SeqIO.parse(fqfile1,'fastq')
	fq2 = SeqIO.parse(fqfile2,'fastq')
	pdict = {}
	infseq = csv.reader(open(configure,'rU'),delimiter='\t')
	for rec in infseq:
		if rec[0].startswith("#"):
			continue
		if rec[0].startswith("A"):
			pdict[rec[1][2:]] = rec[0]
		elif rec[0].startswith("N"):
			pdict[rec[1][2:]] = rec[0]
		else:
			pdict[rec[1]] = rec[0]
	out = csv.writer(open("%s/%s"%(outdir,outname),'wb'),delimiter='\t')
	num0 = 0
	try:
		while True:
			num0+=1
			if num0 % 1000 == 0:
				print "rec %d process ... "%num0
			rec1 = fq1.next()
			rec2 = fq2.next()
			list1,list2 = [],[]
			seqlist1,seqlist2 = [],[]
			for seq in pdict.keys():
				ide1,seq1 = smooth_step(seq,str(rec1.seq)[:len(seq)+window])
				ide2,seq2 = smooth_step(seq,str(rec2.seq)[:len(seq)+window])
				list1.append(ide1)
				list2.append(ide2)
				seqlist1.append(seq1)
				seqlist2.append(seq2)
			Tseq1 = seqlist1[list1.index(min(list1))]
			Tseq2 = seqlist2[list2.index(min(list2))]
			Rlabel1 = pdict[list(pdict.keys())[list1.index(min(list1))]]
			Rlabel2 = pdict[list(pdict.keys())[list2.index(min(list2))]]
			if min(list1)==100 and min(list2)==100:
				result = [rec1.id,min(list1),min(list2),"-","-",Tseq1,Tseq2]
			elif min(list1)==100 and min(list2)!=100:
				if "P5" in Rlabel2:
					result = [rec1.id,min(list2),min(list1),Rlabel2,"-",Tseq2,Tseq1]
				else:
					result = [rec1.id,min(list1),min(list2),"-",Rlabel2,Tseq1,Tseq2]
			elif min(list1)!=100 and min(list2)==100:
				if "P5" in Rlabel1:
					result = [rec1.id,min(list1),min(list2),Rlabel1,"-",Tseq1,Tseq2]
				else:
					result = [rec1.id,min(list2),min(list1),"-",Rlabel1,Tseq2,Tseq1]
			else:
				if "P5" in Rlabel1:
					result = [rec1.id,min(list1),min(list2),Rlabel1,Rlabel2,Tseq1,Tseq2]
				else:
					result = [rec1.id,min(list2),min(list1),Rlabel2,Rlabel1,Tseq2,Tseq1]
			out.writerow(result)
			misnum1 = result[1]
			misnum2 = result[2]
			if misnum1 <= mismatch and misnum2 <= mismatch:
				out1 = open('%s/%s_%s_R1.fastq'%(outdir,result[3],result[4]),'a')
				out2 = open('%s/%s_%s_R2.fastq'%(outdir,result[3],result[4]),'a')
				SeqIO.write(rec1,out1,'fastq')
				SeqIO.write(rec2,out2,'fastq')
			elif misnum1 <= mismatch and misnum2 > mismatch:
				out1 = open('%s/%s_NO_R1.fastq'%(outdir,result[3]),'a')
				out2 = open('%s/%s_NO_R2.fastq'%(outdir,result[3]),'a')
				SeqIO.write(rec1,out1,'fastq')
				SeqIO.write(rec2,out2,'fastq')
			elif misnum1 > mismatch and misnum2 <= mismatch:
				out1 = open('%s/NO_%s_R1.fastq'%(outdir,result[4]),'a')
				out2 = open('%s/NO_%s_R2.fastq'%(outdir,result[4]),'a')
				SeqIO.write(rec1,out1,'fastq')
				SeqIO.write(rec2,out2,'fastq')
			else:
				out1 = open('%s/NO_NO_R1.fastq'%outdir,'a')
				out2 = open('%s/NO_NO_R2.fastq'%outdir,'a')
				SeqIO.write(rec1,out1,'fastq')
				SeqIO.write(rec2,out2,'fastq')
				
	except StopIteration:
		print "Iteration Done !"

if __name__=='__main__':
	parser = argparse.ArgumentParser(prog='my - program',usage='%(prog)s [options]',description = 'Classify barcode to different samples',epilog = 'Created by WangCR')
	parser.add_argument('-f1','--fqfile1',help='Input fastq format file R1.')
	parser.add_argument('-f2','--fqfile2',help='Input fastq format file R2.')
	parser.add_argument('-d','--outdir',default=".",help='Output file director')
	parser.add_argument('-o','--outfile',help='Output tab format file, which record search results.')
	parser.add_argument('-w','--window',type=int,default=10,help="[optional] Search Windowsize. Default=10")
	parser.add_argument('-m','--mismatch',type=int,default=1,help="[optional] Tolerance of mismatches. Default=1")
	parser.add_argument('-c','--configure',type=str,help="Inquiry Sequences file, splited by tab. Like")
	parser.add_argument('-v','--version', action='version', version='Copyright (c) 29/3/2018, created by WangChengrui, version 1.0')
	args = parser.parse_args()
	fqfile1 = args.fqfile1
	fqfile2 = args.fqfile2
	outdir = args.outdir
	outname = args.outfile
	window = args.window
	mismatch = args.mismatch
	configure = args.configure
	os.system("mkdir -p %s"%outdir)
	main()
