#!/usr/bin/env python
import sys,time,argparse
import numpy as np

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	convert_gpd_to_fasta(args.input_fasta,args.input_gpd,args.output_fasta)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def convert_gpd_to_fasta(input_fa_fl,input_gpd_fl,output_fa_fl):
	# --- parse genome fasta file ---
	dic_chr_seq = {}
	chr = ""
	seq_list = [""]
	for line in input_fa_fl:
		if line.startswith(">"):
			dic_chr_seq[chr] = "".join(seq_list)
			seq_list = []
			chr = line.strip().split(" ")[0][1:]
		else:
			seq_list.append(line.strip())
	dic_chr_seq[chr] = "".join(seq_list)
	del dic_chr_seq[""]
	# --- parse gpd annotation and generate fasta for each isoform ---
	dic_gatc_com = {"G":"C","A":"T","T":"A","C":"G"}
	dic_nogatc_code = {"U":["T"],"I":["G"],"R":["A","G"],"Y":["C","T"],"S":["G","C"],"W":["A","T"],"K":["G","T"],"M":["A","C"],"B":["C","G","T"],"D":["A","G","T"],"H":["A","C","T"],"V":["A","C","G"],"N":["A","C","G","T"],"X":["A","C","G","T"]}
	for line in input_gpd_fl:
		# parse first isoform
		gene1,iso1,chr1,strand1,tss1,tts1,cds_start1,cds_end1,exon_number1,exon_start1,exon_end1 = line.strip().split("\t")[:11]
		seq_list1 = []
		for i in range(0,int(exon_number1)):
			start1 = int(exon_start1.split(",")[i])
			end1 = int(exon_end1.split(",")[i])
			seq_list1.append(dic_chr_seq[chr1][start1:end1])
		seq1 = "".join(seq_list1)
		seq1 = seq1.upper()
		up_seq_list1 = list(seq1) # non AGCT nucleotide
		for i in xrange(len(up_seq_list1)):
			if up_seq_list1[i] in dic_nogatc_code:
				up_seq_list1[i] = np.random.choice(dic_nogatc_code[up_seq_list1[i]])
		if strand1 == "-": # minus strand, get reverse complementary sequence
			fusion_site1 = str(int(tss1)+1)
			seq_com1 = ""
			for nuc in up_seq_list1:
				seq_com1 += dic_gatc_com[nuc]
			seq1 = seq_com1[::-1]
		else:
			seq1 = "".join(up_seq_list1)
			fusion_site1 = tts1
		# parse second isoform	
		gene2,iso2,chr2,strand2,tss2,tts2,cds_start2,cds_end2,exon_number2,exon_start2,exon_end2 = line.strip().split("\t")[11:]
		seq_list2 = []
		for i in range(0,int(exon_number2)):
			start2 = int(exon_start2.split(",")[i])
			end2 = int(exon_end2.split(",")[i])
			seq_list2.append(dic_chr_seq[chr2][start2:end2])
		seq2 = "".join(seq_list2)
		seq2 = seq2.upper()
		up_seq_list2 = list(seq2) # non AGCT nucleotide
		for i in xrange(len(up_seq_list2)):
			if up_seq_list2[i] in dic_nogatc_code:
				up_seq_list2[i] = np.random.choice(dic_nogatc_code[up_seq_list2[i]])
		if strand2 == "-": # minus strand, get reverse complementary sequence
			fusion_site2 = tts2
			seq_com2 = ""
			for nuc in up_seq_list2:
				seq_com2 += dic_gatc_com[nuc]
			seq2 = seq_com2[::-1]
		else:
			seq2 = "".join(up_seq_list2)
			fusion_site2 = str(int(tss2)+1)
		iso = "+".join([iso1,fusion_site1,iso2,fusion_site2])
		seq = seq1 + seq2
		print >>output_fa_fl, ">" + iso
		print >>output_fa_fl, seq

	input_fa_fl.close()
	input_gpd_fl.close()
	output_fa_fl.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Generate fasta sequence for each isoform. Specific for fusion transcript. GPD file with 22 columns, first 11 columns are 5'end fusion fragment information and last 11 columns are 3'end fragment information. Fusion transcript is named as 'Isoform ID (5'end) + fusion site in 5'end isoform  + Isoform ID (3'end) + fusion site in 3'end isoform", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a','--input_fasta',type=argparse.FileType('r'),required=True,help="Input: genome fasta file")
	parser.add_argument('-g','--input_gpd',type=argparse.FileType('r'),required=True,help="Input: annotation gpd (22 columns) file")
	parser.add_argument('-o','--output_fasta',type=argparse.FileType('w'),required=True,help="Output: isoform fasta file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
