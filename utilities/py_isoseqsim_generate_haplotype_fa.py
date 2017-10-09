#!/usr/bin/env python
import sys,re,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	generate_haplotype_fasta(args.genome,args.haplotype,args.id,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def generate_haplotype_fasta(input_fa_fl,input_vcf_fl,individual_id,output_fa_fl):
	# --- parse vcf ---
	base4 = ["A","T","C","G"]
	dic_vcf_info = {}
	for line in input_vcf_fl:
		if line.startswith("##"):
			pass
		elif line.startswith("#CHROM\t"):
			individual_id_list = line.rstrip("\n").split("\t")[9:]
			if individual_id in individual_id_list:
				iv_id_index = individual_id_list.index(individual_id)
			else:
				sys.stderr.write("Individual ID does not exist, please check your VCF file")
				sys.exit()
		else:
			chrom,pos,id,ref,alt_set,qual,filter,info,format_set = line.rstrip("\n").split("\t")[:9]
			if ref not in base4:
				continue
			alt_flag = 0
			for alt in alt_set.split(","):
				if alt not in base4:
					alt_flag = 1
					break
			if alt_flag == 1:
				continue
			snp =  ref + "," + alt_set
			format_index = format_set.split(":").index("GT")
			iv_gt = line.rstrip("\n").split("\t")[9:][iv_id_index].split(":")[format_index]
			if "|" in iv_gt and "/" not in iv_gt:
				vcf_info = pos + "\t" + snp + "\t" + iv_gt
				ploidy = len(iv_gt.split("|"))
				if chrom not in dic_vcf_info.keys():
					dic_vcf_info[chrom] = []
					dic_vcf_info[chrom].append(vcf_info)
				else:
					dic_vcf_info[chrom].append(vcf_info)
	input_vcf_fl.close()
#	sys.stdout.write("The number of ploid is "+str(ploidy)+"\n")

	# --- parse genome file ---
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
	chr_list = dic_chr_seq.keys()
	chr_list.sort() # sort chromosome
	input_fa_fl.close()

	# --- generate haplotypes ---
	for chr in chr_list:
		seq = dic_chr_seq[chr].upper()
		if chr in dic_vcf_info.keys():
			dic_ploidy = {}
			for i in range(0,ploidy):
				dic_ploidy[i] = list(seq)
			for vcf_info in dic_vcf_info[chr]:
				pos,snp,iv_gt = vcf_info.split("\t")
				for i in range(0,ploidy):
					dic_ploidy[i][int(pos)-1] = snp.split(",")[int(iv_gt.split("|")[i])]
			for i in range(0,ploidy):
				print >>output_fa_fl, ">" + chr + "\tAllele" + str(i+1)
				new_seq = "".join(dic_ploidy[i])
				for j in range(0,len(new_seq),80):
					print >>output_fa_fl, new_seq[j:j+80]
		else:
			for i in range(0,ploidy):
				print >>output_fa_fl, ">" + chr + "\tAllele" + str(i+1)
				new_seq = seq
				for j in range(0,len(new_seq),80):
					print >>output_fa_fl, new_seq[j:j+80]
	output_fa_fl.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Generate fasta file for each haplotype of multi-ploid samples. Note: VCF file must have the header lines '#CHROM' for showing the individual ID to be used; Only the phased (split by '|' not '/' in GT format of VCF file) single nucleotide substitution site is used.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-g','--genome',type=argparse.FileType('r'),required=True,help="Input: genome fasta file. Genome sequence")
	parser.add_argument('-v','--haplotype',type=argparse.FileType('r'),required=True,help="Input: genotype vcf file. Genotype")
	parser.add_argument('-i','--id',type=str,required=True,help="Individual/Sample ID located in the header line startswith '#CHROM'")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: multi-ploid fasta file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
