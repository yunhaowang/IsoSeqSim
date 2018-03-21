#!/usr/bin/env python
import sys,time,argparse
import numpy as np

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	generate_ats_event(args.input,args.output,args.genome,args.distance)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def generate_ats_event(input_gpd_fl,output_gpd_fl,genome_fa_fl,dis_ts):
	# parse genome fasta file
	dic_chr_size = {}
	chr = ""
	seq_list = [""]
	for line in genome_fa_fl:
		if line.startswith(">"):
			dic_chr_size[chr] = len("".join(seq_list))
			seq_list = []
			chr = line.strip().split()[0][1:]
		else:
			seq_list.append(line.strip())
	dic_chr_size[chr] = len("".join(seq_list))
	del dic_chr_size[""]

	# generate ATSS event per transcript
	dic_iso_chr = {}
	dic_iso_strand = {}
	dic_iso_tss = {}
	dic_iso_tts = {}
	dic_iso_ec = {}
	dic_iso_info = {}
	for line in input_gpd_fl:
		iso_id,chr,strand,tss,tts,cds_s,cds_e,exon_count = line.rstrip("\n").split("\t")[1:9]
		dic_iso_chr[iso_id] = chr
		dic_iso_strand[iso_id] = strand
		dic_iso_tss[iso_id] = int(tss)
		dic_iso_tts[iso_id] = int(tts)
		dic_iso_ec[iso_id] = int(exon_count)
		dic_iso_info[iso_id] = line.rstrip("\n")
		
	ts_count_list = np.random.choice([1,2,3,4,5],len(dic_iso_chr.keys()),p=[0.35,0.25,0.15,0.1,0.15]).tolist()
	i = 0
	for iso in dic_iso_chr.keys():
		if dic_iso_strand[iso] == "+":
			for j in range(0,ts_count_list[i]):
				if dic_iso_tss[iso]-j*int(dis_ts) >= 0:
					info_list = dic_iso_info[iso].split("\t")
					info_list[4] = str(dic_iso_tss[iso]-j*int(dis_ts))
					info_list[-2] = info_list[4] + "," + ",".join(info_list[-2].split(",")[1:])
					print >>output_gpd_fl, "\t".join(info_list)
		else:
			for j in range(0,ts_count_list[i]):
				if dic_iso_tts[iso]+j*int(dis_ts) <= dic_chr_size[dic_iso_chr[iso]]:
					info_list = dic_iso_info[iso].split("\t")
					info_list[5] = str(dic_iso_tts[iso]+j*int(dis_ts))
					if dic_iso_ec[iso] == 1:
						info_list[-1] = info_list[5] + ","
					else:
						info_list[-1] = ",".join(info_list[-1].split(",")[:-2]) + "," + info_list[5] + ","
					print >>output_gpd_fl, "\t".join(info_list)
		i += 1
	genome_fa_fl.close()
	input_gpd_fl.close()
	output_gpd_fl.close()


def do_inputs():
	parser = argparse.ArgumentParser(description="Generate ATS event for each transcript. The possibility matrix {polyA number [1,2,3,4,5] and corresponding possibility [0.35,0.25,0.15,0.1,0.15] } is used.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd file")
	parser.add_argument('-g','--genome',type=argparse.FileType('r'),required=True,help="Genome fasta file")
	parser.add_argument('-d','--distance',type=int,default=50,help="Distance between two adjacent TSSs (bp)")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
