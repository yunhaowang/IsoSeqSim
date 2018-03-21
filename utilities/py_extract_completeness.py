#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	extract_end_completeness(args.input_anno,args.input_anno_c,args.input_lr,args.output_cpt5,args.output_cpt3)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def extract_end_completeness(input_anno_fl,input_anno_c_fl,input_lr_fl,output_5cpt_fl,output_3cpt_fl):
	# parse annotation file
	dic_iso_tss = {}
	dic_iso_tts = {}
	for line in input_anno_fl:
		gene,iso,chr,strand,tss,tts = line.strip().split("\t")[:6]
		if strand == "+":
			dic_iso_tss[iso] = int(tss)
			dic_iso_tts[iso] = int(tts)
		else:
			dic_iso_tss[iso] = int(tts)
			dic_iso_tts[iso] = int(tss)
	# parse long read alignment file
	dic_lr_tss = {}
	dic_lr_tts = {}
	for line in input_lr_fl:
		lr_id,lr_id,chr,strand,tss,tts = line.strip().split("\t")[:6]
		if strand == "+":
			dic_lr_tss[lr_id] = int(tss)
			dic_lr_tts[lr_id] = int(tts)
		else:
			dic_lr_tss[lr_id] = int(tts)
			dic_lr_tts[lr_id] = int(tss)
	# parse annotation gpd file with assigned long read
	cpt5_list = []
	cpt3_list = []
	for line in input_anno_c_fl:
		gene,iso,chr,strand,tss,tts,lr_set,lr_count,exon_count,exon_start_set,exon_end_set = line.strip().split("\t")[:11]
		if (lr_set != "") and ("," not in iso) and iso.startswith("refiso_"):
			anno_tss = dic_iso_tss[iso]
			anno_tts = dic_iso_tts[iso]
			if strand == "+":
				for lr in lr_set.split(","):
					if dic_lr_tss[lr] >= anno_tss and dic_lr_tts[lr] <= anno_tts:
						cpt5_list.append(dic_lr_tss[lr]-anno_tss)
						cpt3_list.append(anno_tts-dic_lr_tts[lr])
			else:
				for lr in lr_set.split(","):
					if dic_lr_tss[lr] <= anno_tss and dic_lr_tts[lr] >= anno_tts:
						cpt5_list.append(anno_tss-dic_lr_tss[lr])
						cpt3_list.append(dic_lr_tts[lr]-anno_tts)
	# stat frequency
	cpt5_list.sort()
	cpt3_list.sort()
	dic_cpt5_count = {}
	dic_cpt3_count = {}
	cpt5_total_count = 0
	cpt3_total_count = 0
	cpt5_uniq_list = []
	cpt3_uniq_list = []
	for cpt5 in cpt5_list:
		cpt5_total_count += 1
		if cpt5 not in dic_cpt5_count.keys():
			cpt5_uniq_list.append(cpt5)
			dic_cpt5_count[cpt5] = 1
		else:
			dic_cpt5_count[cpt5] += 1
	for cpt3 in cpt3_list:
		cpt3_total_count += 1
		if cpt3 not in dic_cpt3_count.keys():
			cpt3_uniq_list.append(cpt3)
			dic_cpt3_count[cpt3] = 1
		else:
			dic_cpt3_count[cpt3] += 1

	for cpt5_uniq in cpt5_uniq_list:
		freq = round(dic_cpt5_count[cpt5_uniq]/float(cpt5_total_count),4)
		if cpt5_uniq <= 100 and freq >= 0.0001:
			print >>output_5cpt_fl, str(cpt5_uniq) + "\t" + str(freq)
	for cpt3_uniq in cpt3_uniq_list:
		freq = round(dic_cpt3_count[cpt3_uniq]/float(cpt3_total_count),4)
		if cpt3_uniq <= 100 and freq >= 0.0001:
			print >>output_3cpt_fl, str(cpt3_uniq) + "\t" + str(freq)

	input_anno_fl.close()
	input_anno_c_fl.close()
	input_lr_fl.close()
	output_5cpt_fl.close()
	output_3cpt_fl.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Extract the completeness information at 5' and 3'end of transcripts. Output deleted nucleotide number and its frequency. Note: only output when the deleted nucleotide number is <= 100 and the frequency is >= 0.0001", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-a','--input_anno',type=argparse.FileType('r'),required=True,help="Input: annotation file, gpd format")
	parser.add_argument('-c','--input_anno_c',type=argparse.FileType('r'),required=True,help="Input: annotation file with assigned long read ID(column 7th), gpd format")
	parser.add_argument('-l','--input_lr',type=argparse.FileType('r'),required=True,help="Input: aligned long read file, gpd format")
	parser.add_argument('-5','--output_cpt5',type=argparse.FileType('w'),required=True,help="Output: 5'end completeness stat file")
	parser.add_argument('-3','--output_cpt3',type=argparse.FileType('w'),required=True,help="Output: 3'end completeness stat file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
