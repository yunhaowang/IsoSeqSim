#!/usr/bin/env python
import sys,time,random,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	gpd_list = generate_fusion_event(args.input)
	ouput_fusion_transcript(gpd_list,args.number,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def generate_fusion_event(input_gpd_fl):
	gpd_list = []
	#--- read gpd file ---
	dic_iso_info = {}
	dic_iso_ec = {}
	dic_gene_iso = {}
	dic_gene_ec = {}
	for line in input_gpd_fl:
		gene,iso,chr,strand,tss,tts,cds_s,cds_t,exon_count,exon_start_set,exon_end_set = line.rstrip("\n").split("\t")
		dic_iso_info[iso] = line.rstrip("\n")
		dic_iso_ec[iso] = int(exon_count)
		if gene not in dic_gene_iso.keys():
			dic_gene_iso[gene] = []
			dic_gene_ec[gene] = []
			dic_gene_iso[gene].append(iso)
			dic_gene_ec[gene].append(int(exon_count))
		else:
			dic_gene_iso[gene].append(iso)
			dic_gene_ec[gene].append(int(exon_count))
	#--- choose one isoform with most exons per gene ---
	iso_cdd_set = set()
	for gene in dic_gene_iso.keys():
		max_ec = max(dic_gene_ec[gene])
		for iso in dic_gene_iso[gene]:
			if dic_iso_ec[iso] == max_ec:
				iso_cdd_set.add(iso)
				break
	#--- generate fusion transcript ---
	#--- total number of generated fusion transcripts = (# of gene)*(# of gene - 1) ---
	for iso in iso_cdd_set:
		new_iso_cdd_list = list(iso_cdd_set)
		new_iso_cdd_list.remove(iso)
		info_list = dic_iso_info[iso].split("\t")
		if info_list[3] == "+":
			if int(info_list[8]) == 1:
				info_list[5] = str(int(info_list[4])+(int(info_list[5])-int(info_list[4]))/2)
				info_list[-1] = info_list[5] + ","
			else:
				exon_first = int(info_list[8])/2
				info_list[5] = info_list[-1].split(",")[exon_first-1]
				info_list[-3] = str(exon_first)
				info_list[-2] = ",".join(info_list[-2].split(",")[:exon_first]) + ","
				info_list[-1] = ",".join(info_list[-1].split(",")[:exon_first]) + ","
		else:
			if int(info_list[8]) == 1:
				info_list[4] = str(int(info_list[4])+(int(info_list[5])-int(info_list[4]))/2)
				info_list[-2] = info_list[4] + ","
			else:
				exon_first = int(info_list[8])/2
				info_list[4] = info_list[-2].split(",")[exon_first]
				info_list[-3] = str(exon_first)
				info_list[-2] = ",".join(info_list[-2].split(",")[-(exon_first+1):])
				info_list[-1] = ",".join(info_list[-1].split(",")[-(exon_first+1):])
		for new_iso in new_iso_cdd_list:
			new_info_list = dic_iso_info[new_iso].split("\t")
			if new_info_list[3] == "+":
				if int(new_info_list[8]) == 1:
					new_info_list[4] = str(int(new_info_list[4])+(int(new_info_list[5])-int(new_info_list[4]))/2)
					new_info_list[-2] = new_info_list[4] + ","
				else:
					exon_first = int(new_info_list[8])/2
					new_info_list[4] = new_info_list[-2].split(",")[exon_first]
					new_info_list[-3] = str(exon_first)
					new_info_list[-2] = ",".join(new_info_list[-2].split(",")[-(exon_first+1):])
					new_info_list[-1] = ",".join(new_info_list[-1].split(",")[-(exon_first+1):])
			else:
				if int(new_info_list[8]) == 1:
					new_info_list[5] = str(int(new_info_list[4])+(int(new_info_list[5])-int(new_info_list[4]))/2)
					new_info_list[-1] = new_info_list[5] + ","
				else:
					exon_first = int(new_info_list[8])/2
					new_info_list[5] = new_info_list[-1].split(",")[exon_first-1]
					new_info_list[-3] = str(exon_first)
					new_info_list[-2] = ",".join(new_info_list[-2].split(",")[:exon_first]) + ","
					new_info_list[-1] = ",".join(new_info_list[-1].split(",")[:exon_first]) + ","
	
			gpd_list.append("\t".join(info_list) + "\t" + "\t".join(new_info_list))

	input_gpd_fl.close()
	return gpd_list

def ouput_fusion_transcript(gpd_list,fusion_count,output_gpd_fl):
	if fusion_count == "all":
		for gpd in gpd_list:
			print >>output_gpd_fl, gpd
	else:
		if int(fusion_count) >= len(gpd_list):
			sys.stdout.write("The total count of genereated fusion transcripts is " + str(len(gpd_list)) + ". Even though you set a big number, we only output all fusion transcript for you!!!\n")
			for gpd in gpd_list:
				print >>output_gpd_fl, gpd
		else:
			new_gpd_list = random.sample(gpd_list,int(fusion_count))
			for gpd in new_gpd_list:
				print >>output_gpd_fl, gpd
	output_gpd_fl.close()

def do_inputs():
	output_gpd_format = '''
1-11 columns are gpd format for first(left or 5'end) part of fusion transcript	
1. gene id
2. isoform id
3. chromosome id
4. strand
5. TSS
6. TTS
7. .
8. .
9. exon count
10. exon start set
11. exon end set

12-22 columns are gpd format for second(right or 3'end) part of fusion transcript'''
	parser = argparse.ArgumentParser(description="Generate fusion transcripts using gpd file. One isoform with most exon count per gene is used, and total number of generated fusion transcripts is equal to (total gene count)*(total gene count - 1)", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd file")
	parser.add_argument('-n','--number',type=str,default="all",help="Number of fusion transcript to generate.")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
