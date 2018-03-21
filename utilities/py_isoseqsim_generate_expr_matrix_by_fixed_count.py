#!/usr/bin/env python
import sys,time,argparse

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	generate_expr_matrix(args.input,args.input_expr,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def generate_expr_matrix(input_gpd_fl,input_txt_fl,output_expr_mtx):
	# parse txt
	dic_iso_expr = {}
	for line in input_txt_fl:
		iso_id,expr_v = line.strip().split("\t")[:2] # expr_v is TPM now
		dic_iso_expr[iso_id] = str(int(round(float(expr_v))))
	
	for line in input_gpd_fl:
		iso_id = line.strip().split("\t")[1]
		print >>output_expr_mtx, line.strip() + "\t" + dic_iso_expr[iso_id]

	input_txt_fl.close()
	input_gpd_fl.close()
	output_expr_mtx.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Randomly generate read count for each isoform based on negative binomial (NB) distribution. Read count is shown in last column of output file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-e','--input_expr',type=argparse.FileType('r'),required=True,help="Input: expression txt file (first colunm is isoform ID, second colunmn is expression value)")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd + read count file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
