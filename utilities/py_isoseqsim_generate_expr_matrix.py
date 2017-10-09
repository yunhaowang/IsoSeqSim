#!/usr/bin/env python
import sys,time,argparse
import numpy as np

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	generate_expr_matrix(args.input,args.nb_n,args.nb_p,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def generate_expr_matrix(input_gpd_fl,nb_r,nb_p,output_expr_mtx):
	# parse gpd file
	gpd_list = []
	for line in input_gpd_fl:
		gpd_list.append(line.strip())
	# generate random read count based on negative binomial distribution
	nb_list = np.random.negative_binomial(nb_r,nb_p,len(gpd_list)).tolist()
	i = 0
	for gpd in gpd_list:
		print >>output_expr_mtx, gpd + "\t" + str(nb_list[i])
		i += 1
	input_gpd_fl.close()
	output_expr_mtx.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="Randomly generate read count for each isoform based on negative binomial (NB) distribution. Read count is shown in last column of output file", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-n','--nb_n',type=int,default=10,help="Parameter of the NB distribution, n")
	parser.add_argument('-p','--nb_p',type=float,default=0.5,help="Parameter of the NB distribution, p")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd + read count file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
