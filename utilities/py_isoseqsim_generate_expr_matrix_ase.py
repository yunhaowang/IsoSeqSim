#!/usr/bin/env python
import sys,time,argparse
import numpy as np
import scipy.stats as sp

def main(args):
	sys.stdout.write("Start analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()
	generate_expr_matrix(args.input,args.nb_n,args.nb_p,args.tn_lower,args.tn_upper,args.tn_mu,args.tn_sigma,args.output)
	sys.stdout.write("Finish analysis: " + time.strftime("%a,%d %b %Y %H:%M:%S") + "\n")
	sys.stdout.flush()

def generate_expr_matrix(input_gpd_fl,nb_r,nb_p,tn_lower,tn_upper,tn_mu,tn_sigma,output_expr_mtx):
	# parse gpd file
	gpd_list = []
	for line in input_gpd_fl:
		gpd_list.append(line.strip())
	# generate random read count based on negative binomial distribution
	nb_list = np.random.negative_binomial(nb_r,nb_p,len(gpd_list)).tolist()
	# generate ASE distribution based on truncated normal distribution
	tn_list = sp.truncnorm.rvs((tn_lower-tn_mu)/tn_sigma,(tn_upper-tn_mu)/tn_sigma,loc=tn_mu,scale=tn_sigma,size=len(gpd_list)).tolist()
	i = 0
	for gpd in gpd_list:
		total_read_count = nb_list[i]
		allele1_read_count = int(round(total_read_count*tn_list[i]))
		allele2_read_count = total_read_count - allele1_read_count
		print >>output_expr_mtx, gpd + "\t" + str(nb_list[i]) + "\t" + str(allele1_read_count) + "\t" + str(allele2_read_count)
		i += 1
	input_gpd_fl.close()
	output_expr_mtx.close()

def do_inputs():
	parser = argparse.ArgumentParser(description="First, generate total read count for each isoform based on negative binomial distribution; then, generate allele-specific read count based on truncated normal distribution. Output format = gpd + total read count + allele1 read count + allele2 read count.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument('-i','--input',type=argparse.FileType('r'),required=True,help="Input: gpd file")
	parser.add_argument('-n','--nb_n',type=int,default=10,help="Parameter of the Negative Binomial distribution, n")
	parser.add_argument('-p','--nb_p',type=float,default=0.5,help="Parameter of the Negative Binomial distribution, p")
	parser.add_argument('-l','--tn_lower',type=float,default=0.0,help="Parameter of the Truncated Normal distribution, lower boundary")
	parser.add_argument('-u','--tn_upper',type=float,default=1.0,help="Parameter of the Truncated Normal distribution, upper boundary")
	parser.add_argument('-m','--tn_mu',type=float,default=0.5,help="Parameter of the Truncated Normal distribution, mean value")
	parser.add_argument('-s','--tn_sigma',type=float,default=0.1,help="Parameter of the Truncated Normal distribution, stdandard variation")
	parser.add_argument('-o','--output',type=argparse.FileType('w'),required=True,help="Output: gpd + total read count + allele1 read count + allele 2 read count file")
	args = parser.parse_args()
	return args

if __name__=="__main__":
	args = do_inputs()
	main(args)
