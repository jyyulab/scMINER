import sys, os, tarfile, shutil
from contextlib import closing
import numpy as np

import statistics

# Variable definition
total_support = {}
total_edge = []
total_mi = {}
total_pearson = {}
total_spearman = {}
total_regression = {}
total_pvalue = {}
symbol_set = {}
parameters = ''
header = ''
run_num = 0
mu = 0
sigma = 0

# Processing all bootstrap networks, summarizing them into corresponding variables
for adj_file in os.listdir(sys.argv[1]):
	total_edge.append(0)
	with open(sys.argv[1] + adj_file, 'r') as f:	# Opening each bootstrap file
		for line in f:
			if line[0] == '>' and run_num == 0:
				parameters += line
			if line[0] != '>':			# Processing non header lines representing the network 
				tokens = line.split('\t')	# Tokenizing each non header line with tab delimmiter
				hub_id = tokens[0]		# First token is the hub id
				for i in range(1,len(tokens),2):	# Iterating on all adjacent genes: Odd indexes are the connected genes and even indexes are the corresponding value to the edge between the hub gene and the gene with an odd index appearing before the value in the tokens list
					key = hub_id + '--' + tokens[i]	# Creating a key for the edge
					if key in total_support:	
						total_support[key] += 1	# Updating the total number of edges observed for the particualr key (edge)
						total_mi[key] += float(tokens[i+1])	# Updating total MI between the genes involving in the particular key (edge)
					else:
						total_support[key] = 1	# Initializing the total number of edges observed for the particular key (edge) if the key (edge) is newly generated
						total_mi[key] = float(tokens[i+1])	# Initializing total MI between the genes involving in the particulat key (edge) if the key (edge) is newly generated
					total_edge[run_num] += 1	# Increment the total number of edges processed so far
	run_num += 1	# Increment the bootstrap file index

# Network name
path_tokens = sys.argv[1].split('/')
name = path_tokens[len(path_tokens) - 2]
path = '/'.join(path_tokens[0:len(path_tokens) - 2]) + '/'

# Archiving the bootstrap networks
current_path = os.getcwd()
os.chdir(path)
with closing(tarfile.open(name + '.tar.gz', "w:gz")) as tar:
	for f in os.listdir(name):
	    tar.add(name + '/' + f)
os.chdir(current_path)
shutil.rmtree(sys.argv[1])

# Writing out the summary of all bootstrap files into bootstrap_info.txt file
info_file = open(sys.argv[3] + 'bootstrap_info_.txt' if sys.argv[3].endswith('/') else sys.argv[3] + '/bootstrap_info_.txt', 'w')
info_file.write('Total edge tested: ' + str(len(total_support)) + '\n')
info_file.write('Bonferroni corrected (0.05) alpha: ' + str(0.05/len(total_support)) + '\n')

# Computing mu and sigma accross all bootstrap files
for i in range(0, run_num):
	prob = float(total_edge[i]) / float(len(total_support))
	mu += prob
	sigma += prob * (1 - prob)
sigma = np.sqrt(sigma)

info_file.write('mu: ' + str(mu) + '\n')
info_file.write('sigma: ' + str(sigma) + '\n')

# Setting p_threshold if given to the given value, if not to Bonferroni corrected value
p_threshold = 0.05 / len(total_support)
if sys.argv[2] != None:
	p_threshold = float(sys.argv[2])

# Writing out the parameters that the bootstrap networks are constructed with plus other parameters that is used to create consensus network
parameter_file = open(sys.argv[3] + 'parameter_info_.txt' if sys.argv[3].endswith('/') else sys.argv[3] + '/parameter_info_.txt', 'w')
parameters += '>  Bootstrap No: ' + str(run_num) + '\n'
parameters += '>  Source: sjaracne2\n'
parameters += '>  Output network: ' + (sys.argv[3] + 'consensus_network_3col_.txt' if sys.argv[3].endswith('/') else sys.argv[3] + '/consensus_network_3col_.txt') + '\n'
parameter_file.write(parameters)

# Writing out the consensus network preserving edges with statistically significant support
consensus_network = open(sys.argv[3] + 'consensus_network_3col_.txt' if sys.argv[3].endswith('/') else sys.argv[3] + '/consensus_network_3col_.txt', 'w')
header += 'source\ttarget\tMI\n'
consensus_network.write(header)
current_gene = 'none'

for key in sorted(total_support.iterkeys()):	# Iterating on all edges in a sorted fashion
	gene1 = key.split('--')[0]		# Extracting first gene involving in an edge from the key (edge)
	gene2 = key.split('--')[1]		# Extracting second gene involving in an edge from the key (edge)
	z = float(total_support[key] - mu) / float(sigma) if sigma != 0 else 100	# Computing the z score of normal distribution 
	pval = statistics.uprob(z)		# Computing p-value corresponding to the z score --> Implemented in statistics.py module inspired by Statistics::Distributions::uprob function in perl
	if pval < p_threshold:			# Decision making if the edge has enough support or not and therefore if it has to be remained or removed
		mi = '{0:.4f}'.format(float(total_mi[key]) / float(total_support[key]))	# Computing MI corresponding to an edge remaining in the network
		consensus_network.write(gene1 + '\t' + gene2 + '\t' + mi + '\n')

