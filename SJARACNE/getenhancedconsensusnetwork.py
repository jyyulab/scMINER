import os, sys, csv, re, igraph
from scipy import stats
import xlsxwriter
import pandas as pd

marker_set = {}
exp = pd.read_csv(sys.argv[1], sep='\t')
for i in range(exp.shape[0]):
	samples = exp.iloc[i, 2:].values.astype(float)
	marker_set[exp.iloc[i, 0]] = dict(symbol=exp.iloc[i, 1], samples=samples)

gene_network = {}
net = pd.read_csv(sys.argv[2], sep='\t')
for i in range(net.shape[0]):
	gene_network[net.iloc[i, 0] + '--' + net.iloc[i, 1]] = net.iloc[i, 2]

for src_gene in marker_set.keys():
	for tar_gene in marker_set.keys():
		if src_gene != tar_gene:
			if src_gene + '--' + tar_gene in gene_network:
				src_marker = marker_set[src_gene]
				tar_marker = marker_set[tar_gene]
				src_symbol = src_marker['symbol']
				tar_symbol = tar_marker['symbol']
				src_samples = src_marker['samples']
				tar_samples = tar_marker['samples']
				rho, intercept, r, p, stderr = stats.linregress(src_samples, tar_samples)
				scc, sp = stats.spearmanr(src_samples, tar_samples)
				pcc, pp = stats.pearsonr(src_samples, tar_samples)

				mi = gene_network[src_gene + '--' + tar_gene]

				gene_network[src_gene + '--' + tar_gene] = dict(srcSymbol=src_symbol, tarSymbol=tar_symbol, mi=mi, rho=rho, intercept=intercept, r=r, p=p, stderr=stderr, scc=scc, sp=sp, pcc=pcc, pp=pp)

input_net_path_tokens = sys.argv[2].split('/')
input_net_name = input_net_path_tokens[len(input_net_path_tokens)-1]
input_net_name_tokens = re.split('_|\.', input_net_name)
out_name = '_'.join(input_net_name_tokens[0:2]) + '_ncol_' + '_'.join(input_net_name_tokens[3:len(input_net_name_tokens)-1])

out_subnet = None
subnet_list = []
if len(sys.argv) == 5:
	out_subnet = open(sys.argv[3] + 'sub_' + out_name + '.txt', 'w')
	with open(sys.argv[4]) as subnet_file:
		for _id in subnet_file:
			subnet_list.append(_id.split('\n')[0].strip())

out = open(sys.argv[3] + out_name + '.txt', 'w')
out_xlsx = xlsxwriter.Workbook(sys.argv[3] + out_name + '.xlsx')
out_xlsx_sheet = out_xlsx.add_worksheet()
header = ('source','target','source.symbol','target.symbol','MI','pearson','spearman','rho','p-value')
out.write('\t'.join(header) + '\n')
if out_subnet != None:
	out_subnet.write('\t'.join(header) + '\n')
row_index = 0
out_xlsx_sheet.write_row(row_index, 0, header)
row_index += 1
for edge in gene_network.keys():
	info = gene_network[edge]
	tokens = edge.split('--')
	src = tokens[0]
	tar = tokens[1]
	data = (src, tar, info['srcSymbol'], info['tarSymbol'], '{0:.4f}'.format(info['mi']), '{0:.4f}'.format(info['pcc']), '{0:.4f}'.format(info['scc']), '{0:.4f}'.format(info['rho']), '{0:.4f}'.format(info['p']))
	out.write('\t'.join(data) + '\n')
	out_xlsx_sheet.write_row(row_index, 0, data)
	if out_subnet != None:
		if info['srcSymbol'].strip() in subnet_list or info['tarSymbol'].strip() in subnet_list:
			out_subnet.write('\t'.join(data) + '\n')
	row_index += 1

out_xlsx.close()
out.close()
if out_subnet != None:
	out_subnet.close()

graph_name = '_'.join(input_net_name_tokens[0:2]) + '_' + '_'.join(input_net_name_tokens[3:len(input_net_name_tokens)-1])
graph = igraph.Graph.Read_Ncol(sys.argv[2])
graph.write_graphml(sys.argv[3] + graph_name + '.graphml')