import os, sys

scMINER_PATH = os.environ['scMINER_PATH']
PYTHON_PATH = os.environ['PYTHON_PATH']
sys.path.insert(0, scMINER_PATH)
sys.path.insert(0, PYTHON_PATH)

import utils
import argparse
import numpy as np
import pandas as pd
import shlex, subprocess, shutil

fig_num = 1

def setup(args):
	parser = argparse.ArgumentParser(description='MICA package.')
	parser.add_argument('mode', help='Procedure to select and execute from the package [Clust]')
	parser.add_argument('project_name', help='Project name')
	args_ = parser.parse_args(args[1:3])
	if args_.mode == 'Clust':
		parser.add_argument('input_file', help='Input file in hdf format')
		parser.add_argument('dist_file', help='Distance square matrix file name in hdf format')
		parser.add_argument('--min_dim', type=int, default=19, help='Minimum number of dimensions used in clustering (default: 19)')
		parser.add_argument('--max_dim', type=int, default=19, help='Maximum number of dimensions used in clustering (default: 19)')
		parser.add_argument('--bootstrap', type=int, default=10, help='Maximum number of iterations per dimnesion (default: 10)')
		parser.add_argument('--k', type=int, default=[2], nargs='+', help='Number of clusters to divide the dataset to (default: [2])')
		parser.add_argument('--transformation', default='MDS', help='Transformation method used for dimension reduction [MDS | PCA | LPL | LPCA] (default: MDS)')
		parser.add_argument('--hclust', default='False', help='Whether apply hierarchical clustering or not (default: False')
		parser.add_argument('--retransformation', default=150, help='If None or False no retransformation is used in the final visualization, if a number, retransformation threshold will be set to the value (default: 150)')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, nargs='+', default=[2000]*6, help='Memory allocation for each individual step in clustering pipeline (default: 2GB)')
		parser.add_argument('--queue', default='compbio', help='Queue name for job allocation')
		parser.add_argument('--preclust', default='None', help='Whether the clustering is based on previous clustering result or not (default: None)')
		parser.add_argument('--perplexity', default=30, help='Visualization parameter determining how dense the clusters are')
	if args_.mode == 'Reclust':
		parser.add_argument('--transformation', default='MDS', help='Transformation method used for dimension reduction [MDS | PCA | LPL | LPCA] (default: MDS)')
		parser.add_argument('--retransformation', nargs='+', default=[60, 80, 100, 120, 150], help='Retransformation threshold will be set to the value (default: [60, 80, 100, 120, 150])')
		parser.add_argument('--max_dim', type=int, default=19, help='Maximum number of dimensions used in clustering (default: 19)')
		parser.add_argument('--k', type=int, default=[2], nargs='+', help='Number of clusters to divide the dataset to (default: [2])')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, nargs='+', default=[2000]*2, help='Memory allocation for each individual step in clustering pipeline (default: 2GB)')
		parser.add_argument('--queue', default='compbio', help='Queue name for job allocation')
	args_ = parser.parse_args(args[1:])
	return args_

def setup_directory(out_dir, project_name, method, ks):
	scMINER_path = []
	log_path = []
	out_path = []
	script_path = []
	for k in ks:
		scMINER_path.append(out_dir + 'scMINER_' + project_name + '_' + method + '_' + str(k) + '/')
		log_path.append(scMINER_path[-1] + 'scMINER_MICA_log/')
		out_path.append(scMINER_path[-1] + 'scMINER_MICA_out/')
		script_path.append(scMINER_path[-1] + 'scMINER_MICA_scripts/')
		if not os.path.exists(scMINER_path[-1]):
			os.mkdir(scMINER_path[-1])
		if not os.path.exists(log_path[-1]):
			os.mkdir(log_path[-1])
		if not os.path.exists(out_path[-1]):
			os.mkdir(out_path[-1])
		if not os.path.exists(script_path[-1]):
			os.mkdir(script_path[-1])
	return [scMINER_path, log_path, out_path, script_path]

def prep(args, paths):
	for path in paths[2]:
		path_tmp = path + '.tmp/'
		if not os.path.exists(path_tmp):
			os.mkdir(path_tmp)
		shutil.copyfile(args.dist_file, path_tmp + args.dist_file.split('/')[-1])
		shutil.copyfile(args.input_file, path_tmp + args.input_file.split('/')[-1])

def hclust(args, paths):
	global fig_num
	for i in range(len(paths[2])):
		path_tmp = paths[2][i] + '.tmp/'
		if not os.path.exists(path_tmp):
			exit()
		script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Hclust.py ' + ' ' + str(fig_num) + ' ' + path_tmp + args.dist_file.split('/')[-1] + ' ' + paths[2][i] + ' ' + args.outfilename + ' '
		out_1 = open(paths[3][i] + '01_Hclust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		out_1.write(script + '\n')
		out_1.close()
		fig_num += 1

def transform(args, paths):
	global fig_num
	for i in range(len(paths[2])):
		path_tmp = paths[2][i] + '.tmp/'
		if not os.path.exists(path_tmp):
			exit()
		script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Transform.py ' + ' ' + args.transformation + ' ' + str(fig_num) + ' ' + path_tmp + args.dist_file.split('/')[-1] + ' ' + str(args.max_dim) + ' ' + paths[2][i] + ' ' + args.outfilename + ' ' + str(args.perplexity) + ' ' + args.preclust + ' '
		out_2 = open(paths[3][i] + '02_Transform_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		out_2.write(script + '\n')
		out_2.close()
		fig_num += 1

def kmeans(args, paths):
	for i in range(len(paths[2])):
		path_tmp = paths[2][i] + '.tmp/'
		if not os.path.exists(path_tmp):
			os.mkdir(path_tmp)
		out_3 = open(paths[3][i] + '03_Prep_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		out_4 = open(paths[3][i] + '04_Kmeans_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		iterations = args.bootstrap*(args.max_dim - args.min_dim + 1)
		for j in range(iterations):
			script = 'cp ' + paths[2][i] + args.outfilename + '_clust.h5 ' + path_tmp + args.outfilename + '_clust.h5.tmp.' + str(j) + ' '
			out_3.write(script + '\n')
			if args.transformation == 'LPCA':
				script = 'cp ' + paths[2][i] + args.outfilename + '_clust.h5 ' + path_tmp + args.outfilename + '_clust.h5.tmp.' + str(iterations + j) + ' '
				out_3.write(script + '\n')
		for j in range(args.bootstrap):
			dim = args.max_dim - args.min_dim + 1
			for k in range(dim):
				if args.transformation == 'MDS':
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.outfilename + '_clust.h5.tmp.' + str(j * dim + k) + ' mds ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.outfilename + '_kmeans.h5.tmp.' + str(j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
				elif args.transformation == 'PCA':
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.outfilename + '_clust.h5.tmp.' + str(j * dim + k) + ' pca ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.outfilename + '_kmeans.h5.tmp.' + str(j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
				elif args.transformation == 'LPL':
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.outfilename + '_clust.h5.tmp.' + str(j * dim + k) + ' lpl ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.outfilename + '_kmeans.h5.tmp.' + str(j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
				elif args.transformation == 'LPCA':
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.outfilename + '_clust.h5.tmp.' + str(j * dim + k) + ' pca ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.outfilename + '_kmeans.h5.tmp.' + str(j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.outfilename + '_clust.h5.tmp.' + str(iterations + j * dim + k) + ' lpl ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.outfilename + '_kmeans.h5.tmp.' + str(iterations + j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
		out_3.close()
		out_4.close()

def cclust(args, paths):
	global fig_num
	for i in range(len(paths[2])):
		path_tmp = paths[2][i] + '.tmp/'
		if not os.path.exists(path_tmp):
			exit()
		script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Cclust.py ' + ' ' + str(fig_num) + ' ' + paths[2][i] + args.outfilename + '_clust.h5 ' + path_tmp + ' ' + str(args.k[i]) + ' ' + args.project_name + ' kmeans.h5.tmp ' + args.transformation.lower() + ' ' + paths[2][i] + ' ' + args.outfilename + ' ' + str(args.max_dim) + ' ' + str(args.retransformation) + ' ' 
		out_5 = open(paths[3][i] + '05_CClust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		out_5.write(script + '\n')
		out_5.close()
		fig_num += 3

def ggplot(args, paths):
	for i in range(len(paths[2])):
		script = 'Rscript ' + scMINER_PATH + 'MICA/ggplot.cc.r ' + ' ' + paths[2][i] + args.outfilename + '.ggplot.txt 1 5 ' + args.outfilename + ' ' + paths[2][i] + args.outfilename + '_clust_k' + str(args.k[i]) + '.rplot.pdf '
		out_6 = open(paths[3][i] + '06_GGplot_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		out_6.write(script + '\n')
		out_6.close()

def clust(args, paths):
	for i in range(len(paths[3])):
		out_0 = open(paths[3][i] + '00_Clust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		if args.host == 'LSF':
			if args.hclust == 'True':
				script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_HClust -q ' + args.queue + ' -M ' + str(args.resource[0]) + ' -i ' + paths[3][i] + '01_Hclust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh -oo ' + paths[1][i] + args.project_name + '_MICA_Hclust.%J.%I.out -eo ' + paths[1][i] + args.project_name + '_MICA_Hclust.%J.%I.err \n'
				out_0.write(script)
			script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Transform -q ' + args.queue + ' -M ' + str(args.resource[1]) + ' -i ' + paths[3][i] + '02_Transform_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh -oo ' + paths[1][i] + args.project_name + '_MICA_Transform.%J.%I.out -eo ' + paths[1][i] + args.project_name + '_MICA_Transform.%J.%I.err \n'
			out_0.write(script)
			script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Kmeans.0 -q ' + args.queue + ' -M ' + str(args.resource[2]) + ' -i ' + paths[3][i] + '03_Prep_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh -oo ' + paths[1][i] + args.project_name + '_MICA_Kmeans.0.%J.%I.out -eo ' + paths[1][i] + args.project_name + '_MICA_Kmeans.0.%J.%I.err \n'
			out_0.write(script)
			script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Kmeans.1 -q ' + args.queue + ' -M ' + str(args.resource[3]) + ' -i ' + paths[3][i] + '04_Kmeans_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh -oo ' + paths[1][i] + args.project_name + '_MICA_Kmeans.1.%J.%I.out -eo ' + paths[1][i] + args.project_name + '_MICA_Kmeans.1.%J.%I.err \n'
			out_0.write(script)
			script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_CClust -q ' + args.queue + ' -M ' + str(args.resource[4]) + ' -i ' + paths[3][i] + '05_CClust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh -oo ' + paths[1][i] + args.project_name + '_MICA_Cclust.%J.%I.out -eo ' + paths[1][i] + args.project_name + '_MICA_Cclust.%J.%I.err \n'
			out_0.write(script)
			script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_GGplot -q ' + args.queue + ' -M ' + str(args.resource[5]) + ' -i ' + paths[3][i] + '06_GGplot_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh -oo ' + paths[1][i] + args.project_name + '_MICA_Ggplot.%J.%I.out -eo ' + paths[1][i] + args.project_name + '_MICA_Ggplot.%J.%I.err \n'
			out_0.write(script)
		elif args.host == 'LOCAL':
			if args.hclust == 'True':
				script = 'sh ' + paths[3][i] + '01_Hclust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh >> ' + paths[1][i] + args.project_name + '_MICA_Hclust.out \n'
				out_0.write(script)
			script = 'sh ' + paths[3][i] + '02_Transform_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh >> ' + paths[1][i] + args.project_name + '_MICA_Transform.out \n'
			out_0.write(script)
			script = 'sh ' + paths[3][i] + '03_Prep_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh >> ' + paths[1][i] + args.project_name + '_MICA_Kmeans.0.out \n'
			out_0.write(script)
			script = 'sh ' + paths[3][i] + '04_Kmeans_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh >> ' + paths[1][i] + args.project_name + '_MICA_Kmeans.1.out \n'
			out_0.write(script)
			out_0.write('jobs=$(ps -ef | grep \"' + args.project_name + '_' + args.transformation + '_' + str(args.k[i]) + '\" | grep Kmeans.py -c)\n')
			out_0.write('while [ $jobs -gt 0 ]\ndo\n\tsleep 30\n\tjobs=$(ps -ef | grep \"' + args.project_name + '_' + args.transformation + '_' + str(args.k[i]) + '\" | grep Kmeans.py -c)\ndone\n')
			script = 'sh ' + paths[3][i] + '05_CClust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh >> ' + paths[1][i] + args.project_name + '_MICA_Cclust.out \n'
			out_0.write(script)			
			script = 'sh ' + paths[3][i] + '06_GGplot_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh >> ' + paths[1][i] + args.project_name + '_MICA_Ggplot.out \n'
			out_0.write(script)
			out_0.write('jobs=$(ps -ef | grep \"' + args.project_name + '_' + args.transformation + '_' + str(args.k[i]) + '\" | grep ggplot.cc.r -c)\n')
			out_0.write('while [ $jobs -gt 0 ]\ndo\n\tsleep 30\n\tjobs=$(ps -ef | grep \"' + args.project_name + '_' + args.transformation + '_' + str(args.k[i]) + '\" | grep ggplot.cc.r -c)\ndone\n')
		out_0.write('rm -rf ' + paths[2][i] + '.tmp/ \n')
		out_0.close()

def reduce(args, paths):
	path_tmp = args.outdir + '.tmp/'
	if not os.path.exists(path_tmp):
		os.mkdir(path_tmp)
	out_1 = open(path_tmp + '01_Hclust_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	out_2 = open(path_tmp + '02_Transform_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	out_3 = open(path_tmp + '03_Share_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	out_4 = open(path_tmp + '04_Prep_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	out_5 = open(path_tmp + '05_Kmeans_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	out_6 = open(path_tmp + '06_CClust_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	out_7 = open(path_tmp + '07_GGplot_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	out_8 = open(path_tmp + '08_Clean_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	for i in range(len(paths[3])):
		with open(paths[3][i] + '01_Hclust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'r') as f1:
			out_1.write(f1.read())
		if i == 0:
			with open(paths[3][i] + '02_Transform_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'r') as f2:
				out_2.write(f2.read())
		else:
			script = 'cp ' + paths[2][0] + args.outfilename + '_clust.h5 ' + paths[2][i] + '\n'
			out_3.write(script)
			if args.transformation != 'LPCA':
				script = 'cp ' + paths[2][0] + args.outfilename + '_' + args.transformation.lower() + '.pdf ' + paths[2][i] + '\n'
				out_3.write(script)
			else:
				script = 'cp ' + paths[2][0] + args.outfilename + '_pca.pdf ' + paths[2][i] + '\n'
				out_3.write(script)
				script = 'cp ' + paths[2][0] + args.outfilename + '_lpl.pdf ' + paths[2][i] + '\n'
				out_3.write(script)
		with open(paths[3][i] + '03_Prep_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'r') as f3:
			out_4.write(f3.read())
		with open(paths[3][i] + '04_Kmeans_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'r') as f4:
			out_5.write(f4.read())
		with open(paths[3][i] + '05_CClust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'r') as f5:
			out_6.write(f5.read())
		with open(paths[3][i] + '06_GGplot_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'r') as f6:
			out_7.write(f6.read())
		out_8.write('rm -rf ' + paths[2][i] + '.tmp/ \n')	
	out_1.close()
	out_2.close()
	out_3.close()
	out_4.close()
	out_5.close()
	out_6.close()
	out_7.close()
	out_8.close()

def reduce_clust(args, paths):
	path_tmp = args.outdir + '.tmp/'
	path_tmp_log = path_tmp + '.log/'
	if not os.path.exists(path_tmp):
		exit()
	if not os.path.exists(path_tmp_log):
		os.mkdir(path_tmp_log)
	out_0 = open(path_tmp + '00_Clust_' + args.project_name + '_' + args.transformation.lower() + '.sh', 'w')
	if args.host == 'LSF':
		if args.hclust == 'True':
			script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_HClust -q ' + args.queue + ' -M ' + str(args.resource[0]) + ' -i ' + path_tmp + '01_Hclust_' + args.project_name + '_' + args.transformation.lower() + '.sh -oo ' + path_tmp_log + args.project_name + '_MICA_Hclust.%J.%I.out -eo ' + path_tmp_log + args.project_name + '_MICA_Hclust.%J.%I.err \n'
			out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Transform.0 -q ' + args.queue + ' -M ' + str(args.resource[1]) + ' -i ' + path_tmp + '02_Transform_' + args.project_name + '_' + args.transformation.lower() + '.sh -oo ' + path_tmp_log + args.project_name + '_MICA_Transform.0.%J.%I.out -eo ' + path_tmp_log + args.project_name + '_MICA_Transform.0.%J.%I.err \n'
		out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Transform.1 -q ' + args.queue + ' -M ' + str(2000) + ' -i ' + path_tmp + '03_Share_' + args.project_name + '_' + args.transformation.lower() + '.sh -oo ' + path_tmp_log + args.project_name + '_MICA_Transform.1.%J.%I.out -eo ' + path_tmp_log + args.project_name + '_MICA_Transform.0.%J.%I.err \n'
		out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Kmeans.0 -q ' + args.queue + ' -M ' + str(args.resource[2]) + ' -i ' + path_tmp + '04_Prep_' + args.project_name + '_' + args.transformation.lower() + '.sh -oo ' + path_tmp_log + args.project_name + '_MICA_Kmeans.0.%J.%I.out -eo ' + path_tmp_log + args.project_name + '_MICA_Kmeans.0.%J.%I.err \n'
		out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Kmeans.1 -q ' + args.queue + ' -M ' + str(args.resource[3]) + ' -i ' + path_tmp + '05_Kmeans_' + args.project_name + '_' + args.transformation.lower() + '.sh -oo ' + path_tmp_log + args.project_name + '_MICA_Kmeans.1.%J.%I.out -eo ' + path_tmp_log + args.project_name + '_MICA_Kmeans.1.%J.%I.err \n'
		out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_CClust -q ' + args.queue + ' -M ' + str(args.resource[4]) + ' -i ' + path_tmp + '06_CClust_' + args.project_name + '_' + args.transformation.lower() + '.sh -oo ' + path_tmp_log + args.project_name + '_MICA_Cclust.%J.%I.out -eo ' + path_tmp_log + args.project_name + '_MICA_Cclust.%J.%I.err \n'
		out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_GGplot -q ' + args.queue + ' -M ' + str(args.resource[5]) + ' -i ' + path_tmp + '07_GGplot_' + args.project_name + '_' + args.transformation.lower() + '.sh -oo ' + path_tmp_log + args.project_name + '_MICA_Ggplot.%J.%I.out -eo ' + path_tmp_log + args.project_name + '_MICA_Ggplot.%J.%I.err \n'
		out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Clean -q ' + args.queue + ' -M ' + str(args.resource[5]) + ' -i ' + path_tmp + '08_Clean_' + args.project_name + '_' + args.transformation.lower() + '.sh -oo ' + path_tmp_log + args.project_name + '_MICA_Clean.%J.%I.out -eo ' + path_tmp_log + args.project_name + '_MICA_Clean.%J.%I.err \n'
		out_0.write(script)
	elif args.host == 'LOCAL':
		if args.hclust == 'True':
			script = 'sh ' + path_tmp + '01_Hclust_' + args.project_name + '_' + args.transformation.lower() + '.sh >> ' + path_tmp_log + args.project_name + '_MICA_Hclust.out \n'
			out_0.write(script)
		script = 'sh ' + path_tmp + '02_Transform_' + args.project_name + '_' + args.transformation.lower() + '.sh >> ' + path_tmp_log + args.project_name + '_MICA_Transform.0.out \n'
		out_0.write(script)
		script = 'sh ' + path_tmp + '03_Share_' + args.project_name + '_' + args.transformation.lower() + '.sh >> ' + path_tmp_log + args.project_name + '_MICA_Transform.1.out \n'
		out_0.write(script)
		script = 'sh ' + path_tmp + '04_Prep_' + args.project_name + '_' + args.transformation.lower() + '.sh >> ' + path_tmp_log + args.project_name + '_MICA_Kmeans.0.out \n'
		out_0.write(script)
		script = 'sh ' + path_tmp + '05_Kmeans_' + args.project_name + '_' + args.transformation.lower() + '.sh >> ' + path_tmp_log + args.project_name + '_MICA_Kmeans.1.out \n'
		out_0.write(script)
		out_0.write('jobs=$(ps -ef | grep \"' + args.project_name + '_' + args.transformation + '\" | grep Kmeans.py -c)\n')
		out_0.write('while [ $jobs -gt 0 ]\ndo\n\tsleep 30\n\tjobs=$(ps -ef | grep \"' + args.project_name + '_' + args.transformation + '\" | grep Kmeans.py -c)\ndone\n')
		script = 'sh ' + path_tmp + '06_CClust_' + args.project_name + '_' + args.transformation.lower() + '.sh >> ' + path_tmp_log + args.project_name + '_MICA_Cclust.out \n'
		out_0.write(script)			
		script = 'sh ' + path_tmp + '07_GGplot_' + args.project_name + '_' + args.transformation.lower() + '.sh >> ' + path_tmp_log + args.project_name + '_MICA_Ggplot.out \n'
		out_0.write(script)
		script = 'sh ' + path_tmp + '08_Clean_' + args.project_name + '_' + args.transformation.lower() + '.sh >> ' + path_tmp_log + args.project_name + '_MICA_Clean.out \n'
		out_0.write(script)
	out_0.close()

def retransform(args, paths):
	global fig_num
	for i in range(len(paths[2])):
		if not os.path.exists(paths[2][i] + '.tmp/'):
			os.mkdir(paths[2][i] + '.tmp/')
		out_7 = open(paths[3][i] + '01_Retransform_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		for j in range(len(args.retransformation)):
			script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Retransform.py ' + ' ' + str(fig_num) + ' ' + paths[2][i] + args.outfilename + '_clust.h5 ' + args.transformation.lower() + ' ' + str(args.max_dim) + ' ' + paths[2][i] + ' ' + args.outfilename + ' ' + str(args.retransformation[j]) + ' '
			out_7.write(script + '\n')
			fig_num += 1
		out_7.close()

def reggplot(args, paths):
	for i in range(len(paths[2])):
		out_8 = open(paths[3][i] + '02_Reggplot_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		for j in range(len(args.retransformation)):
			script = 'Rscript ' + scMINER_PATH + 'MICA/ggplot.cc.r ' + ' ' + paths[2][i] + args.outfilename + '_reclust_re' + str(args.retransformation[j]) + '.ggplot.txt 1 5 ' + args.outfilename + ' ' + paths[2][i] + args.outfilename + '_re' + str(args.retransformation[j]) + '_clust_k' + str(args.k[i]) + '.rplot.pdf '
			out_8.write(script + '\n')
		out_8.close()

def reclust(args, paths):
	for i in range(len(paths[3])):
		out_0 = open(paths[3][i] + '00_Reclust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		if args.host == 'LSF':
			script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Reclust -q ' + args.queue + ' -M ' + str(args.resource[0]) + ' -i ' + paths[3][i] + '01_Retransform_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh -oo ' + paths[1][i] + args.project_name + '_MICA_Reclust.%J.%I.out -eo ' + paths[1][i] + args.project_name + '_MICA_Reclust.%J.%I.err \n'
			out_0.write(script)
			script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MICA_Reggplot -q ' + args.queue + ' -M ' + str(args.resource[1]) + ' -i ' + paths[3][i] + '02_Reggplot_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh -oo ' + paths[1][i] + args.project_name + '_MICA_Reggplot.%J.%I.out -eo ' + paths[1][i] + args.project_name + '_MICA_Reggplot.%J.%I.err \n'
			out_0.write(script)
			out_0.write('rm -rf ' + paths[2][i] + '.tmp/ \n')
		elif args.host == 'LOCAL':
			script = 'sh ' + paths[3][i] + '01_Retransform_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh >> ' + paths[1][i] + args.project_name + '_MICA_Retransform.out \n'
			out_0.write(script)
			script = 'sh ' + paths[3][i] + '02_Reggplot_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh >> ' + paths[1][i] + args.project_name + '_MICA_Reggplot.out \n'
			out_0.write(script)
			out_0.write('rm -rf ' + paths[2][i] + '.tmp/ \n')
		out_0.close()

def run(args):
	args_ = setup(args)
	if args_.mode == 'Clust':
		paths = setup_directory(args_.outdir, args_.project_name, args_.transformation, args_.k)
		if args_.preclust != 'None':
			args_.hclust = 'False'
		prep(args_, paths)		
		hclust(args_, paths)
		transform(args_, paths)
		kmeans(args_, paths)
		cclust(args_, paths)
		ggplot(args_, paths)
		clust(args_, paths)
		if len(args_.k) > 1:
			reduce(args_, paths)
			reduce_clust(args_, paths)
			path_tmp = args_.outdir + '.tmp/'
			path_tmp_log = path_tmp + '.log/'
			if args_.host == 'LSF':
				script = 'bsub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MICA_Clust -q ' + args_.queue + ' -R \"rusage[mem=2000]\" -oo ' + path_tmp_log + args_.project_name + '_MICA_Clust.out -eo ' + path_tmp_log + args_.project_name + '_MICA_Clust.err sh ' + path_tmp + '00_Clust_' + args_.project_name + '_' + args_.transformation.lower() + '.sh \n'
				subprocess.Popen(shlex.split(script))
			elif args_.host == 'LOCAL':
				script = 'sh ' + path_tmp + '00_Clust_' + args_.project_name + '_' + args_.transformation.lower() + '.sh >> ' + path_tmp_log + args_.project_name + '_MICA_Clust.out \n'
				subprocess.Popen(shlex.split(script))
		else:
			""" Old method when redundant transformation ran for multiple k on an input file """
			for i in range(len(paths[3])):
				if args_.host == 'LSF':
					script = 'bsub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MICA_Clust -q ' + args_.queue + ' -R \"rusage[mem=2000]\" -oo ' + paths[1][i] + args_.project_name + '_MICA_Clust.out -eo ' + paths[1][i] + args_.project_name + '_MICA_Clust.err sh ' + paths[3][i] + '00_Clust_' + args_.project_name + '_' + args_.transformation.lower() + '_' + str(args_.k[i]) + '.sh \n'
					subprocess.Popen(shlex.split(script))
				elif args_.host == 'LOCAL':
					script = 'sh ' + paths[3][i] + '00_Clust_' + args_.project_name + '_' + args_.transformation.lower() + '_' + str(args_.k[i]) + '.sh >> ' + paths[1][i] + args_.project_name + '_MICA_Clust.out \n'
					subprocess.Popen(shlex.split(script))
			""""""
	elif args_.mode == 'Reclust':
		paths = setup_directory(args_.outdir, args_.project_name, args_.transformation, args_.k)
		retransform(args_, paths)
		reggplot(args_, paths)
		reclust(args_, paths)
		for i in range(len(paths[3])):
			if args_.host == 'LSF':
				script = 'bsub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MICA_Reclust -q ' + args_.queue + ' -R \"rusage[mem=2000]\" -oo ' + paths[1][i] + args_.project_name + '_MICA_Reclust.out -eo ' + paths[1][i] + args_.project_name + '_MICA_Reclust.err sh ' + paths[3][i] + '00_Reclust_' + args_.project_name + '_' + args_.transformation.lower() + '_' + str(args_.k[i]) + '.sh \n'
				subprocess.Popen(shlex.split(script))
			elif args_.host == 'LOCAL':
				script = 'sh ' + paths[3][i] + '00_Reclust_' + args_.project_name + '_' + args_.transformation.lower() + '_' + str(args_.k[i]) + '.sh >> ' + paths[1][i] + args_.project_name + '_MICA_Reclust.out \n'
				subprocess.Popen(shlex.split(script))
	else:
		print('[EROR] --> [MICA] Unsupported command.')
		exit()

if __name__ == '__main__':
	run(sys.argv)
	print('[INFO] --> [MICA] Finished.')