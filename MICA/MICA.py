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
		parser.add_argument('--retransformation', default='True', help='Whether apply re-transformation for final visualization or not (default: True)')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, nargs='+', default=[2000]*5, help='Memory allocation for each individual step in clustering pipeline (default: 2GB)')
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
		script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Transform.py ' + ' ' + args.transformation + ' ' + str(fig_num) + ' ' + path_tmp + args.dist_file.split('/')[-1] + ' ' + str(args.max_dim) + ' ' + paths[2][i] + ' ' + args.outfilename + ' '
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
			script = 'cp ' + paths[2][i] + args.project_name + '_clust.h5 ' + path_tmp + args.project_name + '_clust.h5.tmp.' + str(j) + ' '
			out_3.write(script + '\n')
			if args.transformation == 'LPCA':
				script = 'cp ' + paths[2][i] + args.project_name + '_clust.h5 ' + path_tmp + args.project_name + '_clust.h5.tmp.' + str(iterations + j) + ' '
				out_3.write(script + '\n')
		for j in range(args.bootstrap):
			dim = args.max_dim - args.min_dim + 1
			for k in range(dim):
				if args.transformation == 'MDS':
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.project_name + '_clust.h5.tmp.' + str(j * dim + k) + ' mds ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.project_name + '_kmeans.h5.tmp.' + str(j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
				elif args.transformation == 'PCA':
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.project_name + '_clust.h5.tmp.' + str(j * dim + k) + ' pca ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.project_name + '_kmeans.h5.tmp.' + str(j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
				elif args.transformation == 'LPL':
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.project_name + '_clust.h5.tmp.' + str(j * dim + k) + ' lpl ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.project_name + '_kmeans.h5.tmp.' + str(j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
				elif args.transformation == 'LPCA':
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.project_name + '_clust.h5.tmp.' + str(j * dim + k) + ' pca ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.project_name + '_kmeans.h5.tmp.' + str(j * dim + k) + ' '
					if args.host == 'LOCAL':
						script += '&'
					out_4.write(script + '\n')
					script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Kmeans.py ' + ' ' + path_tmp + args.project_name + '_clust.h5.tmp.' + str(iterations + j * dim + k) + ' lpl ' + str(args.k[i]) + ' ' + str(k + args.min_dim) + ' ' + path_tmp + ' ' + args.project_name + '_kmeans.h5.tmp.' + str(iterations + j * dim + k) + ' '
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
		script = PYTHON_PATH + ' ' + scMINER_PATH + 'MICA/Cclust.py ' + ' ' + str(fig_num) + ' ' + paths[2][i] + args.project_name + '_clust.h5 ' + path_tmp + ' ' + str(args.k[i]) + ' ' + args.project_name + ' kmeans.h5.tmp ' + args.transformation.lower() + ' ' + paths[2][i] + ' ' + args.outfilename + ' ' + str(args.max_dim) + ' ' + args.retransformation + ' ' 
		out_5 = open(paths[3][i] + '05_CClust_' + args.project_name + '_' + args.transformation.lower() + '_' + str(args.k[i]) + '.sh', 'w')
		out_5.write(script + '\n')
		out_5.close()
		fig_num += 3

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
			out_0.write('jobs=$(ps -ef | grep \"' + args.project_name + '_' + args.transformation + '_' + str(args.k[i]) + '\" | grep Cclust.py -c)\n')
			out_0.write('while [ $jobs -gt 0 ]\ndo\n\tsleep 30\n\tjobs=$(ps -ef | grep \"' + args.project_name + '_' + args.transformation + '_' + str(args.k[i]) + '\" | grep Cclust.py -c)\ndone\n')
		#out_0.write('rm -rf ' + paths[2][i] + '.tmp/ \n')
		out_0.close()

def run(args):
	args_ = setup(args)
	if args_.mode == 'Clust':
		paths = setup_directory(args_.outdir, args_.project_name, args_.transformation, args_.k)
		prep(args_, paths)		
		hclust(args_, paths)
		transform(args_, paths)
		kmeans(args_, paths)
		cclust(args_, paths)
		clust(args_, paths)
		for i in range(len(paths[3])):
			if args_.host == 'LSF':
				script = 'bsub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MICA_Clust -q ' + args_.queue + ' -R \"rusage[mem=2000]\" -oo ' + paths[1][i] + args_.project_name + '_MICA_Clust.out -eo ' + paths[1][i] + args_.project_name + '_MICA_Clust.err sh ' + paths[3][i] + '00_Clust_' + args_.project_name + '_' + args_.transformation.lower() + '_' + str(args_.k[i]) + '.sh \n'
				subprocess.Popen(shlex.split(script))
			elif args_.host == 'LOCAL':
				script = 'sh ' + paths[3][i] + '00_Clust_' + args_.project_name + '_' + args_.transformation.lower() + '_' + str(args_.k[i]) + '.sh >> ' + paths[1][i] + args_.project_name + '_MICA_Clust.out \n'
				subprocess.Popen(shlex.split(script))
	else:
		print('[EROR] --> [MICA] Unsupported command.')
		exit()

if __name__ == '__main__':
	run(sys.argv)
	print('[INFO] --> [MICA] Finished.')