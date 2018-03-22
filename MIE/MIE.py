import os, sys

scMINER_PATH = os.environ['scMINER_PATH']
PYTHON_PATH = os.environ['PYTHON_PATH']
sys.path.insert(0, scMINER_PATH)
sys.path.insert(0, PYTHON_PATH)

import utils
import argparse
import numpy as np
import pandas as pd
import shlex, subprocess


def setup(args):
	parser = argparse.ArgumentParser(description='MIE package.')
	parser.add_argument('mode', help='Procedure to select and execute from the package [Pipeline | Read | Slice | Calc | Merge | Norm]')
	parser.add_argument('project_name', help='Project name')
	args_ = parser.parse_args(args[1:3])
	if args_.mode == 'Read':
		parser.add_argument('input_file', help='Input text-based file name')
		parser.add_argument('--delimiter', default='tab', help='Delimiter used in the input file [tab | comma | semicolon | pipe] (default: tab)')
		parser.add_argument('--header', default='True', help='Whether the input contains header line or not [True | False]')
		parser.add_argument('--index_col', type=int, default=0, help='Column index to be used as row index and the starting poing of actual data [ >= 0 ] (default: 0)')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, default=2000, help='Memory allocation to read step (default: 2GB)')
		parser.add_argument('--queue', default='compbio', help='Queue name for job allocation')
	elif args_.mode == 'Slice':
		parser.add_argument('input_file', help='Input hdf file name')
		parser.add_argument('--slice_size', type=int, default=1000, help='Size for each slice (default: 1000)')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, default=2000, help='Memory allocation to slice step (default: 2GB)')
		parser.add_argument('--queue', default='compbio', help='Queue name for job allocation')
	elif args_.mode == 'Calc':
		parser.add_argument('input_file', help='Input hdf file name')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, default=2000, help='Memory allocation to calc step (default: 2GB)')
		parser.add_argument('--queue', default='compbio', help='Queue name for job allocation')
	elif args_.mode == 'Merge':
		parser.add_argument('input_dir', help='Input directory containing all partial mutual information hdf files')
		parser.add_argument('common_name', help='Shared name between sliced matrices')
		parser.add_argument('number_slice', type=int, help='Number of slices produced by slice procedure')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, default=2000, help='Memory allocation to merge step (default: 2GB)')
		parser.add_argument('--queue', default='compbio', help='Queue name for job allocation')
	elif args_.mode == 'Norm':
		parser.add_argument('input_file', help='Input complete hdf mutual information file name')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, default=2000, help='Memory allocation to norm step (default: 2GB)')
		parser.add_argument('--queue', default='compbio', help='Queue name for job allocation')
	elif args_.mode == 'Pipeline':
		parser.add_argument('input_file', help='Input text-based file name')
		parser.add_argument('--delimiter', default='tab', help='Delimiter used in the input file [tab | comma | semicolon | pipe] (default: tab)')
		parser.add_argument('--header', default='True', help='Whether the input contains header line or not [True | False]')
		parser.add_argument('--index_col', type=int, default=0, help='Column index to be used as row index and the starting poing of actual data [ >= 0 ] (default: 0)')
		parser.add_argument('--slice_size', type=int, default=1000, help='Size for each slice (default: 1000)')
		parser.add_argument('outdir', help='Output directory')
		parser.add_argument('outfilename', help='Common name used for all outputs')
		parser.add_argument('--host', default='LSF', help='Computation host of the jobs [LOCAL | LSF] (default: LSF)')
		parser.add_argument('--resource', type=int, nargs='+', default=[2000] * 5, help='Memory allocation to each step in the following order: read, slice, calc, merge, and norm (default: 2GB)')
		parser.add_argument('--queue', default='compbio', help='Queue name for job allocation')
	args_ = parser.parse_args(args[1:])
	return args_

def setup_directory(out_dir, project_name):
	scMINER_path = out_dir + 'scMINER_' + project_name + '/'
	log_path = scMINER_path + 'scMINER_MIE_log/'
	out_path = scMINER_path + 'scMINER_MIE_out/'
	script_path = scMINER_path + 'scMINER_MIE_scripts/'
	if not os.path.exists(scMINER_path):
		os.mkdir(scMINER_path)
	if not os.path.exists(log_path):
		os.mkdir(log_path)
	if not os.path.exists(out_path):
		os.mkdir(out_path)
	if not os.path.exists(script_path):
		os.mkdir(script_path)
	return [scMINER_path, log_path, out_path, script_path]

def prep(args, paths):
	if args.delimiter == 'tab':
		utils.read_file(args.input_file, '\t', args.header, paths[2], args.outfilename, str(args.index_col))
	elif args.delimiter == 'comma':
		utils.read_file(sys.argv[1], ',', args.header, paths[2], args.outfilename, str(args.index_col))
	elif args.delimiter == 'semicolon':
		utils.read_file(sys.argv[1], ';', args.header, paths[2], args.outfilename, str(args.index_col))
	elif args.delimiter == 'pipe':
		utils.read_file(sys.argv[1], '|', args.header, paths[2], args.outfilename, str(args.index_col))
	utils.patch_file(paths[2] + args.outfilename + '.h5.tmp', paths[2], args.outfilename)
	os.remove(paths[2] + args.outfilename + '.h5.tmp')

	args.input_file = paths[2] + args.outfilename + '.whole.h5'
	utils.slice_file(args.input_file, paths[2], args.outfilename, str(args.slice_size))

def read(args, paths):
	script = PYTHON_PATH + ' ' + scMINER_PATH + 'MIE/Read.py ' + ' ' + args.input_file + ' ' + args.delimiter + ' ' + args.header + ' ' + str(args.index_col) + ' ' + paths[2] + ' ' + args.outfilename
	out_1 = open(paths[3] + '01_Read_' + args.project_name + '.sh', 'w')
	out_1.write(script + '\n')
	out_1.close()

def slice(args, paths):
	script = PYTHON_PATH + ' ' + scMINER_PATH + 'MIE/Slice.py ' + ' ' + args.input_file + ' ' + paths[2] + ' ' + args.outfilename + ' ' + str(args.slice_size)
	out_2 = open(paths[3] + '02_Slice_' + args.project_name + '.sh', 'w')
	out_2.write(script + '\n')
	out_2.close()

def calc(args, paths):
	path_tmp = paths[2] + '.tmp/'
	if not os.path.exists(path_tmp):
		os.mkdir(path_tmp)
	in_ = pd.HDFStore(args.input_file, 'r')
	bins = int(np.floor(in_['attr'].loc['row'] ** (1/3.0)))
	b = in_['attr'].loc['slice', 0]
	m = in_['attr'].loc['col', 0]
	digit = int(np.floor(np.log10(b)) + 1)
	total = int((b * (b + 1)) / 2)
	digit1 = int(np.floor(np.log10(total)) + 1)

	out_3 = open(paths[3] + '03_Calc_' + args.project_name + '.sh', 'w')
	for i in range(b):
		key1 = 'slice_' + str(i).zfill(digit)
		mat1 = in_[key1]
		for j in range(i, b):
			key2 = 'slice_' + str(j).zfill(digit)
			mat2 = in_[key2]

			mat1_tmp = path_tmp + args.outfilename + '.slice_' + str(i).zfill(digit) + '_' + str(j).zfill(digit) + '_1.h5.tmp'
			mat2_tmp = path_tmp + args.outfilename + '.slice_' + str(i).zfill(digit) + '_' + str(j).zfill(digit) + '_2.h5.tmp'
			mat1.to_hdf(mat1_tmp, key1)
			mat2.to_hdf(mat2_tmp, key2)

			idx = int(i * b + j - (i * (i + 1)) / 2)
			name = 'mi_' + str(idx).zfill(digit1)
			script = PYTHON_PATH + ' ' + scMINER_PATH + 'MIE/Calc.py ' + ' ' + mat1_tmp + ' ' + mat2_tmp + ' ' + key1 + ' ' + key2 + ' ' + str(bins) + ' ' + str(m) + ' ' + name + ' ' + paths[2] + ' ' + args.outfilename + ' '
			out_3.write(script)
			if args.host == 'LOCAL':
				out_3.write('&')
			out_3.write('\n')

	out_3.close()
	args.number_slice = b

def merge(args, paths):
	script = PYTHON_PATH + ' ' + scMINER_PATH + 'MIE/Merge.py ' + ' ' + args.input_dir + ' ' + args.common_name + ' ' + str(args.number_slice) + ' ' + paths[2] + ' ' + args.outfilename
	out_4 = open(paths[3] + '04_Merge_' + args.project_name + '.sh', 'w')
	out_4.write(script + '\n')
	out_4.close()

def norm(args, paths):
	script = PYTHON_PATH + ' ' + scMINER_PATH + 'MIE/Norm.py ' + ' ' + args.input_file + ' ' + paths[2] + ' ' + args.outfilename
	out_5 = open(paths[3] + '05_Norm_' + args.project_name + '.sh', 'w')
	out_5.write(script + '\n')
	out_5.close()

def pipeline(args, paths):
	prep(args, paths)
	args.input_file = paths[2] + args.outfilename + '.sliced.h5'
	calc(args, paths)
	args.input_dir = paths[2]
	args.common_name = args.outfilename
	merge(args, paths)
	args.input_file = paths[2] + args.outfilename + '_mi.h5'
	norm(args, paths)
	out_0 = open(paths[3] + '00_Pipeline_' + args.project_name + '.sh', 'w')
	if args.host == 'LSF':
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MIE_Calc -q ' + args.queue + ' -M ' + str(args.resource[2]) + ' -i ' + paths[3] + '03_Calc_' + args.project_name + '.sh -oo ' + paths[1] + args.project_name + '_MIE_Calc.%J.%I.out -eo ' + paths[1] + args.project_name + '_MIE_Calc.%J.%I.err \n'
		out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MIE_Merge -q ' + args.queue + ' -M ' + str(args.resource[3]) + ' -i ' + paths[3] + '04_Merge_' + args.project_name + '.sh -oo ' + paths[1] + args.project_name + '_MIE_Merge.%J.%I.out -eo ' + paths[1] + args.project_name + '_MIE_Merge.%J.%I.err \n'
		out_0.write(script)
		script = 'psub -K -P ' + args.project_name + ' -J ' + args.project_name + '_MIE_Norm -q ' + args.queue + ' -M ' + str(args.resource[4]) + ' -i ' + paths[3] + '05_Norm_' + args.project_name + '.sh -oo ' + paths[1] + args.project_name + '_MIE_Norm.%J.%I.out -eo ' + paths[1] + args.project_name + '_MIE_Norm.%J.%I.err \n'
		out_0.write(script)
	elif args.host == 'LOCAL':
		script = 'sh ' + paths[3] + '03_Calc_' + args.project_name + '.sh >> ' + paths[1] + args.project_name + '_MIE_Calc.out \n'
		out_0.write(script)
		out_0.write('jobs=$(ps -ef | grep ' + args.project_name + ' | grep Calc.py -c)\n')
		out_0.write('while [ $jobs -gt 0 ]\ndo\n\tsleep 30\n\tjobs=$(ps -ef | grep ' + args.project_name + ' | grep Calc.py -c)\ndone\n')
		script = 'sh ' + paths[3] + '04_Merge_' + args.project_name + '.sh >> ' + paths[1] + args.project_name + '_MIE_Merge.out \n'
		out_0.write(script)
		script = 'sh ' + paths[3] + '05_Norm_' + args.project_name + '.sh >> ' + paths[1] + args.project_name + '_MIE_Norm.out \n'
		out_0.write(script)
	out_0.close()

def run(args):
	args_ = setup(args)
	paths = setup_directory(args_.outdir, args_.project_name)
	if args_.mode == 'Read':
		read(args_, paths)
		if args_.host == 'LSF':
			script = 'psub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MIE_Read -q ' + args_.queue + ' -M ' + str(args_.resource) + ' -i ' + paths[3] + '01_Read_' + args_.project_name + '.sh -oo ' + paths[1] + args_.project_name + '_MIE_Read.%J.%I.out -eo ' + paths[1] + args_.project_name + '_MIE_Read.%J.%I.err \n'
			subprocess.Popen(shlex.split(script))
		elif args_.host == 'LOCAL':
			script = 'sh ' + paths[3] + '01_Read_' + args_.project_name + '.sh >> ' + paths[1] + args_.project_name + '_MIE_Read.out \n'
			subprocess.Popen(shlex.split(script))
	elif args_.mode == 'Slice':
		slice(args_, paths)
		if args_.host == 'LSF':
			script = 'psub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MIE_Slice -q ' + args_.queue + ' -M ' + str(args_.resource) + ' -i ' + paths[3] + '02_Slice_' + args_.project_name + '.sh -oo ' + paths[1] + args_.project_name + '_MIE_Slice.%J.%I.out -eo ' + paths[1] + args_.project_name + '_MIE_Slice.%J.%I.err \n'
			subprocess.Popen(shlex.split(script))
		elif args_.host == 'LOCAL':
			script = 'sh ' + paths[3] + '02_Slice_' + args_.project_name + '.sh >> ' + paths[1] + args_.project_name + '_MIE_Slice.out \n'
			subprocess.Popen(shlex.split(script))
	elif args_.mode == 'Calc':
		calc(args_, paths)
		if args_.host == 'LSF':
			script = 'psub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MIE_Calc -q ' + args_.queue + ' -M ' + str(args_.resource) + ' -i ' + paths[3] + '03_Calc_' + args_.project_name + '.sh -oo ' + paths[1] + args_.project_name + '_MIE_Calc.%J.%I.out -eo ' + paths[1] + args_.project_name + '_MIE_Calc.%J.%I.err \n'
			subprocess.Popen(shlex.split(script))
		elif args_.host == 'LOCAL':
			script = 'sh ' + paths[3] + '03_Calc_' + args_.project_name + '.sh >> ' + paths[1] + args_.project_name + '_MIE_Calc.out \n'
			subprocess.Popen(shlex.split(script))
	elif args_.mode == 'Merge':
		merge(args_, paths)
		if args_.host == 'LSF':
			script = 'psub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MIE_Merge -q ' + args_.queue + ' -M ' + str(args_.resource) + ' -i ' + paths[3] + '04_Merge_' + args_.project_name + '.sh -oo ' + paths[1] + args_.project_name + '_MIE_Merge.%J.%I.out -eo ' + paths[1] + args_.project_name + '_MIE_Merge.%J.%I.err \n'
			subprocess.Popen(shlex.split(script))
		elif args_.host == 'LOCAL':
			script = 'sh ' + paths[3] + '04_Merge_' + args_.project_name + '.sh >> ' + paths[1] + args_.project_name + '_MIE_Merge.out \n'
			subprocess.Popen(shlex.split(script))
	elif args_.mode == 'Norm':
		norm(args_, paths)
		if args_.host == 'LSF':
			script = 'psub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MIE_Norm -q ' + args_.queue + ' -M ' + str(args_.resource) + ' -i ' + paths[3] + '05_Norm_' + args_.project_name + '.sh -oo ' + paths[1] + args_.project_name + '_MIE_Norm.%J.%I.out -eo ' + paths[1] + args_.project_name + '_MIE_Norm.%J.%I.err \n'
			subprocess.Popen(shlex.split(script))
		elif args_.host == 'LOCAL':
			script = 'sh ' + paths[3] + '05_Norm_' + args_.project_name + '.sh >> ' + paths[1] + args_.project_name + '_MIE_Norm.out \n'
			subprocess.Popen(shlex.split(script))
	elif args_.mode == 'Pipeline':
		pipeline(args_, paths)
		if args_.host == 'LSF':
			script = 'bsub -P ' + args_.project_name + ' -J ' + args_.project_name + '_MIE_Pipeline -q ' + args_.queue + ' -R \"rusage[mem=2000]\" -oo ' + paths[1] + args_.project_name + '_MIE_Pipeline.out -eo ' + paths[1] + args_.project_name + '_MIE_Pipeline.err sh ' +  paths[3] + '00_Pipeline_' + args_.project_name + '.sh \n'
			subprocess.Popen(shlex.split(script))
		elif args_.host == 'LOCAL':
			script = 'sh ' + paths[3] + '00_Pipeline_' + args_.project_name + '.sh >> ' + paths[1] + args_.project_name + '_MIE_Pipeline.out \n'
			subprocess.Popen(shlex.split(script))
	else:
		print('[EROR] --> [MIE] Unsupported command.')
		exit()

if __name__ == '__main__':
	run(sys.argv)
	print('[INFO] --> [MIE] Finished.')
