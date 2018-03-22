import os, sys

scMINER_PATH = os.environ['scMINER_PATH']
PYTHON_PATH = os.environ['PYTHON_PATH']
sys.path.insert(0, scMINER_PATH)
sys.path.insert(0, PYTHON_PATH)

import utils

trans = sys.argv[1]
if trans == 'MDS':
	utils.mds(int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5], sys.argv[6])
elif trans == 'LPL':
	utils.lpl(int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5], sys.argv[6])
elif trans == 'PCA':
	utils.pca(int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5], sys.argv[6])
elif trans == 'LPCA':
	utils.lpl(int(sys.argv[2]), sys.argv[3], int(sys.argv[4]), sys.argv[5], sys.argv[6])
	utils.pca(int(sys.argv[2]) + 1, sys.argv[3], int(sys.argv[4]), sys.argv[5], sys.argv[6])
