import os, sys

scMINER_PATH = os.environ['scMINER_PATH']
PYTHON_PATH = os.environ['PYTHON_PATH']
sys.path.insert(0, scMINER_PATH)
sys.path.insert(0, PYTHON_PATH)

import utils
import pandas as pd

mat1 = pd.HDFStore(sys.argv[1], 'r')[sys.argv[3]]
mat2 = pd.HDFStore(sys.argv[2], 'r')[sys.argv[4]]
bins = int(sys.argv[5])
m = int(sys.argv[6])
name = sys.argv[7]

utils.calc_mi_mat(mat1, mat2, bins, m, name, sys.argv[8], sys.argv[9] + '_' + name)
os.remove(sys.argv[1])
os.remove(sys.argv[2])