import os, sys

scMINER_PATH = os.environ['scMINER_PATH']
PYTHON_PATH = os.environ['PYTHON_PATH']
sys.path.insert(0, scMINER_PATH)
sys.path.insert(0, PYTHON_PATH)

import utils

utils.merge_mi_mats(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5])