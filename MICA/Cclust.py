import os, sys

scMINER_PATH = os.environ['scMINER_PATH']
PYTHON_PATH = os.environ['PYTHON_PATH']
sys.path.insert(0, scMINER_PATH)
sys.path.insert(0, PYTHON_PATH)

import utils

utils.cc(sys.argv[3], int(sys.argv[4]), sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8], sys.argv[9], int(sys.argv[10]), sys.argv[11])
#utils.plot(int(sys.argv[1]), sys.argv[2], sys.argv[8], sys.argv[9])