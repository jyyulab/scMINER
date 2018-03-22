import os, sys

scMINER_PATH = os.environ['scMINER_PATH']
PYTHON_PATH = os.environ['PYTHON_PATH']
sys.path.insert(0, scMINER_PATH)
sys.path.insert(0, PYTHON_PATH)

import utils

if sys.argv[2] == 'tab':
	utils.read_file(sys.argv[1], '\t', sys.argv[3], sys.argv[5], sys.argv[6], sys.argv[4])
elif sys.argv[2] == 'comma':
	utils.read_file(sys.argv[1], ',', sys.argv[3], sys.argv[5], sys.argv[6], sys.argv[4])
elif sys.argv[2] == 'semicolon':
	utils.read_file(sys.argv[1], ';', sys.argv[3], sys.argv[5], sys.argv[6], sys.argv[4])
elif sys.argv[2] == 'pipe':
	utils.read_file(sys.argv[1], '|', sys.argv[3], sys.argv[5], sys.argv[6], sys.argv[4])
utils.patch_file(sys.argv[5] + sys.argv[6] + '.h5.tmp', sys.argv[5], sys.argv[6])
os.remove(sys.argv[5] + sys.argv[6] + '.h5.tmp')