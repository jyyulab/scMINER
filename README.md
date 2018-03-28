# scMINER
## single-cell Mutual Information Network Engineered Ranger

scMINER is a toolbox for single-cell analysis based on mutual information. This package is a combination of several individual tools including but not limited to MIE and MICA. Other tools including network construction, master regultors identification and activity based analysis may be added. MIE is the baseline of the package and the other tools require the output of MIE tool. MICA in this package has been improved and the main advantage of that is parallelism using St. Jude research cluster and Phoenix application.

## Download

<code>git clone https://github.com/jyyulab/scMINER.git</code>

## Starting Point

To start using the package on St. Jude research cluster (by default), follow the below steps:

### Phoenix Setup

In St. Jude's Research Cluster, Phoenix is accessible via:

<code>setcbenv prod</code></br>
<code>set official</code></br>
<code>cbload phoenix</code>

And the following two environment variables for scMINER:

### Setup scMINER

Set the two environment variables as follows:

<code>export PYTHON_PATH=[path_to_target_python3]</code>
</br>
<code>export scMINER_PATH=[path_to_root_directory_of_scMINER]</code>

### Setup R

To generate R plot load R module as follows:

<code>module load R/3.4.0</code>

To use the package in a regular machine, just follow Setup scMINER section and make sure to have R and requirement Python packages installed and when running any command set the --host option to LOCAL.

## Run MIE

To run MIE, simply follow the following example:

<code>python3 scMINER.py MIE Pipeline [project_name] [input_tab_separated_file] [output_dir] [output_file_name]</code>

To avoid memory issue in case that the input file is large, make sure to use bsub on research cluster.

## Run MICA

To run MICA in clustering mode, simply follow the following example:

<code>python3 scMINER.py MICA Clust [project_name] [path_to: project_name.whole.h5] [path to: project_name_mi.h5] [output_dir] [project_name] --k 2 3 4 ... </code>
