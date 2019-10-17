---
layout: default
title: Installation
nav_order: 2

---
# Installation

## Dependencies
The entry point of scMINER is a R package, which includes some of the major functions such as QC, ... etc. It has been
throughly tested in R 3.5.1. So it is highly recommended to have **R version 3.5.1 or higher** in your environment. 
scMINER also includes two independent python packages [MICA](https://github.com/jyyulab/MICA) and 
[SJARACNe](https://github.com/jyyulab/SJARACNe), for which **python 3.6.1 or higher** is required.

## Installing scMINER

### Install [scMINER](https://github.com/jyyulab/scMINER) R package from github
```R
#install dev_tool first install.packages(devtools)
devtools::install_github("jyyulab/scMINER") 
#or
devtools::install_local('scMINER_0.0.1.tar.gz')
```

### Install [MICA](https://github.com/jyyulab/MICA)

* Installing the official package from PyPi:
```bash
pip install MICA
```

* Or you can install from source: 
```bash
git clone https://github.com/jyyulab/MICA
cd MICA
python setup.py install
```

### Install [SJARACNe](https://github.com/jyyulab/SJARACNe)

* Installing the official package from PyPi:
```bash
pip install SJARACNe
```

* Or you can install from source: 
```bash
git clone https://github.com/jyyulab/SJARACNe
cd SJARACNe
python setup.py build
python setup.py install
```
