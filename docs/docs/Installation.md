---
layout: default
title: Installation
nav_order: 2

---
#Installation

### Dependencies
scMINER is depend on [python3](https://www.python.org/downloads/) and [R](https://www.r-project.org/), including three individual modules: scMINER R package, MICA and SJARACNe python package. 


### Local installation: 
1. Install [scMINER](https://github.com/jyyulab/scMINER) from github:

```R
#install dev_tool first install.packages(devtools)
devtools::install_github("jyyulab/scMINER") 

#or
devtools::install_local('scMINER_0.0.1.tar.gz')
```

2. Install [MICA](https://github.com/jyyulab/MICA) from source:
```
$ git clone https://github.com/jyyulab/MICA
$ cd MICA
$ python setup.py install
```


3. Install [SJARACNe](https://github.com/jyyulab/SJARACNe) from source:
```
$ git clone https://github.com/jyyulab/SJARACNe
$ cd SJARACNe
$ python setup.py install
```