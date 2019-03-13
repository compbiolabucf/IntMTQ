# IntMTQ
We provide the source code of paper *Platform-integrated mRNA Isoform Quantification* for bioinformatics submission. 

## Goal
 - IntMTQ is an integrative method combining isoform expressions from NanoString/Exon-array platforms to provide better quantification of RNA-Seq based transcript abundances.
 - The code takes aligned Bam files as input and preprocesses RNA-Seq paired-end reads. Then, it performs integrated penalized model and generates isoform expressions.

## Dependencies
The code is written in python 3.6, the following environment is suggested:

 - Python 3.6
 - CVXPY 1.0.6
 - Numpy
 - Pandas
 - Itertools

To install [cvxpy] with conda, run the following command.
```sh
$ conda config --add channels oxfordcontrol
$ conda install -c cvxgrp cvxpy
```

## Prepare RNA-Seq data
[TopHat] is applied to do RNA-Seq short read alignment. With the aligned Bam files and hg19 annotation, we can run the below command to generate q matrix (read counts table).
```sh
$ python read_counts.py
```

## Run IntMTQ with NanoString platform integrated
```sh
$ python IntMTQ.py
```
`IntMTQ.py` will load RNA-Seq data (obtained from RNA-Seq data preparation), NanoString data (**NanoStringData.xlsx** and **GeneName_NanoString.txt**). An excel file **E1_expression.xlsx** will be generated as output which contains the isoform expressions estimated by IntMTQ.       







 [cvxpy]: <https://www.cvxpy.org/index.html>
 [TopHat]: <https://ccb.jhu.edu/software/tophat/index.shtml>
