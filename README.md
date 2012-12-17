SparSNP
=======

SparSNP fits lasso-penalized linear models to SNP data. Its main features are:

* it can fit squared hinge loss for classification (case/control) and linear regression (quantitative phenotypes)
* takes PLINK BED/FAM files as input
* the amount of memory is bounded - can work with large datasets using little memory (typically <1GB, more for better performance)
* fits a model over a grid of penalties, and writes the estimated coefficients to disk
* it can also do cross-validation, using the estimated coefficients to predict outputs for other datasets
* efficient - it uses warm-restarts plus an active-set approach, the model fitting part of 3-fold cross-validation for a dataset of 2000
   samples by 300,000 SNP dataset takes ~5min, and about 25min for ~6800 samples / ~516,000 SNPs

Contact
-------

Gad Abraham, gad.abraham@unimelb.edu.au

Citation
--------

G. Abraham, A. Kowalczyk, J. Zobel, and M. Inouye, ``SparSNP: Fast and
memory-efficient analysis of all SNPs for phenotype prediction'', BMC
Bioinformatics,
2012, 13:88, [doi:10.1186/1471-2105-13-88](http://www.biomedcentral.com/1471-2105/13/88/)

Copyright (C), Gad Abraham and National ICT Australia (2011-2012), All Rights Reserved. 

Quick setup
-----------

To get the latest version:

    git clone git://github.com/gabraham/SparSNP

To install:

    cd SparSNP
    make

Requirements for the post-analysis scripts:
   R packages ggplot2, scales, grid, abind, ROCR

Documentation: see the document workflow.pdf




