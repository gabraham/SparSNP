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

License
-------

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 3 of the License, or
(at your option) any later version.

Copyright (C) 2011-2012 Gad Abraham and National ICT Australia (NICTA).
All rights reserved.

Requirements
------------

   For the post-analysis scripts:  R packages ggplot2 >=0.9.3, scales, grid, abind, ROCR

   A 64-bit operating system is recommended; we have tested SparSNP on 64-bit
   OSX and Linux.


Quick Start
-----------

To get the latest version:

   ```
   git clone git://github.com/gabraham/SparSNP
   ```

To install:

   ```
   cd SparSNP
   make
   ```

Run (assuming a PLINK BED/BIM/FAM dataset named MYDATA, i.e. MYDATA.bim)

   ```
   export PATH=<PATH_TO_SPARSNP>:$PATH
   crossval.sh MYDATA sqrhinge 2>&1 | tee log
   eval.R
   ```


Documentation: see the document workflow.pdf




