
################################################################################
# SparSNP --- Sparse SNP analysis
#
# Contact: Gad Abraham, gabraham@csse.unimelb.edu.au
#
# Copyright (c) 2011, NICTA
#
################################################################################

(1) Introduction
================

SparSNP is a set of tools for fitting large scale lasso models to SNP data. It
can be used to build classifiers for case/control data.

(2) Limitations
===============

Theoretically, SparSNP can handle up to 2^31 - 1 SNPs and samples,
since it uses signed ints internally for array indicies, but has never
been tested with that. We have tested sample sizes up to about 12,000,
and SNPs up to about 900,000 without any issue.

(3) Installation
================
Operating System:
   * Unix-like environment with Bash, tested on 64-bit Ubuntu 10.04 and OSX 10.6.8

Requirements:
   * bash, unix tools (wc)

Requirements for code to analyse results:
   * R, with packages:
      - ggplot2
      - glmnet

(4) Usage
===============

Compile the code or use the prebuilt binaries.

Assuming you have a PLINK dataset in binary format (BED/BIM/FAM files) called
discovery, i.e., you have
   discovery.bed
   discovery.bim
   discovery.fam

IMPORTANT: the code assumes that the binary data are coded using the
default PLINK encoding (which is in effect major allele dosage). The
definition of which allele is major and which is minor may change
between different datasets.  Therefore, we recommend using PLINK to
either merge them into one
(http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#merge), so
as to make the encoding consistent, and then to split them, or to
change the reference allele in one to match the other
(http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#refallele).

ALSO IMPORTANT: SparSNP currently does NOT handle missing phenotypes
values such as -9. These should be excluded with PLINK prior to
running SparSNP. Missing genotypes are OK, they will be uniformly
randomly imputed (but be aware that missingness that is too high or
not random will impact the results and may produce spurious
associations).


(4.1) Discovery dataset
-----------------------

Cross-validation: use the helper script cv3.sh to run 3-fold cross-validation, 

   ./cv3.sh discovery sqrhinge

where "sqrhinge" is for classification

or for linear regression
   ./cv3.sh discovery linear

Each cross-validation replication will be stored in the directory
"discovery", in it's own crossvalXX directory. By default, the script
will skip existing directories to prevent over-writing old results.

Results: in each crossvalXX directory, there will be files:
   folds.ind, a binary file representing the split into crossval folds
   folds.txt, a text file representing which fold each sample belongs to
   scale.bin.XX, a binary file of means and std devs for each SNP in each fold
	 XX
   lambda1path.csv.XX, a text file of the grid of lambda1 penalties in each
	 fold XX
   beta.csv.YY.XX, a text file of estimated coefficients in fold XX for the
	 YYth penalty
   beta.csv.YY.XX.pred, a text file of predicted outputs (classes/responses)
	 in fold XX for the YYth penalty for the *test data*
   y.XX: the *test* output/classes/responses for the XXth fold

(4.2) Validation dataset
------------------------

   ./validation.sh validation

Will produce a directory called "validation" with the predictions of
the discovery models on the validation dataset.


(5) Analysis after Model Fitting
=================================

(5.1) Discovery dataset
-----------------------

To produce plots of AUC for classification and R^2 for regression in
the discovery data (cross-validated), use

   ./eval.R

which will produce a plot "discovery_AUC.pdf" in the directory "discovery"

Optionally, prevalence (0 to 1) can be supplied to plot the explained
phenotypic variance for case/control data

   ./eval.R prev=0.01

which will also produce a plot "discovery_VarExp.pdf".

(5.2) Validation dataset
-----------------------

To plot AUC/R^2 and/or phenotypic variance for the validation dataset, use

   ./eval.R mode=validation

and

   ./eval.R mode=validation prev=0.01




