#!/bin/bash
#PBS -l nodes=1:ppn=15,walltime=6:00:00

module load R
# module load R/3.2.3

# w/o parameter
  Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/wgs2/Manu/R/gwas_qtlrel_vs.R

# w/ parameter
# Rscript --no-restore --quiet /home/xwang/Dropbox/GitHub/wgs2/Manu/R/myQtlRel.R ${chr}

