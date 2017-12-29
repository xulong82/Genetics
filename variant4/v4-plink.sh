#!/bin/bash
#PBS -l nodes=1:ppn=1,walltime=1:00:00

module load plink/1.07

v4="$HOME/Dropbox/GitHub/wgs2/v4/v4.txt"

adsp="$HOME/Dropbox/GitHub/wgs2/v4/adsp"
jaxcs="$HOME/Dropbox/GitHub/wgs2/v4/jaxcs"

# by JAX CS

cd /data/xwang/adsp/plink2

plink --noweb --bfile autosome --extract $v4 --range --recodeA --out ${jaxcs} 

# by ADSP

cd /data/xwang/adsp3/plink

plink --noweb --bfile wgs --extract $v4 --range --recodeA --out ${adsp} 

