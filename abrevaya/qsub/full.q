#!/bin/bash
#$ -N abe_full
#$ -l mf=64G        
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/pie/abrevaya/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R/3.2.1

R CMD BATCH --no-save --no-restore "--args 2 $SGE_TASK_ID" submit.R full/rep_full_$SGE_TASK_ID.rout

