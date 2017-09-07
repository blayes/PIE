#!/bin/bash
#$ -N abe_part
#$ -l mf=16G        
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/pie/abrevaya/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-200
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R/3.2.1

R CMD BATCH --no-save --no-restore "--args 4 $SGE_TASK_ID" submit.R part/rep_part_$SGE_TASK_ID.rout


