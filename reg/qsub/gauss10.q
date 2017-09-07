#!/bin/bash
#$ -N reg_bary10
#$ -l mf=16G
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/wasp/reg/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-120
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R/3.2.1

R CMD BATCH --no-save --no-restore "--args 7 $SGE_TASK_ID" submit.R lap/bary10_$SGE_TASK_ID.rout
