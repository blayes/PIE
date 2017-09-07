#!/bin/bash
#$ -N reg_wasp20
#$ -l mf=16G
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/wasp/reg/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-2400
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R/3.2.1

R CMD BATCH --no-save --no-restore "--args 5 $SGE_TASK_ID" submit.R wasp/wasp20_$SGE_TASK_ID.rout
