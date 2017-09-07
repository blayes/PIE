#!/bin/bash
#$ -N reg_sim
#$ -l mf=256G
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/wasp/reg/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load R/3.2.1

R CMD BATCH --no-save --no-restore simulate_data.R sim_wasp.rout
