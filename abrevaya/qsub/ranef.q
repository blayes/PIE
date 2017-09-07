#!/bin/bash
#$ -N abe_ranef
#$ -l mf=16G        
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/pie/abrevaya/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load gurobi

module load matlab/R2015a

matlab -nojvm -nodisplay -singleCompThread -r "calc_wasp_ranef($SGE_TASK_ID, 20, 6, [15 16 17 19 20 23], '/nfsscratch/Users/ssrivastva/pie/abrevaya/')"






