#!/bin/bash
#$ -N mix_wasp_fixef
#$ -l mf=16G        
#$ -l h_rt=320:00:00
#$ -l s_rt=320:00:00
#$ -wd /Users/ssrivastva/pie/mixef/code/
#$ -m a
#$ -M sanvesh-srivastava@uiowa.edu
#$ -t 1-10
#$ -V
#$ -e /Users/ssrivastva/err/
#$ -o /Users/ssrivastva/out/

module load gurobi

module load matlab/R2015a

matlab -nojvm -nodisplay -singleCompThread -r "calc_wasp_fixef($SGE_TASK_ID, 20, 4, 1:4, '/nfsscratch/Users/ssrivastva/pie/mixef/')"





