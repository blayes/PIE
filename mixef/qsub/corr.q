#!/bin/bash
#$ -N mix_wasp_ranef
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

module load gurobi/6.5.1

module load matlab/R2015b

matlab -nojvm -nodisplay -singleCompThread -r "calc_wasp_corr($SGE_TASK_ID, 20, 6, [5 6 7 9 10 13], '/Shared/ssrivastva/pie/mixef/')"







