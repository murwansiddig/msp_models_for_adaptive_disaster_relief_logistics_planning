#/bin/bash
 
#PBS -N hurricane
#PBS -l select=1:ncpus=24:mem=1000gb,walltime=168:00:00
#PBS -q bigmem
#PBS -M msiddig@clemson.edu

module add gnu-parallel
module add julia/1.6.1-gcc/8.3.1
module add gurobi

cd $PBS_O_WORKDIR

parallel < commands.txt
