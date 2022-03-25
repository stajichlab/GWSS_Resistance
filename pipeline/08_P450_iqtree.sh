#!/usr/bin/bash
#SBATCH --ntasks 24 -N 1 --mem 48gb --time 12-0:00:00 -p intel --out iqtree.%A.log

CPU=$SLURM_CPUS_ON_NODE
if [ -z $CPU ]; then
	CPU=1
fi

module load iqtree/2.1.3 
module load muscle/3.8.425
module load  trimal/1.4.1


muscle -in p450_all_T1_and_any_T2_expressed.fa -out p450_all_T1_and_any_T2_expressed.aln
trimal -in p450_all_T1_and_any_T2_expressed.aln -out p450_all_T1_and_any_T2_expressed_trim.aln -automated1
iqtree2 -s p450_all_T1_and_any_T2_expressed_trim.aln -bb 1000 -nt AUTO


