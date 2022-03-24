#!/usr/bin/bash
#SBATCH --nodes 1 --ntasks 1 --mem 40gb --out structure.%a.log
#SBATCH -p intel,batch

module unload miniconda2
module load faststructure

#for different clustering K (trying 1-10)
for kval in {1..30}; do 
  mkdir k_${kval}
  structure.py -K $kval --input=plink --output=k_${kval}/plink
done 

export MPLBACKEND="agg"

#generating plots based on fs.indiv = just lists our individuals (e.g. A.1, A.2, etc)
for kval in {1..30}; do
  distruct.py -K $kval --input=k_${kval}/plink --output=k_${kval}_distruct.svg --popfile=fs.indiv.txt
done

#best k
chooseK.py --input=k_*/plink

