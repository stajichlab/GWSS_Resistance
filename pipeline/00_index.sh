#!/usr/bin/bash
#SBATCH -p intel,batch --mem 32gb -N 1 -n 8 --out logs/STARindex.%a.log
#SBATCH -D /rhome/cassande/shared/projects/SharpShooter/GWSS_ragtag_genome/RNASeq
#SBATCH -J gwss_index
#SBATCH --mail-type=END # notifications for job done & fail
#SBATCH --mail-user=cassande@ucr.edu # send-to address
#SBATCH --array 1

module load samtools/1.11
module load bwa/0.7.17

if [ -f config.txt ]; then
	source config.txt
fi


if [[ ! -f $REFGENOME.fai || $$REFGENOME -nt $REFGENOME.fai ]]; then
	samtools faidx $REFGENOME
fi
if [[ ! -f $REFGENOME.bwt || $REFGENOME -nt $REFGENOME.bwt ]]; then
	bwa index $REFGENOME
fi


DICT=$(basename $REFGENOME .fasta)".dict"

if [[ ! -f $DICT || $REFGENOME -nt $DICT ]]; then
	rm -f $DICT
	samtools dict $REFGENOME > $DICT
	ln -s $DICT $REFGENOME.dict 
fi

popd
