#!/usr/bin/bash
#SBATCH -p intel -N 1 -n 8 --mem 32gb --out logs/realign.%a.log
#SBATCH --array 1-213

module load picard
module load samtools
module load java


MEM=32g
TEMP=temp
ALN=results/STAR_int
GVCFDIR=gvcf_int
mkdir -p $GVCFDIR
SAMPLEFILE=samples.tsv
GENOME=$(realpath data/Homalodisca_vitripennis_A6A7A9_masurca_v1_ragtag_v1.scaffolds.fa)
GENOMEDICT=$(echo $GENOME | perl -p -e 's/\.fasta/.dict/')
if [[ ! -f $GENOME.fai || $GENOME -nt $GENOME.fai ]]; then
  samtools faidx $GENOME
fi

if [[ ! -f $GENOMEDICT || $GENOME -nt $GENOMEDICT ]]; then
  picard CreateSequenceDictionary R=$GENOME O=$GENOMEDICT
fi

GFF=$(realpath data/Homalodisca_vitripennis_A6A7A9_masurca_v1_ragtag_v1.gff3)
GTF=Homalodisca_vitripennis_A6A7A9_masurca_v1_ragtag_v1.gtf

CPU=1
if [ $SLURM_CPUS_ON_NODE ]; then
    CPU=$SLURM_CPUS_ON_NODE
fi

N=${SLURM_ARRAY_TASK_ID}

if [ -z $N ]; then
    N=$1
    if [ -z $N ]; then
        echo "Need an array id or cmdline val for the job"
        exit
    fi
fi

hostname
date
tail -n +2 $SAMPLEFILE |  sed -n ${N}p | while read SAMPLE NAME REP LOCATION STATUS
do
 INNAME=$NAME.${REP}
 INBAM=$ALN/${INNAME}.Aligned.sortedByCoord.out.bam
 RGBAM=$TEMP/${INNAME}.RG.bam
 DDBAM=$TEMP/${INNAME}.DD.bam
 SPLITBAM=$ALN/${INNAME}.split.bam
 GVCF=$GVCFDIR/${INNAME}.g.vcf
 if [[ ! -f $SPLITBAM || $INBAM -nt $SPLITBAM ]]; then
   if [[ ! -f $DDBAM || $INBAM -nt $DDBAM ]]; then
     if [[ ! -f $RGBAM || $INBAM -nt $RGBAM ]]; then
       picard AddOrReplaceReadGroups I=$INBAM O=$RGBAM SO=coordinate RGID=$INNAME RGLB=$SAMPLE RGPL=illumina RGPU=$SAMPLE.$REP RGSM=$INNAME
     fi
     picard MarkDuplicates I=$RGBAM O=$DDBAM CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=logs/${INNAME}.output.metrics
   fi
   module load gatk/4.1.8.1
   time gatk SplitNCigarReads --reference $GENOME --input $DDBAM --output $SPLITBAM
   #-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
   module unload gatk
   if [ -f $SPLITBAM ]; then
     samtools index $SPLITBAM
     #rm -f $RGBAM $DDBAM $DDBAM.bai
   fi
   if [ ! -f $SPLITBAM.bai ]; then
     echo "Error or failure in creating Split BAM $SPLITBAM $SPLITBAM.bai don't exist"
     exit
   fi
   # make gvcf
 fi
 echo "$SPLITBAM $GENOME $GVCF"
 if [[ ! -f $GVCF || $SPLITBAM -nt $GVCF ]]; then
   module load gatk/4
#   time java -Xmx${MEM} -jar $GATK  -T HaplotypeCaller \
#   	    -ERC GVCF -ploidy 2 \
#   	    -I $SPLITBAM -R $GENOME \
#   	    -o $GVCF -nct $CPU
   time gatk --java-options "-Xmx${MEM}" HaplotypeCaller -ERC GVCF -I  $SPLITBAM -R $GENOME \
	   -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation \
	   -O $GVCF --native-pair-hmm-threads $CPU
 fi
done
date
