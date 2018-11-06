#!/bin/bash


#SBATCH --job-name="array_samtools"
#SBATCH -p backfill
#SBATCH --mail-type="ALL"


cd ~/LUSTRE_HOME_DIRECTORY/RNA_Macrophage/

for f in $(find . -name '*accepted_hits.bam*'); do
    cp "$f" "$(dirname "${f:1}" | tr "/" "_")_$(basename $f)"
done

mkdir samtools
mv *.bam samtools
cd samtools

ls *.bam | parallel "samtools view -bq 1 {} | samtools sort - {.}; samtools index {.}.bam"



