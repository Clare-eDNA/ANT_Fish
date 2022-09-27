#!/bin/bash -e
#SBATCH -A uoo03666
#SBATCH -J CutAdapt
#SBATCH --time 1:00:00
#SBATCH -N 1
#SBATCH -c 8
#SBATCH -n 1
#SBATCH --mem=24G
#SBATCH --partition=large
#SBATCH --mail-user=c.adams@otago.ac.nz
#SBATCH --output cutAdapt.%j.out # CHANGE each run
#SBATCH --error cutAdapt.%j.err # CHANGE each run

module load cutadapt/3.5-gimkl-2020a-Python-3.8.2

cd /nesi/nobackup/uoo03666/ASP_GBS/samplesTrimmed

for sample in *.trim.fq.gz
do
base=$(basename ${sample} .trim.fq.gz)
echo "${base}"
cutadapt -l 67 -o "${base}".fq.gz "${base}".trim.fq.gz -m 67 -j 8
done

