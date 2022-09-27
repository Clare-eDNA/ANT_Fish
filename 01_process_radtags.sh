#!/bin/bash -e

#SBATCH --job-name			radtags_inline_null			# the name of your job
#SBATCH --output			radtags1.%j.out			# Include the job ID in the names of
#SBATCH --error				radtags1.%j.err			# the output and error files
#SBATCH --time				20:00:00			# 60 min, this is the MAX time your job will run, (HH:MM:SS)
#SBATCH --mem				64GB				# Memory request #U need to specify for DADA2
#SBATCH --cpus-per-task			4
#SBATCH --export			NONE				# This will stop opening Unix Gui system X11
#SBATCH --chdir				/nesi/project/uoo03666/Clare/ASP_GBS/  # your work directory
#SBATCH --account			uoo03666

# load your program 
module purge
module load Stacks/2.61-gimkl-2022a

echo "loaded Stacks, starting run"

process_radtags -p /nesi/project/uoo03666/Clare/ASP_GBS/RAW/ -o /nesi/project/uoo03666/Clare/ASP_GBS/Plate1/ -b /nesi/project/uoo03666/Clare/ASP_GBS/barcodes/Plate_1A.txt -e pstI -r -c -q -P --threads 4 -D --inline-inline --filter-illumina

echo "finished plate one, starting on plate two"

process_radtags -p /nesi/project/uoo03666/Clare/ASP_GBS/RAW/ -o /nesi/project/uoo03666/Clare/ASP_GBS/Plate2/ -b /nesi/project/uoo03666/Clare/ASP_GBS/barcodes/Plate_2A.txt -e pstI -r -c -q -P --threads 4 -D --inline-inline --filter-illumina

echo "done"
