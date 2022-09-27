# ANT_Fish
Antarctic Fish "ANT" Antimora rostrata, Blue antimora

These scripts take you through some of the basic data processing for GBS data
Please note that you'll want to run FastQC and MultiQC before and after demultiplexing to get a good idea of the samples' quality. 
You may also want to optimize little m as well as Big M in the Stacks program, although here we only have a script for optimization of Big M. 

All of these files are tailored to be run under the uoo03666 project on NZ's NeSI server under Clare Adams' username (adacl879). 
Thus, to run these on your own personal computer or your own NeSI account, you will have to change these names.

The GBS visualization analyses are not quite complete, as they would benefit from a true LEA analysis or otherwise.

Tutorials that may be useful or have useful parts include: 

Stacks2:
https://catchenlab.life.illinois.edu/stacks/manual/ (manual)
 https://doi.org/10.1111/2041-210X.12775  (rationale behind choosing M variations)

SNP filtering (although this data may find minor allele count filtering harsh):
https://www.ddocent.com/filtering/

Population genomic analyses once you have a VCF file:
https://grunwaldlab.github.io/Population_Genetics_in_R/gbs_analysis.html
https://adegenet.r-forge.r-project.org/files/tutorial-basics.pdf
