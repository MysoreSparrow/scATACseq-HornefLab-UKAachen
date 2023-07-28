
#!/bin/bash

#Request cores. In this case - 2 cores will be requested
#SBATCH --cpus-per-task=6

#Request memory per thread. Here it's 1 GB per requested core
#SBATCH --mem-per-cpu=256G

# Set a job name
#SBATCH --job-name=Run_scATAcseq_AGHornef

#Specify time before the job is killed by scheduler (in case is hangs). In this case - 24 hours
#SBATCH --time=72:00:00 

#Declare the merged STDOUT/STDERR file
#SBATCH --output=/data/keshav/sc-ATACseq/output/atac.%J.out

#Run something
#module load R
#echo $PATH #only diplays the current PATH-variable

#Please note: you have to load the modules and every other variable you might have set in yout .bashrc again in the execution file, in order for it to be loaded properly on the compute nodes.
pwd
conda info --envs
source /home/keshav/anaconda3/bin/activate
conda activate base

echo "Slurm Starting now!!!."

# Call the first script for preprocessing and filtering
Rscript	1_PreProcessing_DoubletFiltering.R
Rscript	2_Integration_Clustering_PeakCalling.R
Rscript	3_DownstreamAnalysis.R

echo "End of Slurm Run!"
