#!/bin/sh

#SBATCH -J bigsnpr_single           # Job name
#SBATCH -o ./single_nodes_log/bigsnpr_single_test_1k.out         # Name of stdout output file (%j expands to jobId)
##SBATCH -t 02:00:00                # Run time (hh:mm:ss) - 0.5 hours
#SBATCH -p gp1d                   # partition
#SBATCH -A GOV109113               # iService Project id
#SBATCH -n 1                       # ask for n srun processes
#SBATCH --cpus-per-task 1
#SBATCH --ntasks-per-node 8
#SBATCH --gres=gpu:1             # Number of GPUs per node

pwd; hostname; date

SNPS_COL=1
SNPS_ROW=14
CORES=4
MODE=2
BED_FILE=/home/u8294235/bigsnpr_medicine/snp_fastImpute_modify/fake_snps/${SNPS_ROW}k_${SNPS_COL}k_0.05.bed
TEMP_FILE=/home/u8294235/bigsnpr_medicine/snp_fastImpute_modify/rds_temp/temp${SNPS_ROW}_${SNPS_COL}

echo "Running example Rmpi script. Using $SLURM_JOB_NUM_NODES nodes with $SLURM_NTASKS 
tasks, each with $SLURM_CPUS_PER_TASK cores."

# module purge
# module load miniconda3
# conda activate R_36
module load nvidia/cuda/10.0
module load cudnn/latest
module load gnu8/8.3.0
scontrol show hostname > hostlist

cmd="time Rscript --vanilla fastImpute_GPU_MPI_SingleNode.R $CORES $MODE $BED_FILE $TEMP_FILE"

Rscript --vanilla create_temp_files.R $BED_FILE $TEMP_FILE

time srun $cmd

Rscript --vanilla imputation_check.R $BED_FILE $TEMP_FILE
