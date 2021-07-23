# #!/bin/zsh

CORE=4
MODE=1
NUMA_CORE=`expr $CORE - 1`
export CUDA_MPS_PIPE_DIRECTORY=/home/bigsnpr_medicine/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/home/bigsnpr_medicine/nvidia-log

# 1k 4 core CPU with GPU
SNPS=14k_1k
SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

# NVVP_FILENAME=${SNPS}_4cores_${CORES}process_%p
# NVVP_ROUTE=/home/bigsnpr_medicine/nvprof_result/${NVVP_FILENAME}.nvvp



# time -v numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial_SingleChr.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE

time numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE


# nvprof -f --print-gpu-trace --profile-child-processes --export-profile $NVVP_ROUTE $cmd
# nvprof -f --print-gpu-trace --profile-all-processes --export-profile $NVVP_ROUTE 



# nsys profile -o ./nvprof_result/4processes_nomps ./fastImpute_GPU_serial_job_nvprof.sh