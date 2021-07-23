# #!/bin/bash

# para
export CUDA_MPS_PIPE_DIRECTORY=/home/bigsnpr_medicine/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/home/bigsnpr_medicine/nvidia-log

# MODE1 = GPU, MODE2 =CPU
CORE=8
MODE=3
NUMA_CORE=`expr $CORE - 1`

# 1k 16 core ALL CPU
SNPS=14k_1k

SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) 2>&1 | tee $LOG_FILE
# (time -v Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) 2>&1 | tee $LOG_FILE
echo "job $SNPS done..."
echo "==========="
#=============

# 2k 16 core ALL CPU
SNPS=14k_2k

SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) 2>&1 | tee $LOG_FILE
echo "job $SNPS done..."
echo "==========="
#=============

# 4k 16 core ALL CPU
SNPS=14k_4k

SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) 2>&1 | tee $LOG_FILE
echo "job $SNPS done..."
echo "==========="
#=============

# 8k 16 core ALL CPU
SNPS=14k_8k

SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) 2>&1 | tee $LOG_FILE
echo "job $SNPS done..."
echo "==========="
#=============

# 16k 16 core ALL CPU
SNPS=14k_16k

SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) 2>&1 | tee $LOG_FILE
echo "job $SNPS done..."
echo "==========="
#=============

# 32k 16 core ALL CPU
SNPS=14k_32k

SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) 2>&1 | tee $LOG_FILE
echo "job $SNPS done..."
echo "==========="
#=============


#=============

# 64k 16 core ALL CPU
SNPS=14k_64k

SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v numactl -C 0-$NUMA_CORE Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) 2>&1 | tee $LOG_FILE
echo "job $SNPS done..."
echo "==========="
#=============