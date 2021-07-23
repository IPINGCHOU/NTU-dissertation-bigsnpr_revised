# 1k 16 core ALL CPU
SNPS=14k_1k
CORE=16
MODE=1
SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) > $LOG_FILE 2>&1
echo "job $SNPS done..."
echo "==========="
#=============

# 2k 16 core ALL CPU
SNPS=14k_2k
CORE=16
MODE=1
SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) > $LOG_FILE 2>&1
echo "job $SNPS done..."
echo "==========="
#=============

# 4k 16 core ALL CPU
SNPS=14k_4k
CORE=16
MODE=1
SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) > $LOG_FILE 2>&1
echo "job $SNPS done..."
echo "==========="
#=============

# 8k 16 core ALL CPU
SNPS=14k_8k
CORE=16
MODE=1
SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) > $LOG_FILE 2>&1
echo "job $SNPS done..."
echo "==========="
#=============

# 16k 16 core ALL CPU
SNPS=14k_16k
CORE=16
MODE=1
SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) > $LOG_FILE 2>&1
echo "job $SNPS done..."
echo "==========="
#=============

# 32k 16 core ALL CPU
SNPS=14k_32k
CORE=16
MODE=1
SNPS_ROUTE=/home/bigsnpr_medicine/fake_snps/twcc_fake_snps/${SNPS}_0.05.bed
RDS_ROUTE=/home/bigsnpr_medicine/rds_temp/$SNPS
LOG_FILE=/home/bigsnpr_medicine/fastImpute_serial_log/${SNPS}_core${CORE}_mode${MODE}.txt

echo "job $SNPS with core $CORE mode $MODE start..."
(time -v Rscript --vanilla fastImpute_GPU_serial.R $CORE $MODE $SNPS_ROUTE $RDS_ROUTE) > $LOG_FILE 2>&1
echo "job $SNPS done..."
echo "==========="
#=============