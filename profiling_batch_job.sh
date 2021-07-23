RDS_TEMP=/home/bigsnpr_medicine/rds_temp

rm -rf $RDS_TEMP
echo "rds_temp deleted"
mkdir $RDS_TEMP
echo "rds_temp created"

echo "Profiling 14k * 1k ~ 32k in all CPU 16process"
sh fastImpute_CPU_serial_job.sh

echo "========  CPU DONE ========="

rm -rf $RDS_TEMP
echo "rds_temp deleted"
mkdir $RDS_TEMP
echo "rds_temp created"

echo "Profiling 14k * 1k ~ 32k in CPU 16process + GPU"
sh fastImpute_GPU_serial_job.sh

echo "========  CPU + GPU DONE ========="