#!/bin/bash
# the following must be performed with root privilege
export CUDA_VISIBLE_DEVICES=0
export CUDA_MPS_PIPE_DIRECTORY=/home/bigsnpr_medicine/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/home/bigsnpr_medicine/nvidia-log

echo "step1 done"
nvidia-smi -i 0 -c EXCLUSIVE_PROCESS
echo "step2 done"
nvidia-cuda-mps-control -d
echo "all done"
