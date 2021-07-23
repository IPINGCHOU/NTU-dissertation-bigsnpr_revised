#!/bin/bash
export CUDA_VISIBLE_DEVICES=0
export CUDA_MPS_PIPE_DIRECTORY=/home/bigsnpr_medicine/nvidia-mps
export CUDA_MPS_LOG_DIRECTORY=/home/bigsnpr_medicine/nvidia-log
echo quit | nvidia-cuda-mps-control
nvidia-smi -i 0 -c DEFAULT
