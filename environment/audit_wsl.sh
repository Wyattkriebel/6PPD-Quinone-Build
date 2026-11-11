#!/usr/bin/env bash
# Read-only: run inside the user's WSL Ubuntu, redirect to a new report file.
set -u
uname -a
cat /etc/os-release
lscpu
free -h
df -h .
command -v python python3 conda docker nvidia-smi nvcc || true
python3 --version
nvidia-smi || true
nvidia-smi --query-gpu=name,driver_version,memory.total,compute_cap --format=csv || true
nvcc --version || true
conda env list || true
# Detection alone is not a kernel-compatibility check.
