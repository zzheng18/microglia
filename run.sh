#!/bin/bash
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_long
#SBATCH --job-name=multimil
#SBATCH --nodes=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=20
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --output=/home/icb/zihe.zheng/projects/microglia/err/output_%A_%a.txt
#SBATCH --error=/home/icb/zihe.zheng/projects/microglia/err/err_%A_%a.txt
#SBATCH --nice=1000
#SBATCH --array=0

source /home/icb/zihe.zheng/.bashrc

conda activate multimil
srun python /home/icb/zihe.zheng/projects/microglia/integration/multimil_all.py

# conda activate scvi
# srun python /home/icb/zihe.zheng/projects/microglia/peakvi.py
# srun python /home/icb/zihe.zheng/projects/microglia/scvi_rna.py