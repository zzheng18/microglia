#!/bin/bash
#SBATCH --partition=gpu_p
#SBATCH --qos=gpu_priority
#SBATCH --job-name=multigrate
#SBATCH --nodes=1
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=28
#SBATCH --mem=500G
#SBATCH --time=12:00:00
#SBATCH --output=/home/icb/zihe.zheng/projects/microglia/err/output_%A_%a.txt
#SBATCH --error=/home/icb/zihe.zheng/projects/microglia/err/err_%A_%a.txt
#SBATCH --nice=0
#SBATCH --array=0

source /home/icb/zihe.zheng/.bashrc

# conda activate multimil
# srun python /home/icb/zihe.zheng/projects/microglia/integration/multimil_all.py

# conda activate multigrate
# srun python /home/icb/zihe.zheng/projects/microglia/integration/multigrate_all_cell_type.py

# conda activate multimil2
# srun python /home/icb/zihe.zheng/projects/microglia/disease_prediction/integrated_prediction_ad.py

# conda activate sysvi
# srun python /home/icb/zihe.zheng/projects/microglia/integration/sysvi.py

conda activate sysvi
srun python /home/icb/zihe.zheng/projects/microglia/integration/scvi_all_cell_type_rna.py