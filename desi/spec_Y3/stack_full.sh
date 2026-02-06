#! /bin/bash -l

#SBATCH -A desi
#SBATCH -C cpu
#SBATCH --qos=regular
#SBATCH --time=02:00:00
#SBATCH --nodes=1
## SBATCH --ntasks-per-node=1
#SBATCH -o ../../Outputs_Perlmutter/slurm-%j.out # STDOUT
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=r.henryliu@berkeley.edu

# salloc --nodes 1 --qos interactive --time 04:00:00 --constraint cpu --account m3058

# source /global/common/software/desi/desi_environment.sh 23.1 # inherits it
# source /global/common/software/desi/users/adematti/cosmodesi_environment.sh main
# module unload desiutil
# module load desiutil/3.2.6
conda activate rhliu_tSZ

python -u dsigma_thumbstack_yaml.py -p ./configs/diskring_full.yaml
python -u dsigma_thumbstack_yaml.py -p ./configs/dsigma_full.yaml