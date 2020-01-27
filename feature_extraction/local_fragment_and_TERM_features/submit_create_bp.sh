#!/bin/bash

#SBATCH -n1
#SBATCH -N1
#SBATCH -p short
#SBATCH --mem=3G
#SBATCH -a1-8753

export OPENBLAS_NUM_THREADS=1
n=$SLURM_ARRAY_TASK_ID
pdb=$(sed "${n}q;d" input.lst ) ;
pdb_name=$( echo $(basename $pdb) | cut -d'.' -f1)
/home/basantab/scripts/make_bp.sh $(echo $pdb) raw_bps/$(echo $pdb_name).bp
sed -i '/^#/ d' raw_bps/$(echo $pdb_name).bp
sed -i '/^S/ d' raw_bps/$(echo $pdb_name).bp
