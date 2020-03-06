#!/bin/bash

cd $1
export OMP_NUM_THREADS=32
/home/bcov/rifdock/latest/rif_dock_test -scaffold_res $(cat /home/basantab/NTF2_project/20191216_deNovoNTF2DockingBenchmark_RIFDOCK/RIFGEN/POSFILE.list.$(printf %02d $2)) -scaffolds $(cat /home/basantab/NTF2_project/20191216_deNovoNTF2DockingBenchmark_RIFDOCK/RIFGEN/SCAFFOLD.list.$(printf %02d $2)) -rif_dock:dokfile all.dok.$2 @local_rifdock.flags
cd ..

