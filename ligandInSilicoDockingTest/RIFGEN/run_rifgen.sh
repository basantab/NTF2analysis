#!/bin/bash

export OMP_NUM_THREADS=32

cd $1;
donres=$(echo THR TYR TRP GLN ASN SER LYS ARG);
accres=$(echo THR TYR SER GLN ASN ASP GLU);
params=$(ls ???.params)

if [ -f has_don.txt ]; then donres=$(echo $donres $1) ; fi
if [ -f has_acc.txt ]; then accres=$(echo $accres $1) ; fi

/home/bcov/rifdock/latest/rifgen -rifgen:target $(ls ???_rot_min_0001.pdb) -extra_res_fa $params -rifgen:outdir ./rifgen_init -rifgen:data_cache_dir ./rifgen_init -rifgen:outfile rifgen_init.rif.gz -rifgen:apores PHE ILE LEU VAL MET -rifgen:donres $donres -rifgen:accres $accres @./template_flags &> output.log

cd ..
