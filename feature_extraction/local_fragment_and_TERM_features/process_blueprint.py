#!/software/miniconda3/envs/pyrosetta3/bin/python
import Blueprint
import sys

ss_fname = '/net/scratch/basantab/NTF2_project/20190621_filter_and_select_longNTF2/TERMS/TERMANAL_digs/ss.tab'
ss_handle = open(ss_fname,'r')
lines = [ i[:-1] for i in ss_handle.readlines()]
ss_handle.close()
ss_dict = { i.split()[0].split('.')[0]:i.split()[1] for i in lines }

bp_fnames = sys.argv[1:]

for bp_fname in bp_fnames:
	key_name = bp_fname.split('/')[-1].split('.')[0]
	bp = Blueprint.Blueprint(bp_fname)
	ss = ss_dict[ key_name ]
	new_bp_vec = [ [ pos[0], pos[1], s+pos[2][1], '.' ] for s,pos in zip(ss,bp.bp_data) ]
	new_bp = Blueprint.Blueprint(data=new_bp_vec)
	key_name_bp = key_name + '.bp'
	new_bp.dump_blueprint('./processed_bps/%s'%key_name_bp)
