#!/software/miniconda3/envs/pyrosetta3/bin/python

import numpy as np
import pandas as pnd

def ProcessSS(input_ss):
	'''
	Return a dictionary where keys are H, E and L, and values are lists of indexes where elements start and end.
	E.g.: 'E':[ [1,7], [14,19] ] there are two strands, one starts at index 1 and ends at 7, and the other one starts at 14 and
	ends at 19.
	'''
	ss_string = list(input_ss)
	container_dict = {'E':[],'H':[],'L':[]}
	switch = False
	index_container = [0,0]
	for i,n_pos in enumerate(ss_string):
		prev_index = max(i-1,0)
		prev = ss_string[prev_index]
		if prev != n_pos: switch = True
		else: switch = False
		if switch:
			container_dict[prev].append(index_container)
			index_container = [i,i]
		else: index_container = [ index_container[0],i ]
		if ( i == len(ss_string) - 1 ):
			index_container = [ index_container[0],i ]
			container_dict[n_pos].append(index_container) 
	return container_dict

def get_des_metrics(folder,ss):
	frag_qual_handle = open('%s/frag_qual.dat'%folder,'r')
	frag_qual_lines = [ i[:-1] for i in frag_qual_handle.readlines() ]
	frag_qual_handle.close()
	rms_pos_dict = { int( line.split()[1] ):[] for line in frag_qual_lines }
	for line in frag_qual_lines: rms_pos_dict[int( line.split()[1] )].append(float(line.split()[3]))
	indexes = ProcessSS(ss)
	# Loops:
	n_term = 0
	c_term = len(ss)-1
	# remove n and c term loops: 
	indexes['L']  = [ loop for loop in indexes['L'] if (n_term not in loop) and (c_term not in loop) ]
	# extract loop numbers:
	loop_pos_and_flank =  [ [ i+1 for i in range(loop[0]-1,loop[1]+2) if i+1 in rms_pos_dict.keys() ] for loop in indexes['L'] ]
	av_all_loop = np.average( [ np.average( [ np.average(rms_pos_dict[pos]) for pos in loop ] ) for loop in loop_pos_and_flank ] )
	av_best_loop = np.average( [ np.average([ min(rms_pos_dict[pos]) for pos in loop ]) for loop in loop_pos_and_flank ] )
	max_av_loop = max( [ np.average([ np.average(rms_pos_dict[pos]) for pos in loop ]) for loop in loop_pos_and_flank ] )
	max_av_best_loop = max( [ np.average([ min(rms_pos_dict[pos]) for pos in loop ]) for loop in loop_pos_and_flank ] )
	point_loop_av_all = np.average( [ np.average(rms_pos_dict[loop[0]]) for loop in loop_pos_and_flank ] )
	point_loop_av_worst = max( [ np.average(rms_pos_dict[loop[0]]) for loop in loop_pos_and_flank ] )
	# extract strand numbers:
	strand_pos_and_flank =	[ [ i+1 for i in range(strand[0]-1,strand[1]+2) if i+1 in rms_pos_dict.keys() ] for strand in indexes['E'] ]
	av_all_strand = np.average( [ np.average( [ np.average(rms_pos_dict[pos]) for pos in strand ] ) for strand in strand_pos_and_flank ] )
	av_best_strand = np.average( [ np.average([ min(rms_pos_dict[pos]) for pos in strand ]) for strand in strand_pos_and_flank ] )
	max_av_strand = max( [ np.average([ np.average(rms_pos_dict[pos]) for pos in strand ]) for strand in strand_pos_and_flank ] )
	max_av_best_strand = max( [ np.average([ min(rms_pos_dict[pos]) for pos in strand ]) for strand in strand_pos_and_flank ] )
	# extract helix numbers:
	helix_pos_and_flank =  [ [ i+1 for i in range(helix[0]-1,helix[1]+2) if i+1 in rms_pos_dict.keys() ] for helix in indexes['H'] ]
	av_all_helix = np.average( [ np.average( [ np.average(rms_pos_dict[pos]) for pos in helix ] ) for helix in helix_pos_and_flank ] )
	av_best_helix = np.average( [ np.average([ min(rms_pos_dict[pos]) for pos in helix ]) for helix in helix_pos_and_flank ] )
	max_av_helix = max( [ np.average([ np.average(rms_pos_dict[pos]) for pos in helix ]) for helix in helix_pos_and_flank ] )
	max_av_best_helix = max( [ np.average([ min(rms_pos_dict[pos]) for pos in helix ]) for helix in helix_pos_and_flank ] )

	metrics = { 'av_all_loop':av_all_loop,\
		    'av_best_loop':av_best_loop,\
		    'max_av_loop':max_av_loop,\
		    'max_av_best_loop':max_av_best_loop,\
		    'point_loop_av_all':point_loop_av_all,\
		    'point_loop_av_worst':point_loop_av_worst,\
		    'av_all_strand':av_all_strand,\
		    'av_best_strand':av_best_strand,\
		    'max_av_strand':max_av_strand,\
		    'max_av_best_strand':max_av_best_strand,\
		    'av_all_helix':av_all_helix,\
		    'av_best_helix':av_best_helix,\
		    'max_av_helix':max_av_helix,\
		    'max_av_best_helix':max_av_best_helix }
	
	return metrics

ss_handle = open('./NTF2_chip2.dssp','r')
folders_handle = open('./all_fragment_metrics_correct.tab','r')
desname_handle = open('./mapped_all_fragment_metrics_correct.list','r')

ss_lines = [ i[:-1] for i in ss_handle.readlines() ]
folders_lines = [ i[:-1] for i in folders_handle.readlines() ]
desname_lines = [ i[:-1] for i in desname_handle.readlines() ]

ss_handle.close()
folders_handle.close()
desname_handle.close()

# Create a dictionary from design name to folder:

designs = [ i.split()[0] for i in desname_lines ]
folders = [ i.split()[0].split('/')[0] for i in folders_lines ]
design_folder_dict = { des:folder for des,folder in zip(designs,folders) }
print( "Done creating dictionary from design name to folder")
# Create a dictionary from design to SS:

design_ss_dict = { line.split()[0]:line.split()[1] for line in ss_lines }
print( "Done creating dictionary from design name to SS")
# Now get the fragment quality at each ss element:
keys = [des for des in design_folder_dict.keys()]
des_metrics_storage_dict = { des:get_des_metrics(design_folder_dict[des],design_ss_dict[des]) for des in keys }

keys = sorted([ i for i in des_metrics_storage_dict[designs[0]].keys()])

output = open('SS_fragment_metrics.csv','w')
output.write(','+','.join(keys)+'\n')
for des in designs :
	metrics = des_metrics_storage_dict[des]
	str_ready = [ '%0.3f'%(metrics[key]) for key in keys]
	output.write(des+','+','.join(str_ready)+'\n')
	print(des+','+','.join(str_ready))
output.close()
