#!/software/miniconda3/envs/pyrosetta3/bin/python

import numpy as np
import pandas as pnd
import Blueprint

def get_segment_dict(bp_fname,for_fragments=False):
	'''
	Return a dictionary where keys are all the NTF2 subdomains:
	[ 'N-term_helix', 'H1H2_link', 'loop3_flank', 'hairpin', 'H3_n' , 'H3', 'H3C_str3', 'str3_4', 'str4_5', 'str_5_6', 'str6c_ch' ]
	and the values are the protein position indexes.
	'''
	bp = Blueprint.Blueprint(bp_fname)
	subunits = 1
	if len(bp.bp_data) > 200: subunits = 2
	container_dict = {}
	container_dict[ 'N-term_helix' ] = [ pos[0] for pos in bp.segment_dict['H1'].bp_data[:-4] ]
	container_dict[ 'H1H2_link' ] = [ pos[0] for pos in bp.segment_dict['H1'].bp_data[-4:]+bp.segment_dict['L2'].bp_data+bp.segment_dict['H2'].bp_data[:5] ]
	container_dict[ 'loop3_flank' ] = [ pos[0] for pos in bp.segment_dict['H2'].bp_data[5:]+bp.segment_dict['L3'].bp_data+bp.segment_dict['E1'].bp_data[:3] ]
	container_dict[ 'hairpin' ] = [ pos[0] for pos in bp.segment_dict['E1'].bp_data+bp.segment_dict['L4'].bp_data+bp.segment_dict['E2'].bp_data ]
	container_dict[ 'H3_n' ] = [ pos[0] for pos in bp.segment_dict['E2'].bp_data[-2:]+bp.segment_dict['L5'].bp_data+bp.segment_dict['H3'].bp_data[:5] ]
	container_dict[ 'H3' ] = [ pos[0] for pos in bp.segment_dict['H3'].bp_data ]
	container_dict[ 'H3C_str3' ] = [ pos[0] for pos in bp.segment_dict['H3'].bp_data[-4:]+bp.segment_dict['L6'].bp_data+bp.segment_dict['E3'].bp_data[:6] ]
	container_dict[ 'str3_4' ] = [ pos[0] for pos in bp.segment_dict['E3'].bp_data[6:]+bp.segment_dict['L7'].bp_data+bp.segment_dict['E4'].bp_data[:6] ]
	container_dict[ 'str4_5' ] = [ pos[0] for pos in bp.segment_dict['E4'].bp_data[6:]+bp.segment_dict['L8'].bp_data+bp.segment_dict['E5'].bp_data[:6] ]
	if for_fragments:
		str5_6_indexes = bp.segment_dict['E5'].bp_data[6:]+bp.segment_dict['L9'].bp_data+bp.segment_dict['E6'].bp_data[:6]
		container_dict[ 'str5_6' ] = [ pos[0] for pos in str5_6_indexes if pos[0] < int(len(bp.bp_data)/subunits)-9  ]
		if (len(bp.segment_dict['E6'].bp_data) > 6) and (bp.segment_dict['E6'].bp_data[6:][0][0] < int(len(bp.bp_data)/subunits)-9) :
			container_dict[ 'str6c_ch' ] = [ i for i in range( bp.segment_dict['E6'].bp_data[6:][0][0],int(len(bp.bp_data)/subunits)-9 ) ]
	else:
		container_dict[ 'str5_6' ] = [ pos[0] for pos in bp.segment_dict['E5'].bp_data[6:]+bp.segment_dict['L9'].bp_data+bp.segment_dict['E6'].bp_data[:6] ]
		if bp.segment_dict['E6'].bp_data[6:][0][0] < int(len(bp.bp_data)/subunits):
			container_dict[ 'str6c_ch' ] = [ i for i in range( bp.segment_dict['E6'].bp_data[6:][0][0],int(len(bp.bp_data)/subunits) ) ]
	return container_dict

def get_general_seq_features(bp_fname):
	bp = Blueprint.Blueprint(bp_fname)
	subunits = 1
	if len(bp.bp_data) > 200: subunits = 2
	try: test = [ i for i,pos in enumerate(bp.segment_dict['E6'].bp_data) if pos[2][1]=='A' ][0]
	except IndexError: return [ 'X' for i in range(48) ]
	def long_short_arm(bp):
		bulge2_relative_pos = [ i for i,pos in enumerate(bp.segment_dict['E6'].bp_data) if pos[2][1]=='A' ][0]
		return (bulge2_relative_pos>2)

	H1_H2 = bp.segment_dict['H1'].bp_data[-8:]+bp.segment_dict['L2'].bp_data+bp.segment_dict['H2'].bp_data+bp.segment_dict['L3'].bp_data
	#hairpin = bp.segment_dict['E1'].bp_data+bp.segment_dict['L4'].bp_data+bp.segment_dict['E2'].bp_data+bp.segment_dict['L5'].bp_data+bp.segment_dict['H3'].bp_data[:9]
	hairpin = bp.segment_dict['E1'].bp_data[:4]+bp.segment_dict['L4'].bp_data+bp.segment_dict['E2'].bp_data[-4:]+bp.segment_dict['L5'].bp_data+bp.segment_dict['H3'].bp_data[:9]
	hairpin_index = [  i for i in range(len(hairpin)) ]
	E6_pos = []
	E5_pos = []
	E4_pos = []
	if long_short_arm(bp):
		E6_pos = [2,3,4,5,6,7,8]
		#E5_pos = [-3,-4,-5,-6,-7,-8,-9,-10]
		E5_pos = [-5,-6,-7,-8,-9,-10,-11,-12]
		E4_pos = [2,3,4,5,6,7,8,9]
	else:
		E6_pos = [0,1,2,3,4,5,6]
		#E5_pos = [-1,-2,-3,-4,-5,-6,-7,-8]
		E5_pos = [-3,-4,-5,-6,-7,-8,-9,-10]
		E4_pos = [0,1,2,3,4,5,6,7]
	# Now get all the positions in one string:
	try: identitites = [ H1_H2[pos][1] for pos in range( len(H1_H2) ) ] + \
		      [ hairpin[pos][1] for pos in hairpin_index ] + \
		      [ bp.segment_dict['E4'].bp_data[pos][1] for pos in E4_pos ] + \
		      [ bp.segment_dict['E5'].bp_data[pos][1] for pos in E5_pos ] + \
		      [ bp.segment_dict['E6'].bp_data[pos][1] for pos in E6_pos ]
	except IndexError: return [ 'X' for i in range(65) ]
	if len(identitites) == 65:
		return identitites
	else:
		return [ 'X' for i in range(65) ]

def get_frag_metrics(folder,bp_fname):
	
	frag_qual_handle = open('/home/basantab/NTF2_project/20170707_filtering_chip_designs/fragments2/%s/frag_qual.dat'%folder,'r')
	frag_qual_lines = [ i[:-1] for i in frag_qual_handle.readlines() ]
	frag_qual_handle.close()
	rms_pos_dict = { int( line.split()[1] ):[] for line in frag_qual_lines }
	for line in frag_qual_lines: rms_pos_dict[int( line.split()[1] )].append(float(line.split()[3]))
	segments = get_segment_dict(bp_fname,for_fragments=True)
	metrics = {}
	for key in segments.keys():
		metrics[key] = { 'av_allfr':np.average( [ np.average(rms_pos_dict[pos]) for pos in segments[key] ] ) }
		metrics[key]['av_bestfr'] = np.average( [ min(rms_pos_dict[pos]) for pos in segments[key] ] )
		metrics[key]['av_worstfr'] = np.average( [ max(rms_pos_dict[pos]) for pos in segments[key] ] )
		metrics[key]['best_at_worstfr'] = max( [ min(rms_pos_dict[pos]) for pos in segments[key] ] )
	if 'str6c_ch' not in segments.keys():
		metrics['str6c_ch'] = { 'av_allfr':0, 'av_bestfr':0, 'av_worstfr':0, 'best_at_worstfr':0 }
	renamed_metrics = {}
	for key in metrics.keys():
		for feature in metrics[key].keys():
			new_name = key + '_' + feature
			renamed_metrics[new_name] = metrics[key][feature]
	return renamed_metrics

def get_term_metrics(terms_dict,bp_fname):
	
	segments = get_segment_dict(bp_fname,for_fragments=False)
	metrics = {}
	for key in segments.keys():
		for term in terms_dict.keys():
			metrics[key] = { term+'_av':np.average( [ terms_dict[term][i-1] for i in segments[key] ] ) }
			metrics[key][term+'_worst'] = max( [ terms_dict[term][i-1] for i in segments[key] ] )
	if 'str6c_ch' not in segments.keys():
		metrics['str6c_ch'] = { 'abd50_av':0, 'dsc50_av':0, 'ssc50_av':0, 'abd50_worst':0, 'dsc50_worst':0, 'ssc50_worst':0 }
	renamed_metrics = {}
	for key in metrics.keys():
		for feature in metrics[key].keys():
			new_name = key + '_' + feature
			renamed_metrics[new_name] = metrics[key][feature]
	return renamed_metrics


bp_handle = open('/home/basantab/NTF2_project/20180917_MiSeqData/sequence_features/all_processed_bps.list','r')
#folders_handle = open('/home/basantab/NTF2_project/20170707_filtering_chip_designs/fragments2/all_scores.txt','r')
#desname_handle = open('/home/basantab/NTF2_project/20170707_filtering_chip_designs/fragments2/mapped_all_fragment_metrics.list','r')
#TERMANAL_files = open('/home/basantab/NTF2_project/20170707_filtering_chip_designs/local_and_sequence_features/TERMANAL_file.lst','r')

bp_lines = [ i[:-1] for i in bp_handle.readlines() ]
#folders_lines = [ i[:-1] for i in folders_handle.readlines() ]
#desname_lines = [ i[:-1] for i in desname_handle.readlines() ]
#TERMANAL_lines = [ i[:-1] for i in TERMANAL_files.readlines() ]

bp_handle.close()
#folders_handle.close()
#desname_handle.close()
#TERMANAL_files.close()

# Create a dictionary from design name to folder:

designs = [ i.split()[0] for i in bp_lines ]
#folders = [ i.split()[0].split('/')[0] for i in folders_lines ]
#design_folder_dict = { des:folder for des,folder in zip(designs,folders) }

# Create TERMANAL feature dict

def get_scores_from_file(fname):
	with open( fname ) as f:
		vals = [ float(line.strip().split()[-1]) for line in f ]
	return vals

#TERMANAL_fname_dict = { fname.split('/')[-1].split('.')[0]:{ 'dsc':fname, 'abd':fname.split('.')[0]+'.abd50', 'ssc':fname.split('.')[0]+'.ssc50' } for fname in TERMANAL_lines if '.dsc' in fname }

#TERMANAL_dict = { des:{ 'abd50':get_scores_from_file(TERMANAL_fname_dict[des]['abd']),\
#			'dsc50':get_scores_from_file(TERMANAL_fname_dict[des]['dsc']),\
#			'ssc50':get_scores_from_file(TERMANAL_fname_dict[des]['ssc'])\
#			} for des in designs if des in TERMANAL_fname_dict.keys() }

# Create a dictionary from design to SS:

design_bp_dict = { line.split()[0]:line.split()[1] for line in bp_lines }

# Now get the fragment quality at each segment:
keys = [des for des in design_bp_dict.keys()]

#term_metrics_storage_dict = { des:get_term_metrics(TERMANAL_dict[des],design_bp_dict[des]) for des in TERMANAL_dict.keys() }
#df_tm = pnd.DataFrame.from_dict(term_metrics_storage_dict,orient='index')
#df_tm.to_csv('output_tm_test1.csv')
sequence_met = { des:get_general_seq_features(design_bp_dict[des]) for des in keys }
sequence_met_noX = { key:value for key,value in sequence_met.items() if 'X' not in value }
df_sm = pnd.DataFrame.from_dict(sequence_met_noX,orient='index')
df_sm.to_csv('output_seq_met_extended.csv')


