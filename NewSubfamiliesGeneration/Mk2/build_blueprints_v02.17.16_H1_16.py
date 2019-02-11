""
#Script for writing blueprint files from a string defining the protein topology
#H[n-m]L[n-m]
#Length: xxx
""
# USE: python build_blueprints.v2.py -xml template_bb+design.xml -blueresfile
import itertools
from argparse import ArgumentParser
import re
import sys
import os
import copy
from Blueprint import Blueprint
import numpy as np
from BPBf_av_mod import *
import subprocess

#==============================
# INPUT PARAMETERS
#==============================
parser = ArgumentParser()
parser.add_argument('-xml', type=str, help="xml template")
parser.add_argument('-pdb', type=str, help="input pdb")
parser.add_argument('-wts', type=str, help="wts file for design")
parser.add_argument('-motif', type=str, help="motif file")
parser.add_argument('-blueresfile', action='store_true')
parser.add_argument('-abego', type=str)
parser.add_argument('-resfile', type=str)
args = parser.parse_args()

template_xml = args.xml
motif_filename = args.motif
abego = args.abego # for L3
#topol = "L[1-1]H[15-15]L[2-2]H[7-7]L[5-5]E[4-4]L[2-2]E[4-4]L[1-1]H[10-15]L[1-2]E[11-11]L[2-2]E[10-10]L[2-2]E[12-12]L[2-2]E[10-10]L[1-1]"
#topol = "L[1-1]H[18-18]L[1-2]H[6-6]L[5-5]E[4-4]L[2-2]E[4-4]L[1-1]H[11,13]L[2-3]E[11-11]L[2-2]E[10-10]L[2-2]E[12-12]L[2-2]E[10-10]L[1-1]"
topol = "L[1-1]H[18-18]L[2-2]H[7-7]L[5-5]E[4-4]L[2-2]E[4-4]L[1-1]H[14-14]L[2-2]E[12-12]L[2-2]E[12-12]L[2-2]E[14-14]L[2-2]E[9-9]L[1-1]" # Original H3 length: 14

if motif_filename:

	ss,combinations, dic_abego = MotifTopology(topol,motif_filename)
else:
	ss,combinations = GetCombinations(topol)
#print combinations

print topol
#==============================
# ncenter: Bulge N-ter is free of restrictions 

# ccenter: Bulge C-ter is free of restrictions



bulges =  {'E3':7,'E6':5} # Specific definition. Bulges are located at positions defined by the user


# Make directory and bluprints for each combination along with bulge positions
for comb in combinations:
	## Build directories and add bulge combinations
	#---------------------------------------------------------
	# Make the directory name
	filename = 'r0-all_rev-' ; strand=0 ; resnum=0
	for i,s in enumerate(ss):
		filename+='%s%i' %(ss[i],comb[i])		

	# First build a plain blueprint file to after creat a bp object to make easier the bulge combinations

	MakePlainBlueprint(ss,comb,'bp')
	blue = Blueprint('bp')

	# Bulges
	keys = bulges.keys() # bulged strand names
	keys.sort()
	bpos_dic = {} # all positions considered for each bulged strand
	for key in keys:
		for seg in blue.segments:
			if seg.id == key:
				st_len = len(seg.bp_data)
				if bulges[key] == 'n-center':		
					if st_len %2 == 0: # even strand length. Bulge must be at odd position
						bposs = range(1,st_len,2)[1:-1]
					else:
						bposs = range(2,st_len,2)[1:-1]
				elif bulges[key] == 'c-center': # it depends on the orientation of the first residue. In this case...
					bposs = range(1,st_len,2)[1:-1]
				else:
					bposs=[bulges[key]]

				bpos_dic[key] = bposs

	# Take all bulge combinations
	bcomb=[]
	for key in keys:
	        bcomb.append( bpos_dic[key] )

	bcombinations = list(itertools.product(*bcomb))

	for j,bulcomb in enumerate(bcombinations):
			picked_bulges={}
			pathname=filename
			for k,key in enumerate(keys):
				pathname += '-b%s.%i' %(key,bulcomb[k]) # this is the new filename
				picked_bulges[key]=bulcomb[k]

			# Make directory of the topology
			if not os.path.exists(pathname):
				os.mkdir(pathname)
			os.chdir(pathname)
			# copy files provided by the first step
			#os.system('cp ../input.pdb .')
			#os.system('cp ../bp1 .')

			## Build blueprints
			MakeRefBlueprint(ss,comb,picked_bulges,refblue = 'bp')			
			# split blueprint according to different phases in backbone generation
			# pairing according to current blueprint
			###########
			# ATTENTION!!! to make more efficient the biasing with the SSPAIR
			# headers of blueprints it is better only write the pairing that is 
			# actually built with this blueprint, then 100% of biasing is focused on the region being constructed.
			# Pairings must be indexed according to the current blueprint not the global one in contrast to segments
			#########
			# step1
			MakeFirstBlueprint(refblue = 'bp', segments = ['E4','L8','E5'], \
					   newblue = 'bp1', ss_pairing={'E1':['E2.A']},adapt={'E5':'E4'},specific_abego={'L2':[(0,'G'),(1,'G')]})
			
			# step2
			AddSegmentToBlueprint(refblue = 'bp', segments = ['E3','L7'], blue0 = 'bp1', newblue='bp2', \
						      append=False, ss_pairing_shift="SSPAIR 1-2.A.1", ss_pairing={'E1':['E2.A'],'E2':['E3.A']})
			# step3
			shift = Shift(refblue = 'bp', seg1 = 'E5', seg2 = 'E4')
			AddSegmentToBlueprint(refblue = 'bp', segments = ['L9','E6'], \
					      blue0 = 'bp2', newblue='bp3', append=True, \
					      ss_pairing={'E1':['E2.A'],'E2':['E3.A'],'E3':['E4.A']},insert={'E3':shift})

			# step4
			AddSegmentToBlueprint(refblue = 'bp', segments = ['H2','L3','E1','L4','E2','L5','H3','L6'], \
					      blue0 = 'bp3', newblue='bp4', append=False,\
					      ss_pairing_shift="SSPAIR 1-2.A.0;1-6.P.-5;3-4.A.99", \
					      hs_pairing='HSSTRIPLET 2,1-2;1,5-6',seg_abego={'L3':'BBAAB','L4':'EA','L5':'E','L6':'BA'}) # Added B to L6
			# step 5
			#AddSegmentToBlueprint(refblue = 'bp', segments = ['L10','H4'], blue0 = 'bp4', \
                        #                      newblue='bp5', append=True,  ss_pairing_shift="SSPAIR 5-6.A.99;3-4.A.99", \
                        #                      hs_pairing="HSSTRIPLET 4,3-4", hh_pairing="HHPAIR 2-4.A", seg_abego={'L10':'B'} )
			# step6
			AddSegmentToBlueprint(refblue = 'bp', segments = ['H1','L2'], blue0 = 'bp4', \
					      newblue='bp6', append=False,  ss_pairing_shift="SSPAIR 1-2.A.0;1-6.P.-5;3-4.A.99;4-5.A.99;5-6.A.99", \
					      hh_pairing="HHPAIR 1-2.A",hs_pairing='HSSTRIPLET 2,1-2',only_remodel=['H2','L3'], seg_abego={'L2':'GB'}) # seg_abego={'L2':'GB'}
			#---------------------------------------------------------
			# Blueprints are already written. Now XMLs, cst, condor files...
			# Write constraints for good HB pairing of built strands
			#for bp in ['bp1.b','bp2.b','bp3.b']:
			#	os.system('generate_sheetcst2.py %s > %s.cst' %(bp,bp))
			# input pdb
			if args.pdb:
				os.system('cp ../%s input.pdb' %(args.pdb))
			else:
				write_dummy_pdb('input.pdb')

			if args.wts:
				os.system('cp ../%s .' %(args.wts))
			# take the sheet pdb reset to one
			# RESFILE
			if args.resfile:
				os.system( 'python ../manual_resfile.py bp6.b %s' %( args.resfile ) )

			# XML
		    	os.system('cp ../%s foo.xml' %(template_xml))
			xml_lines = open('../%s' %(template_xml),'r').readlines()
   			# names of blueprints must be consistent between here and the template.xml
			# Move above dir for another topology

			################################################
			## Define bulged hairping shifts to find hbond partner
			blue = Blueprint('bp6') ; blue.reindex_blueprint(start=1) # Get the global blueprint.
			# First bulged hairpin: E3-E4. I define it as a "N-bulged hairpin"
			# For a N-bulged hairpin
                        bulged_strand = 'E3'
                        seg = blue.segment_dict[bulged_strand]
                        bulge1_shift = ( len(seg.bp_data) - bulges[bulged_strand] ) - 1

                        # Second bulged hairpin: E5-E6. I define it as a "C-bulged hairpin"
                        # For a C-bulged hairpin
                        bulged_strand = 'E6'
                        seg = blue.segment_dict[bulged_strand]
                        bulge2_shift = bulges[bulged_strand]  
			 
			#----------------------
			# Write flags file
			#----------------------
			flags_out = open('flags','w')
			################################################

			# step1
			blue = Blueprint('bp1') ; blue.reindex_blueprint(start=1)
			fileout = open('cst1','w') # E2//E3
			fileout_min = open('cst1_min','w') # E2//E3
                        #st11 = RegularStrandCurvature(level=3,blueprint=blue,strand='E1',\
                        #global_twist=40.0, global_twist_tol=30.0 )			
                        #st22 = RegularStrandCurvature(level=3,blueprint=blue,strand='E2',\
                        #global_twist=50.0, global_twist_tol=35.0 )

			
                        st1 = RegularStrandCurvature(level=2,blueprint=blue,strand='E1',\
                        global_bend=55.0,global_bend_tol=10.0,\
			global_twist=35.0,global_twist_tol=10.0,\
			constraint_type='bounded' )

                        st2 = RegularStrandCurvature(level=2,blueprint=blue,strand='E2',\
			global_twist=35.0,global_twist_tol=10.0,\
                        global_bend=55.0,global_bend_tol=10.0,constraint_type='bounded' )
			'''
			st1 = RegularStrandCurvature(level=2,blueprint=blue,strand='E1',\
                        global_bend=45.0,global_bend_tol=5.0,\
                        global_twist=55.0,global_twist_tol=5.0,\
                        constraint_type='bounded' )

                        st2 = RegularStrandCurvature(level=2,blueprint=blue,strand='E2',\
                        global_twist=55.0,global_twist_tol=5.0,\
                        global_bend=45.0,global_bend_tol=5.0,constraint_type='bounded' )
			# Original: 70 twist, tol 5.0
			'''
			# Vars for Stage #1
			SS_s1 = "-parser:script_vars s1_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp1"])) ; flags_out.write(SS_s1)
			bp1_name = "-parser:script_vars bp1=../bp1\n" ; flags_out.write(bp1_name)
			bp1b_name = "-parser:script_vars bp1.b=../bp1.b\n" ; flags_out.write(bp1b_name)
			bp1 = Blueprint('bp1'); r_res = [ i+1 for i,pos in enumerate(bp1.bp_data) if pos[-1] == 'R' ]
			R_in_bp1= "-parser:script_vars R_in_bp1=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp1)
			cst1_name = "-parser:script_vars cst1=../cst1\n" ; flags_out.write(cst1_name)
			cst1_min_name = "-parser:script_vars cst1_min=../cst1_min\n" ; flags_out.write(cst1_min_name)
			
			sthb = HbondsRegularHairpin(strand1="E1",strand2="E2",blueprint=blue) ; fileout.write(sthb)			

			fileout.write(st1) ; fileout.write(st2) # ; fileout.write(st11) ; fileout.write(st22)
			fileout_min.write(st1) ; fileout_min.write(st2)

			fileout.close()	

			#############################################

                        # step2
                        fileout = open('cst2','w') # E1

			# bulge constraints
			blue = Blueprint('bp2.b') ; blue.reindex_blueprint(start=1)
			s1 = blue.segment_dict['E1']	
			s2 = blue.segment_dict['E2']
			
			st = HbondsBulgedStrand(strand1="E1",strand2="E2",blueprint=blue) ; fileout.write(st)
			#st = BulgedStrandCurvature(strand="E1",bend2=50,bend2_tol=5.0,blueprint=blue,constraint_type='bounded') ; fileout.write(st)
			
                        # Specifiy movemap for minimization in XML
                        pos1 = s1.bp_data[0]
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='MoveMap2 Flexible',xxx=1, yyy=s2.bp_data[0][0]-1)
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='MoveMap2 Rigid',xxx=s2.bp_data[0][0], yyy=blue.bp_data[-1][0])
                        fileout.close()
			
			# Vars for Stage #2
			SS_s2 = "-parser:script_vars s2_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp2"])) ; flags_out.write(SS_s2)
			bp2_name = "-parser:script_vars bp2=../bp2\n" ; flags_out.write(bp2_name)
			bp2b_name = "-parser:script_vars bp2.b=../bp2.b\n" ; flags_out.write(bp2b_name)
                        bp2 = Blueprint('bp2'); r_res = [ i+1 for i,pos in enumerate(bp2.bp_data) if pos[-1] == 'R' ]
			R_in_bp2= "-parser:script_vars R_in_bp2=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp2)
			cst2_name = "-parser:script_vars cst2=../cst2\n" ; flags_out.write(cst2_name)

                        #############################################

                        # step3 # in this test i dont reead cst3 in xml.  try without these new constraints only for E4
                        fileout = open('cst3','w') # E4

                        # bulge constraints
                        blue = Blueprint('bp3.b') ; blue.reindex_blueprint(start=1)
                        s1 = blue.segment_dict['E1']
                        s2 = blue.segment_dict['E2']
                        s3 = blue.segment_dict['E3']
                        s4 = blue.segment_dict['E4']	

			# hbonds
			st = HbondsBulgedStrand(strand1="E3",strand2="E4",blueprint=blue) ; fileout.write(st)
			st = HbondsRegularHairpin(strand1="E2",strand2="E3",blueprint=blue) ; fileout.write(st)
		

			# csts of bulged E4
			#st = BulgedStrandCurvature(strand="E4",bend1=60.0,bend1_tol=20.0, blueprint=blue,constraint_type='bounded')  ; fileout.write(st)

			# cst of remodeled E3
                        #st1 = RegularStrandCurvature(level=2,blueprint=blue,strand='E3',\
                        #global_twist=50.0,global_twist_tol=35.0 )
			#fileout.write(st1)

			# Inter-strand twist: E2-E3
			dic_pairs = AllSheetSegmentPairs(blue)
			p1 = s3.bp_data[-1][0]
			p2 = p1-4
			p3 = dic_pairs[p2]['E2']
			p4 = p3-2
			twist = -35.0 ; twist_tol=10.0
			#st = CircularHarmonicCaDihedralConstraints(p1,p2,p3,p4,twist,twist_tol) ; fileout.write(st)
			#st = CaDihedralConstraints(p1,p2,p3,p4,twist,twist_tol) ; fileout.write(st)
			################ DIHEDRAL CST FOR BETTER TWIST IN E5 and E6 ###########
			## 20140415 changed residues and angle for E3 to better reproduce native twist
			cE3 = s3.bp_data[-1][0]
			bulge_list = Bulges(blue)
                        bulge2 = bulge_list[1]
			#st = HarmonicDihedralConstraints(cE3-4,cE3,82,10) ; fileout.write(st)
			#st = HarmonicDihedralConstraints(bulge2,bulge2+4,125,10) ; fileout.write(st)
			######################################################################

                        # Specifiy movemap for minimization in XML
                        pos1 = s1.bp_data[0]
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='MoveMap3 Rigid',xxx=blue.bp_data[0][0], yyy=s2.bp_data[-2][0])
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='MoveMap3 Flexible',xxx=s2.bp_data[-1][0], yyy=blue.bp_data[-1][0])

                        fileout.close()
			
			# Vars for Stage #3
			SS_s3 = "-parser:script_vars s3_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp3"])) ; flags_out.write(SS_s3)
			bp3_name = "-parser:script_vars bp3=../bp3\n" ; flags_out.write(bp3_name)
			bp3b_name = "-parser:script_vars bp3.b=../bp3.b\n" ; flags_out.write(bp3b_name)
			bp3 = Blueprint('bp3'); r_res = [ i+1 for i,pos in enumerate(bp3.bp_data) if pos[-1] == 'R' ]
			R_in_bp3= "-parser:script_vars R_in_bp3=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp3)
			cst3_name = "-parser:script_vars cst3=../cst3\n" ; flags_out.write(cst3_name)
			
			##### STEP 4 #####
                        blue = Blueprint('bp4')
			blue.reindex_blueprint(start=1)
                        bulge_list = Bulges(blue)
			bulge2 = bulge_list[1]
			e1 = blue.segment_dict['E1']
			e1_hb = e1.bp_data[0][0]
			seg1 = blue.segment_dict['H1']
			seg2 = blue.segment_dict['H2']
			h2c = int(seg1.bp_data[-1][0])
			h3n = int(seg2.bp_data[0][0])
			fileout = open('cst4','w')
			fileoutb = open('cst4b','w')
			fileoutc = open('cst4c','w')
			# 2 pairs of hydrogen bonding between E1 and E6
			st = CircularHBondConstraints(e1_hb,bulge2) ; fileout.write(st) ; fileoutb.write(st)
                        st = CircularHBondConstraints(bulge2+2,e1_hb) ; fileout.write(st); fileoutb.write(st) # ensure correct E1-E6 pairing

                        st = CircularHBondConstraints(e1_hb+2,bulge2+2) ; fileout.write(st); fileoutb.write(st) # ensure correct E1-E6 pairing
                        st = CircularHBondConstraints(bulge2+4,e1_hb+2) ; fileout.write(st); fileoutb.write(st) # ensure correct E1-E6 pairing

			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist4a',xxx=e1_hb, yyy=bulge2)

			# loop L3 constraints
			st = CircularHBondConstraints(bulge2-1,h2c) ; fileout.write(st); fileoutb.write(st)
			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist4b',xxx=bulge2-1, yyy=h2c)
			st = CircularHBondConstraints(h2c+2,bulge2-1) ; fileout.write(st); fileoutb.write(st)
			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist4c',xxx=h2c+2, yyy=bulge2-1)
			st = CircularHBondConstraints(h3n,h2c+3) ; fileout.write(st); fileoutb.write(st)
			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist4d',xxx=h3n, yyy=h2c+3) #deactivated in xml
			st = CircularHBondConstraints(h2c+5,h2c+2) ; fileout.write(st); fileoutb.write(st)
			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist4e',xxx=h2c+5, yyy=h2c+2)
			# loop L5 constraints
			st = CircularHBondConstraints(h3n+3,h3n-1) ; fileout.write(st); fileoutb.write(st)
			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist4f',xxx=h3n+3, yyy=h3n-1) #deactivated in xml
			st = CircularHBondConstraints(h3n-1,h2c+5) ; fileout.write(st); fileoutb.write(st)
			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist4g',xxx=h3n-1, yyy=h2c+5) #deactivated in xml
			## CSTs for H3-E3 loop:
			E4 = blue.segment_dict['E4']
                        e4c = int(E4.bp_data[-1][0])
			L6 = blue.segment_dict['L5']
                        l6c = int(L6.bp_data[-1][0])
			h3c = int(seg2.bp_data[-1][0])
			l6n = int(L6.bp_data[0][0])
			st = "AtomPair N %i O %i HARMONIC 2.5 0.3\n" %(l6n,h3c-2) ; fileout.write(st) ; fileoutb.write(st)
			st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(l6n,h3c-3); fileout.write(st) ; fileoutb.write(st)
			st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(l6n+1,e4c-1); fileout.write(st) ; fileoutb.write(st)
			st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(l6n+2,e4c-1); fileout.write(st) ; fileoutb.write(st)
			#st = "AtomPair N %i O %i HARMONIC 3.0 0.5\n" %(e4c,l6n); fileout.write(st) ; fileoutb.write(st) #Additional Hbond between H3 C term and L8 N term
			## H3-H4 distance from hairpin:     
                        e2 = blue.segment_dict['E2']
                        e2n = e2.bp_data[0][0]
                        l5 = blue.segment_dict['L5']
                        l5n = int(l5.bp_data[0][0])
			
			# E1-E2 Hairpin constraints
			st = HbondsRegularHairpin(strand1="E1",strand2="E2",blueprint=blue) ; fileout.write(st); fileoutb.write(st);
			
			# Perfect Helices
			st = PerfectHelixCst('bp4',1); fileoutb.write(st); fileout.write(st);
			st = PerfectHelixCst('bp4',2); fileoutb.write(st); fileout.write(st);
			
			# Angle between E4 C, L6_2 and H3C to avoid weird loop
			st = AngleConstraints(e4c-1,l6n+1,h3c,95,3,"E6H4") ; fileout.write(st); fileoutb.write(st);
			st = AngleConstraints(e4c-1,l6n+1,h3c,95,15,"E6H4b") ; fileoutc.write(st); fileoutb.write(st);
			
			#Enforce ABEGO at positions 10 and 9:
			l2 = blue.segment_dict['L2']
			l2N = int(l2.bp_data[0][0])
			st = AbegoPhiPsiConstraints(l2N,blue) ; fileoutb.write(st);
			st = AbegoPhiPsiConstraints(l2N+1,blue) ; fileoutb.write(st);

			fileout.close()
			fileoutb.close()
			fileoutc.close()
			
			# Vars for Stage #4
			SS_s4 = "-parser:script_vars s4_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp4"])) ; flags_out.write(SS_s4)
			#bp4_name = "-parser:script_vars bp4=../bp4\n" ; flags_out.write(bp4_name)
			bp4b_name = "-parser:script_vars bp4.b=../bp4.b\n" ; flags_out.write(bp4b_name)
			bp4 = Blueprint('bp4'); r_res = [ i+1 for i,pos in enumerate(bp4.bp_data) if pos[-1] == 'R' ]
                        R_in_bp4= "-parser:script_vars R_in_bp4=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp4)
			cst4_name = "-parser:script_vars cst4=../cst4\n" ; flags_out.write(cst4_name)
			cst4b_name = "-parser:script_vars cst4b=../cst4b\n" ; flags_out.write(cst4b_name)
			cst4c_name = "-parser:script_vars cst4c=../cst4c\n" ; flags_out.write(cst4c_name)
			H34dist_r1 = "-parser:script_vars H34dist_r1=%s\n"%(str(e2n)) ; flags_out.write(H34dist_r1)
			H34dist_r2 = "-parser:script_vars H34dist_r2=%s\n"%(str(l5n)) ; flags_out.write(H34dist_r2)
			H34dist_dist = "-parser:script_vars H34dist_dist=%s\n"%("18") ; flags_out.write(H34dist_dist)
			
			#### STEP 5: Adding the C term Helix:
                        #### Distance between bulge2 and helix2
			'''
			blue = Blueprint('bp5') ; blue.reindex_blueprint(start=1)
                        seg0 = blue.segment_dict['H3']
                        seg1 = blue.segment_dict['E6']
			seg2 = blue.segment_dict['E4']
			seg3 = blue.segment_dict['E6']
			seg4 = blue.segment_dict['L9']
			seg5 = blue.segment_dict['H2']
			L10 = int(seg4.bp_data[0][0])
			E6c = int(seg3.bp_data[-1][0])
			H4c = int(seg0.bp_data[-1][0])
			E4c = int(seg2.bp_data[-1][0])
			H3c = int(seg5.bp_data[-1][0])
                        bulge_list = Bulges(blue) ; bulge2 = bulge_list[1]; bulge1 = bulge_list[0]

                        fileout = open('cst5','w')
			fileoutb = open('cst5b','w')
			## H4 csts:
                        st = PerfectHelixCst('bp5',3); fileout.write(st);
                        st = PairConstraints(H4c,H3c,6,2,"H4H3") ; fileout.write(st);
			st = PairConstraints(H4c,H3c,6,2.5,"H4H3b") ; fileoutb.write(st);
			st = PairConstraints(H4c-1,E4c,6,2,"H4E4a") ; fileout.write(st);
			st = PairConstraints(H4c-1,E4c,6,2.5,"H4Eb") ; fileoutb.write(st); 
			st = PairConstraints(H4c-1,E4c+3,6,2,"H4E5") ; fileout.write(st);
			st = PairConstraints(H4c-1,E4c+3,6,2.5,"H4E5b") ; fileoutb.write(st);
			#st = AngleConstraints(E6c-2,L10,H4c,80,10,"E6H4") ; fileout.write(st); fileoutb.write(st); # Changed angle form 80 to 110
			## Hbonds of linker-bulge between E6 and H5
			
			E5 = blue.segment_dict['E5']
			E5n = int(E5.bp_data[0][0])
			st = CircularHBondConstraints(L10,E5n+5) ; fileout.write(st)
                        st = CircularHBondConstraints(L10+1,E5n+5) ; fileout.write(st)
			
			fileout.close()
			fileoutb.close()

			# Vars for Stage #5
			SS_s5 = "-parser:script_vars s5_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp5"])) ; flags_out.write(SS_s5)
                        #bp5_name = "-parser:script_vars bp5=../bp5\n" ; flags_out.write(bp5_name)
                        bp5b_name = "-parser:script_vars bp5.b=../bp5.b\n" ; flags_out.write(bp5b_name)
                        bp5 = Blueprint('bp5'); r_res = [ i+1 for i,pos in enumerate(bp5.bp_data) if pos[-1] == 'R' ]
                        R_in_bp5= "-parser:script_vars R_in_bp5=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp5)
                        cst5_name = "-parser:script_vars cst5=../cst5\n" ; flags_out.write(cst5_name)
			cst5b_name = "-parser:script_vars cst5b=../cst5b\n" ; flags_out.write(cst5b_name)
                        distH4Sheet_r1 = "-parser:script_vars distH4Sheet_r1=%s\n"%(str(H4c-1)) ; flags_out.write(distH4Sheet_r1)
			distH4Sheet_r2 = "-parser:script_vars distH4Sheet_r2=%s\n"%(str(E4c)) ; flags_out.write(distH4Sheet_r2)
			distH4Sheet_dist = "-parser:script_vars distH4Sheet_dist=8.5\n" ; flags_out.write(distH4Sheet_dist)
			'''
			
                        ##### STEP 6 #####
                        ## Distance between bulge2 and helix2
                        blue = Blueprint('bp6') ; blue.reindex_blueprint(start=1)
			seg0 = blue.segment_dict['H1']
                        seg1 = blue.segment_dict['H2']
                        seg2 = blue.segment_dict['H3']
                        h2c = int(seg1.bp_data[-1][0])
                        h3n = int(seg2.bp_data[0][0])
                        e1 = blue.segment_dict['E1']
                        e1_hb = e1.bp_data[0][0]

                        bulge_list = Bulges(blue) ; bulge2 = bulge_list[1]; bulge1 = bulge_list[0]

                        fileout = open('cst6','w')
			fileoutb = open('cst6b','w')
			fileoutc = open('cst6c','w')
                        # HB contraints of bulge2 and H2
                        st = CircularHBondConstraints(e1_hb,bulge2) ; fileout.write(st) ; fileoutb.write(st)
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist6a',xxx=e1_hb, yyy=bulge2)
                        st = CircularHBondConstraints(bulge2-1,h2c) ; fileout.write(st) ; fileoutb.write(st)
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist6b',xxx=bulge2-1, yyy=h2c)
                        st = CircularHBondConstraints(h2c+2,bulge2-1) ; fileout.write(st) ; fileoutb.write(st)
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist6c',xxx=h2c+2, yyy=bulge2-1)
                        st = CircularHBondConstraints(h3n,h2c+3) ; fileout.write(st) ; fileoutb.write(st)
                        XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist6d',xxx=h3n, yyy=h2c+3)
                        st = CircularHBondConstraints(h2c+5,h2c+2) ; fileout.write(st) ; fileoutb.write(st)
			XMLReplaceXXXYYY(xml_lines=xml_lines,identifier='dist6e',xxx=h2c+5, yyy=h2c+2)
			
			## Perfect helices:
			st = PerfectHelixCst('bp6',1); fileout.write(st); fileoutb.write(st);
                        st = PerfectHelixCst('bp6',2); fileout.write(st); fileoutb.write(st);
			
			## H1 csts:
			h2n = int(seg1.bp_data[0][0])
			h2c = int(seg1.bp_data[-1][0])
			h1n = int(seg0.bp_data[0][0])
			h1c = int(seg0.bp_data[-1][0])
                        e1 = blue.segment_dict['E1']
                        e1_hb = e1.bp_data[0][0]
			e3 = blue.segment_dict['E3']
			e3c = int(e3.bp_data[-1][0])
			loop6 = int(blue.segment_dict['L6'].bp_data[0][0])
			
			st = PairConstraints(h1c-2,loop6,8,3,"H1L6a") ; fileout.write(st); fileoutb.write(st); fileoutc.write(st);
			st = PairConstraints(h1c,loop6,8,3,"H1L6b") ; fileout.write(st); fileoutb.write(st);
			#st = PairConstraints(h1c-2,loop6,8,3,"H1L6c") ; fileout.write(st); fileoutb.write(st);
			st = PairConstraints(h1n,e3c,8,3,"H1E3n") ; fileout.write(st); fileoutb.write(st);
			st = PairConstraints(h1n,e3c+3,8,2,"H1E3") ; fileout.write(st); fileoutb.write(st); fileoutc.write(st);
			st = PairConstraints(h1c-1,loop6+2,7.7,1.5,"H1L6b") ; fileout.write(st); fileoutb.write(st); fileoutc.write(st);
			st = AngleConstraints(h1c-6,h1c+1,h2c,30,8,"E6H4") ; fileout.write(st); fileoutb.write(st); fileoutc.write(st);
			
			#H1-capping Hbond
			LOOP2 = blue.segment_dict['L2']
			l2n = int(LOOP2.bp_data[0][0])

			st = "AtomPair N %i O %i HARMONIC 2.5 0.3\n" %(l2n,h1c-2) ; fileout.write(st) ; fileoutb.write(st)
			
			#Enforce ABEGO at positions 27 and 28:
                        l3 = blue.segment_dict['L3']
                        l3N = int(l3.bp_data[0][0])
                        st = AbegoPhiPsiConstraints(l3N,blue) ; fileout.write(st) ; fileoutb.write(st);
                        st = AbegoPhiPsiConstraints(l3N+1,blue) ; fileout.write(st) ; fileoutb.write(st);

			## Done with constraints
                        fileout.close()
			fileoutb.close()
			fileoutc.close()
			
			# Vars for Stage #6
			SS_s6 = "-parser:script_vars s6_SS=%s\n" %(subprocess.check_output([ "bash", "/work/basantab/scripts/ss_from_bp.sh", "bp6"])) ; flags_out.write(SS_s6)
                        bp6b_name = "-parser:script_vars bp6.b=../bp6.b\n" ; flags_out.write(bp6b_name)
                        bp6 = Blueprint('bp6'); r_res = [ i+1 for i,pos in enumerate(bp6.bp_data) if pos[-1] == 'R' ]
                        R_in_bp6= "-parser:script_vars R_in_bp6=%d-%d\n"%(r_res[0],r_res[-1]) ; flags_out.write(R_in_bp6)
                        cst6_name = "-parser:script_vars cst6=../cst6\n" ; flags_out.write(cst6_name)
			cst6b_name = "-parser:script_vars cst6b=../cst6b\n" ; flags_out.write(cst6b_name)
			cst6c_name = "-parser:script_vars cst6c=../cst6c\n" ; flags_out.write(cst6c_name)
			H1Na_r1 = "-parser:script_vars H1Na_r1=%s\n"%(str(h1n)) ; flags_out.write(H1Na_r1)
			H1Na_r2 = "-parser:script_vars H1Na_r2=%s\n"%(str(e3c+3)) ; flags_out.write(H1Na_r2)
			H1Na_dist = "-parser:script_vars H1Na_dist=%s\n"%("13.5") ; flags_out.write(H1Na_dist)
			H1Nb_r1= "-parser:script_vars H1Nb_r1=%s\n"%(str(h1n)) ; flags_out.write(H1Nb_r1)
			H1Nb_r2= "-parser:script_vars H1Nb_r2=%s\n"%(str(e3c)) ; flags_out.write(H1Nb_r2)
			H1Nb_dist = "-parser:script_vars H1Nb_dist=%s\n"%("12.5") ; flags_out.write(H1Nb_dist)
			H2N_r1 = "-parser:script_vars H2N_r1=%s\n"%(str(h2n)) ; flags_out.write(H2N_r1)
			H2N_r2 = "-parser:script_vars H2N_r2=%s\n"%(str(h3n+8)) ; flags_out.write(H2N_r2)
			H2N_dist = "-parser:script_vars H2N_dist=%s\n"%("8") ; flags_out.write(H2N_dist)
			H1C_r1 = "-parser:script_vars H1C_r1=%s\n"%(str(h1c)) ; flags_out.write(H1C_r1)
			H1C_r2 = "-parser:script_vars H1C_r2=%s\n"%(str(loop6)) ;flags_out.write(H1C_r2)
			H1C_dist = "-parser:script_vars H1C_dist=%s\n"%("11.0") ;flags_out.write(H1C_dist)
			'''
			<AtomicDistance name="H1Na" residue1="%%H1Na_r1%%" atomname1="CA" residue2="%%H1Na_r2%%" atomname2="CA" distance="%%H1Na_dist%%" confidence="1" />
			<AtomicDistance name="H1Nb" residue1="%%H1Nb_r1%%" atomname1="CA" residue2="%%H1Nb_r2%%" atomname2="CA" distance="%%H1Nb_dist%%" confidence="1" />
			<AtomicDistance name="H2N" residue1="%%H2N_r1%%" atomname1="CA" residue2="%%H2N_r2%%" atomname2="CA" distance="%%H2N_dist%%" confidence="0" />
			<AtomicDistance name="H1C" residue1="%%H1C_r1%%" atomname1="CA" residue2="%%H1C_r2%%" atomname2="CA" distance="%%H1C_dist%%" confidence="1" />
			'''
			flags_out.close()
			#-----------------------
			# Write modified xml
			#-----------------------
			xml_out = open('input.xml','w')
			for line in xml_lines:
				xml_out.write(line)
			xml_out.close()
			
			os.chdir('../')


