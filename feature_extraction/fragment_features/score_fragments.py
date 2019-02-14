# for i in `seq 0 n` ; do python score_fragments.py  $i & done (where n is the number of processors)
# Script written by Enrique Marcos for fragment picking and scoring.

import glob,os,re,sys
from numpy import average
from Blueprint import Blueprint
from math import ceil

def CheckSS(psipredfile, bluefile):
        blue = Blueprint(bluefile)
        blue.reindex_blueprint(start=1)
        filein = open(psipredfile,'r')
        psipred_lines = filein.readlines()

        cons=0 ; e1_cons=0
        for line in psipred_lines[2:]:
                ss = line.split()[2]
                index = int(line.split()[0])
                bluess = blue.bp_data[index-1][2][0]
                if bluess =='L':
                        bluess = 'C'
                if ss == bluess:
                        cons+=1

        nres = blue.bp_data[-1][0]
        conserv = cons/float(nres)*100. 

        return conserv





block = int(sys.argv[1])
designs = glob.glob('????')
nblocks = int(sys.argv[2])
s = int(ceil(len(designs)/float(nblocks)))
seldesigns = designs[block*s:(block+1)*s]
print 'Designing in these dirs ', seldesigns

fileout = open('score_tjfrags_%i.txt' %(block),'w')
for fragdir in seldesigns:
	        os.chdir(fragdir)
		# Get scores in pdb
		for line in open('00001.pdb'):
			if len(line.split())>0 and line.split()[0]=='score_res':
				score_res = float(line.split()[-1])
			if len(line.split())>0 and line.split()[0]=='holes':
				holes = float(line.split()[1])
	                if len(line.split())>0 and line.split()[0]=='pack':
        	                pack = float(line.split()[1])
                        if len(line.split())>0 and line.split()[0]=='sspred':
                                sspred = float(line.split()[1])
			if len(line.split())>0 and line[:4]=='ATOM':
				nres = int(line.split()[5])
                        if len(line.split())>0 and line.split()[0]=='uhb':
                                uhb = float(line.split()[1])
                        if len(line.split())>0 and line.split()[0]=='rotamer_boltz_core':
                                rotamer_boltz_core = float(line.split()[1])
                        if len(line.split())>0 and line.split()[0]=='ss_sc':
                                ss_sc = float(line.split()[1])				



		if not os.path.exists('frag_qual.dat'):
			os.system('/work/basantab/scripts/make_fragments_0614.py  -pdbs_folder . -n_proc 1')

		# fragment rmsd of 3-mers
		os.system('/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/Rosetta/main/source/bin/r_frag_quality.linuxgccrelease -in:file:native 00001.pdb -f  00001.200.3mers -out:qual frag_qual3.dat -ignore_unrecognized_res')
		frag_qual = map(float, open('good_fragments_count.txt','r').readlines()[0].split()[1:])
		nfragsok = frag_qual[1]/frag_qual[0]*100.

                os.system('/work/basantab/scripts/lrmsdwr.py frag_qual.dat > lrmsdwr') 
                lrmsdwr =  float( open('lrmsdwr','r').readlines()[0].split()[0] )

		# check SS
		bluefile = 'combined.bp'
		ss_conserv = 0.0#CheckSS('00001.psipred_ss2', bluefile)
		# Add Gabes' metrics:
		xs, ys=[],[]
		best = {}
		pdb = '%s/00001.pdb' %(fragdir)
		average_all_frags = []
		with open('frag_qual.dat') as file:
			lines = file.readlines()
			for line in lines:
				xs.append(int(line.split()[1]))
				ys.append(float(line.split()[3]))
				if xs[-1] not in best:
					best[xs[-1]] = []
				best[xs[-1]].append(ys[-1])
				average_all_frags.append(ys[-1])
			for i in best:
				best[i].sort()
			lowest = [best[x][0] for x in best]
			lowest.sort()
			lowest.reverse()
		with open('00001.fasta') as file:
			seq = file.readlines()[-1].strip()
		##
		## complress fragments
		os.system('gzip *mers')
		os.chdir('../') # come back to fragmens dir
		
		fileout.write('%s %6.3f %6.3f %6.3f %6.3f %6.3f%6.3f ' %(pdb,ss_conserv,nfragsok,frag_qual[2],frag_qual[3],frag_qual[3]/nres,lrmsdwr))
		fileout.write('%s ' % seq)
		fileout.write('%s ' % lowest[0])
		fileout.write('%s ' % lowest[1])
		fileout.write('%s ' % lowest[2])
		fileout.write('%s ' % lowest[3])
		fileout.write('%s ' % lowest[4])
		fileout.write('%s ' % lowest[5])
		fileout.write('sum6 ')
		fileout.write('%s ' % sum(lowest[0:6]))
		fileout.write('sum ')
		fileout.write('%s ' % sum(lowest))
		fileout.write('avBest ')
		fileout.write('%0.04f ' % (sum(lowest)/len(lowest)))
		fileout.write('%0.04f \n' % average(average_all_frags))
fileout.close()

#os.system('python submit_designlist_boinc.py %s-frag-%i.txt' %(dir,block))
