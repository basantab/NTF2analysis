# This script will print out all RMSDs calculated after TMalign,
# and at the end present mean and stddev.
# All vs all
#import matplotlib.pyplot as plt
#import matplotlib.mlab as mlab
from sys import argv
from TMalign import *
import numpy as np
import multiprocessing as mp
import time
import Bio.PDB

def AllTmalign(in1,in2):
    out = TMalign(in1,in2)
    rmsd = get_RMSD(out)
    return rmsd

def AllTmalign_TMscore(in1,in2):
    out = TMalign(in1,in2)
    tmscore = get_tmscore(out)
    return tmscore

def Allscores(in1,in2):
    out = TMalign(in1,in2)
    alldict = get_ftmscore_and_all_mets(out)
    return alldict

def PDBio_superimpose(in1,in2):
    ''' Alignment by CA '''
    parser = Bio.PDB.PDBParser()
    str1 = parser.get_structure('in1', in1)
    str2 = parser.get_structure('in2', in2)
    ref_model = str1[0]['A']
    alt_model = str2[0]['A']
    ref_atoms = []
    alt_atoms = []
    for ref_res in ref_model:
        ref_atoms.append(ref_res['CA'])
    for alt_res in alt_model:
        alt_atoms.append(alt_res['CA'])
    super_imposer = Bio.PDB.Superimposer()
    super_imposer.set_atoms(ref_atoms, alt_atoms)
    #super_imposer.apply(alt_model.get_atoms())
    rmsd = super_imposer.rms
    return rmsd

def BatchProcess(n,tup_list):
    result = []
    for i in tup_list:
        result.append(AllTmalign(i[0],i[1]))
    print("Batch ready")
    return (n,result)

def BatchProcess_tmscore(n,tup_list):
    result = []
    for i in tup_list:
        result.append(1-AllTmalign_TMscore(i[0],i[1]))
    print("Batch ready")
    return (n,result)

def BatchProcessGetAll(n,tup_list):
    result = []
    for i in tup_list:
        result_dict = Allscores(i[0],i[1])
        #result.append( [ result_dict[k] for k in ['fTM','TM','TM1','RMSD'] ] )
        result.append( (result_dict,i[0],i[1]) )
    print("Batch ready")
    return (n,result)

def BatchProcess_PDBio_superimpose(n,tup_list):
    result = []
    for i in tup_list:
        result.append(PDBio_superimpose(i[0],i[1]))
    print("Batch ready")
    return (n,result)

def PrepareBatch(decoy_list,processors,output_dir):
    fd = open(output_dir+'/decoy_list.txt','w')
    job_tuples = []
    decoys = decoy_list
    for n in range(0,len(decoys)):
        fd.write(decoys[n]+'\n')
        for m in range(n+1,len(decoys)):
            job_tuples.append((decoys[n],decoys[m]))
    fd.close()
    
    jobs = []
    length_of_decoys = int(len(decoys))
    final = ((length_of_decoys*length_of_decoys)-length_of_decoys)/2
    per_batch = 0

    if final%processors != 0:
        per_batch = int(final/processors + 1)
    else:
        per_batch = int(final/processors)
    print("%d comparisons will be done per processor" %(per_batch))
    for i in range(0,processors):
        current = []
        if len(job_tuples) < per_batch:
            current = job_tuples
                #current.append(n)
        else:
                current = current + job_tuples[:per_batch]
                job_tuples = job_tuples[per_batch:]
                #for n in range(0,per_batch):
                #    current.append(job_tuples.pop(0))
        jobs.append(current)
    print("Done preparing jobs")
    return jobs

def Make_alignement_tmscore(decoy_list,processors,output_dir):
    RMSDs = []
    print("Comparing all to template and all vs all")
    final = float(((len(decoy_list)*len(decoy_list))-len(decoy_list))/2)
    print("%d processors will be used"%(processors))
    jobs = PrepareBatch(decoy_list,processors,output_dir)
    pool =  mp.Pool(processes=processors)
    print("jobs: %d"%len(jobs))
    partials = [pool.apply_async(BatchProcess, args=(n,i)) for n,i in enumerate(jobs)]
    while True:
        if False in [x.ready() for x in partials]:
            time.sleep(10)
        else:
            break
    partials = [ i[1] for i in sorted( [x.get() for x in partials], key=lambda n: n[0])]
    for i in partials:
        for n in i:
            RMSDs.append(n)
    print("Done comparing, writing data to allxall.txt")
    print("Mean: "+str(np.mean(RMSDs))+" StDev: "+str(np.std(RMSDs)))
    allvall_arr = np.array(RMSDs)
    allvall_arr.tofile(output_dir+'/allxall.txt',sep=" ")
    print("Done writing data, making plots")
    bins = np.linspace(0, 3, 31)
    #plt.figure(figsize=(10,8))
    weights2 = np.zeros_like(RMSDs)+1./len(RMSDs)
    #plt.hist(RMSDs,bins,alpha=0.5, label='Design RMSD to all desings',weights=weights2,histtype='stepfilled')
    #plt.legend(loc='upper right')
    #plt.xlabel(r'RMSD',fontsize=15)
    #plt.ylabel('Frequency',fontsize=15)
    #plt.savefig(output_dir+'RMSD_distributions.png')
    print("Done with all.")

def Make_alignement(decoy_list,processors,output_dir):
    RMSDs = []
    print("Comparing all to template and all vs all")
    final = float(((len(decoy_list)*len(decoy_list))-len(decoy_list))/2)
    print("%d processors will be used"%(processors))
    jobs = PrepareBatch(decoy_list,processors,output_dir)
    pool =  mp.Pool(processes=processors)
    print("jobs: %d"%len(jobs))
    partials = [pool.apply_async(BatchProcess, args=(n,i)) for n,i in enumerate(jobs)]
    while True:
        if False in [x.ready() for x in partials]:
            time.sleep(10)
        else:
            break
    partials_sorted = [ i[1] for i in sorted( [x.get() for x in partials], key=lambda n: n[0])]
    for i in partials_sorted:
        for n in i:
            RMSDs.append(n)
    print("Done comparing, writing data to allxall.txt")
    print("Mean: "+str(np.mean(RMSDs))+" StDev: "+str(np.std(RMSDs)))
    allvall_arr = np.array(RMSDs)
    allvall_arr.tofile(output_dir+'/allxall.txt',sep=" ")
    print("Done writing data, making plots")
    bins = np.linspace(0, 3, 31)
    #plt.figure(figsize=(10,8))
    weights2 = np.zeros_like(RMSDs)+1./len(RMSDs)
    #plt.hist(RMSDs,bins,alpha=0.5, label='Design RMSD to all desings',weights=weights2,histtype='stepfilled')
    #plt.legend(loc='upper right')
    #plt.xlabel(r'RMSD',fontsize=15)
    #plt.ylabel('Frequency',fontsize=15)
    #plt.savefig(output_dir+'RMSD_distributions.png')
    print("Done with all.")



def Make_alignement_quick(decoy_list,processors,output_dir):
    RMSDs = []
    print("Comparing all to template and all vs all")
    final = float(((len(decoy_list)*len(decoy_list))-len(decoy_list))/2)
    print("%d processors will be used"%(processors))
    jobs = PrepareBatch(decoy_list,processors,output_dir)
    pool =  mp.Pool(processes=processors)
    print("jobs: %d"%len(jobs))
    partials = [pool.apply_async(BatchProcess_PDBio_superimpose, args=(n,i)) for n,i in enumerate(jobs)]
    while True:
        if False in [x.ready() for x in partials]:
            time.sleep(10)
        else:
            break
    partials = [ i[1] for i in sorted( [x.get() for x in partials], key=lambda n: n[0])]
    for i in partials:
        for n in i:
            RMSDs.append(n)
    print("Done comparing, writing data to allxall.txt")
    print("Mean: "+str(np.mean(RMSDs))+" StDev: "+str(np.std(RMSDs)))
    allvall_arr = np.array(RMSDs)
    allvall_arr.tofile(output_dir+'/allxall.txt',sep=" ")
    print("Done writing data, making plots")
    bins = np.linspace(0, 3, 31)
    #plt.figure(figsize=(10,8))
    weights2 = np.zeros_like(RMSDs)+1./len(RMSDs)
    #plt.hist(RMSDs,bins,alpha=0.5, label='Design RMSD to all desings',weights=weights2,histtype='stepfilled')
    #plt.legend(loc='upper right')
    #plt.xlabel(r'RMSD',fontsize=15)
    #plt.ylabel('Frequency',fontsize=15)
    #plt.savefig(output_dir+'/RMSD_distributions.png')
    print("Done with all.")

def Make_alignment_all(decoy_list,processors,output_dir):
    RMSDs = []
    fTMs = []
    TMs = []
    TM1s = []
    print("Comparing all to template and all vs all")
    final = float(((len(decoy_list)*len(decoy_list))-len(decoy_list))/2)
    print("%d processors will be used"%(processors))
    jobs = PrepareBatch(decoy_list,processors,output_dir)
    pool =  mp.Pool(processes=processors)
    print("jobs: %d"%len(jobs))
    partials = [pool.apply_async(BatchProcessGetAll, args=(n,i)) for n,i in enumerate(jobs)]
    while True:
        if False in [x.ready() for x in partials]:
            time.sleep(10)
        else:
            break
    partials_sorted = [ i[1] for i in sorted( [x.get() for x in partials], key=lambda n: n[0])]
    out = open(output_dir+'/table_output_allxall.csv','w')
    out.write("mem1,mem2,fTM,TM,TM1,RMSD\n")
    for j in partials_sorted:
        for i in j:
           out.write("%s,%s,%0.3f,%0.3f,%0.3f,%0.3f\n"%(i[2],i[1],i[0]['fTM'],i[0]['TM'],i[0]['TM1'],i[0]['RMSD']))
           print("%s,%s,%0.3f,%0.3f,%0.3f,%0.3f"%(i[2],i[1],i[0]['fTM'],i[0]['TM'],i[0]['TM1'],i[0]['RMSD']))
        '''
        for n in i:
            RMSDs.append(n['RMSD'])
            fTMs.append(n['fTM'])
            TM1s.append(n['TM1'])
            TMs.append(n['TM'])
         '''
    print("Done comparing, writing data")
    out.close()
    '''
    allvall_arr_rmsd = np.array(RMSDs)
    allvall_arr_rmsd.tofile(output_dir+'/RMSDs_allxall.txt',sep=" ")
    allvall_arr_fTM = np.array(fTMs)
    allvall_arr_fTM.tofile(output_dir+'/fTMs_allxall.txt',sep=" ")
    allvall_arr_TM = np.array(TMs)
    allvall_arr_TM.tofile(output_dir+'/TMs_allxall.txt',sep=" ")
    allvall_arr_TM1 = np.array(TM1s)
    allvall_arr_TM1.tofile(output_dir+'/TM1s_allxall.txt',sep=" ")
    '''
    print("Done writing data, done with all!")

if __name__ == '__main__':
    decoys_lst_fname = argv[1]
    decoys_handle = open(decoys_lst_fname,'r')
    decoys = [ i[:-1] for i in decoys_handle.readlines() ]
    decoys_handle.close()
    test = Allscores(decoys[0],decoys[1])
    #Make_alignement_tmscore(decoys,20,'./')
    Make_alignment_all(decoys,20,'./')
