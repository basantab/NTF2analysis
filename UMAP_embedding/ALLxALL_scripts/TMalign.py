from subprocess import check_output
import numpy as np

def get_RMSD(tmout):
    lines = tmout.splitlines()
    rmsd_line = lines[16]
    rmsd = rmsd_line.split()[4].split(',')[0]
    return float(rmsd)

def get_tmscore(tmout):
    lines = tmout.splitlines()
    tmscore1_line = lines[17]
    tmscore2_line = lines[18]
    tmscore1 = float(tmscore1_line.split()[1])
    tmscore2 = float(tmscore2_line.split()[1])
    tmscore = ( tmscore1 + tmscore2 )/2
    return tmscore

def get_tmscore_ch1(tmout):
    lines = tmout.splitlines()
    tmscore1_line = lines[17]
    tmscore1 = float(tmscore1_line.split()[1])
    return tmscore1

def get_ftmscore_and_all_mets(tmout):
    """
    As described in "Evolutionary inaccuracy of pairwise structural alignments" See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3338010/#SEC2
    As used in "Consequences of domain insertion on sequence-structure divergence in a superfold" See https://www.pnas.org/content/110/36/E3381
    """
    lines = tmout.splitlines()
    rmsd_line = lines[16]
    rmsd = float(rmsd_line.split()[4].split(',')[0])
    chain1_line = lines[13]
    chain2_line = lines[14]
    chain1_len = float(chain1_line.split()[3])
    chain2_len = float(chain2_line.split()[3])
    avlen = (chain1_len + chain2_len) / 2
    minlen = min(chain1_len,chain2_len)
    al_len = float(rmsd_line.split()[2].split(',')[0])
    d0 = (1.24 * np.cbrt( minlen - 15 )) - 1.8
    fTM = al_len / ( avlen * ( 1 + ( rmsd / d0 ) * ( rmsd / d0 ) ) )
    tmscore1_line = lines[17]
    tmscore1 = float(tmscore1_line.split()[1])
    tmscore2_line = lines[18]
    tmscore2 = float(tmscore2_line.split()[1])
    tmscore = ( tmscore1 + tmscore2 )/2
    return { 'fTM':fTM, 'TM':tmscore, 'TM1':max(tmscore1,tmscore2), 'RMSD':rmsd }

def get_aligned_pos(tmout,print_result=True):
    lines = tmout.splitlines()
    top_PDB_seq_algn = lines[-4][:-1]
    charact = lines[-3][:-1]
    bottom_PDB_seq_algn = lines[-2][:-1]
    print(top_PDB_seq_algn)
    print(charact)
    print(bottom_PDB_seq_algn)
    aligment_tuple_list = [ i for i in zip(top_PDB_seq_algn,charact,bottom_PDB_seq_algn)]
    nat_counter = 1
    des_counter = 1
    pos_align_tuple_list = []
    for tup in aligment_tuple_list:
        this_new_tup = ['-',tup[1],'-']
        if tup[0] != '-':
            this_new_tup[0] = nat_counter
            nat_counter += 1
        if tup[2] != '-':
            this_new_tup[2] = des_counter
            des_counter += 1
        pos_align_tuple_list.append(tuple(this_new_tup))
    if print_result:
        nat_counter = 1
        des_counter = 1
        pos_align_tuple_list_str = []
        for tup in aligment_tuple_list:
            this_new_tup_str = ['| - ',tup[1],'| - ']
            if tup[0] != '-':
               this_new_tup_str[0] = "|%03d"%nat_counter
               nat_counter += 1
            if tup[2] != '-':
               this_new_tup_str[2] = "|%03d"%des_counter
               des_counter += 1
            pos_align_tuple_list_str.append(tuple(this_new_tup_str))
        print( ''.join([ i[0] for i in pos_align_tuple_list_str ] )+'|' )
        print( ''.join([ '| '+i[1]+' ' for i in pos_align_tuple_list_str ] )+'|' )
        print( ''.join([ i[2] for i in pos_align_tuple_list_str ] )+'|' )
    return pos_align_tuple_list

def TMalign(str1,str2,**kwargs):
    to_exec = "/home/basantab/TMalign/TMalign"
    if 'exe' in kwargs: # get the TMalign executable ult=Trueath from optional args
        to_exec = kwargs['exe']
    TM_out = check_output([to_exec, str1, str2]).decode('ascii')
    return TM_out
        
if __name__ == '__main__':
    from os import path
    from sys import argv
    # test with:
    # python TMalign.py ~/NTF2_project/20151125_3x4_len1_H3x2_seq/only_freq_ABEGOs/hyak_run_20151208/20151208_sequence_on_short_big_pocket/j0/0/1573_L1H18L2H7L5E4L2E4L1H12L1H8L1E11L2E10L2E12L2E10L1_41_4_input_0006_1filtered_0001.pdb ~/NTF2_project/20151125_3x4_len1_H3x2_seq/only_freq_ABEGOs/hyak_run_20151208/20151208_sequence_on_short_big_pocket/j0/0/1573_L1H18L2H7L5E4L2E4L1H12L1H8L1E11L2E10L2E12L2E10L1_41_4_input_0006_1filtered_0001_0015.pdb ~basantab/bin/TMalign
    str1 = argv[1]
    str2 = argv[2]
    path = argv[3]
    out = TMalign(str1,str2)
    out = TMalign(str1,str2,exe=path)
    print(get_RMSD(out))
'''
**************************************************************************
 *                        TM-align (Version 20170708)                     *
 * An algorithm for protein structure alignment and comparison            *
 * Based on statistics:                                                   *
 *       0.0 < TM-score < 0.30, random structural similarity              *
 *       0.5 < TM-score < 1.00, in about the same fold                    *
 * Reference: Y Zhang and J Skolnick, Nucl Acids Res 33, 2302-9 (2005)    *
 * Please email your comments and suggestions to: zhng@umich.edu          *
 **************************************************************************

Name of Chain_1: /home/basantab/NTF2_project/20190710_all_vs_denovo
Name of Chain_2: /home/basantab/NTF2_project/20190710_all_vs_denovo
Length of Chain_1:  160 residues
Length of Chain_2:  144 residues

Aligned length=  104, RMSD=   2.98, Seq_ID=n_identical/n_aligned= 0.106
TM-score= 0.52062 (if normalized by length of Chain_1)
TM-score= 0.56914 (if normalized by length of Chain_2)
(You should use TM-score normalized by length of the reference protein)

(":" denotes aligned residue pairs of d < 5.0 A, "." denotes other aligned residues)
IEQPRWASKDPAAGK-----ASTPDEKIVLEFMDAL-T----SND--AAKLIKYFAEDTMYQNMPLPPAYGRDAVEQTLA---GLFKVM---SIDAVEVFHIGSS-KGLVFTERVDVLRALPTGKSYNLSILGVFQLTD----------------GKITGWRDYFDLREFEEAVDLPLRGKLGPEQKLISEEDLN-----
                    .::::::::::::::: :    :::  ::::::::::::            ::::::::.   ::::::   ::::::::::::: :::::::::::::....::::::::::::::::                ::::::::::.
---------------HETSYGEVVDRYWLNQYVLNRETYDYDTIQLNYDTTALLSTAAV------------QQEFYKIYEGEDARDKVLSNKARITVKVRSIQPNGRGQATVRFTTQQHDSTGAVGVKQHQIATIGYTYVGAPMKSSDRLLNPLGFQVTSYRTDPE-----------------------------ILLNN
'''

