#%%
### coding for K mers
# feature of Kmers was followed from PreTIS
# there will be thoundsand of feature output to the file
# Py version: Python 3.11.4
import pandas as pd
import numpy as np
from tqdm import tqdm
import sys, os, warnings, argparse, time
warnings.filterwarnings("ignore")

amino_acid_and_codon_dict={
    "Alanine" : ["GCT","GCC","GCA","GCG"]
    ,"Leucine" : ["TTA","TTG","CTT","CTC","CTA","CTG"]
    ,"Arginine" : ["CGT","CGC","CGA","CGG","AGA","AGG"]
    ,"Lysine" : ["AAA","AAG"]
    ,"Asparagine" : ["AAT","AAC"]
    ,"Methionine" : ["ATG"]
    ,"Aspartic_acid" : ["GAT","GAC"]
    ,"Phenylalanine" : ["TTT","TTC"]
    ,"Cysteine" : ["TGT","TGC"]
    ,"Proline" : ["CCT","CCC","CCA","CCG"]
    ,"Glutamine" : ["CAA","CAG"]
    ,"Serine" : ["TCT","TCC","TCA","TCG","AGT","AGC"]
    ,"Glutamic_acid" : ["GAA","GAG"]
    ,"Threonine" : ["ACT","ACC","ACA","ACG"]
    ,"Glycine" : ["GGT","GGC","GGA","GGG"]
    ,"Tryptophan" : ["TGG"]
    ,"Histidine" : ["CAT","CAC"]
    ,"Tyrosine" : ["TAT","TAC"]
    ,"Isoleucine" : ["ATT","ATC","ATA"]
    ,"Valine" : ["GTT","GTC","GTA","GTG"]
    ,"Stop_codon" : ["TAA","TGA","TAG"]}
#    ,"Shine–Dalgarno": ["AGGAGG"]   ## 2024/07/03 SD need to be independently processed

from itertools import chain
all_type_of_codon_list = list(chain.from_iterable([*amino_acid_and_codon_dict.values()]))
amino_acid_and_codon_count_dict = {key: len(value) for key, value in amino_acid_and_codon_dict.items()}
#print(amino_acid_and_codon_count_dict)

TPTN_path = 'C:/Users/User/Desktop/course files archive/112/Summer Internship at AS/Kmer_redesign/Tomato_TP_and_TN_data.txt'
seq_PATH = 'C:/Users/User/Desktop/course files archive/112/Summer Internship at AS/Kmer_redesign/Tomato_seq.txt'
m = pd.read_csv(TPTN_path,sep="\t",index_col=0)
n = pd.read_csv(seq_PATH,sep="\t",header=None,names=["seq"],index_col=0)



def sudo_seq(gene_ID, n, ss, up, down):
    seq = n.loc[gene_ID].seq
    end = len(seq)
    #print(str(ss) + ':' + seq[ss:ss+3] + ':' + str(end))

    if ss < up:
        N_up = 'N'*(up-ss)
        seq = N_up + seq
    else:
        seq = seq[ss-up:]
    
    if (end - (ss+2)) < down:
        N_down = 'N'*(down-(end-(ss+1+2)))
        #print(len(N_down))
        seq = seq + N_down
    else:
        seq = seq[:(up+3+down)]
    
    #print(len(seq))
    #print(seq[99:102])
    return seq

def main():
    for i in range(m.shape[0]):
        if not i % 3203:
            tqdm.write('Complete calculations of '+str(i+1)+' genes.\n'+str(i/m.shape[0]*100)+'% were completed.')
        ss = m['py_position'][i]
        ID = m.index[i]
        #print('Processing sequence ' + ID)
        up = 99
        down = 99
        full_seq = sudo_seq(ID, n, ss, up = 99, down = 99)
        full_len = up + 3 + down
        
        feat_dict = {}
        # binary features: K_mers_position_plus/minus_pos_is_NT
        for pos in range(-up,0):
            for nt in 'ATCG': # Loop upstream area
                feat_name = 'K_mers_position_minus_' + str(abs(pos)) + '_is_' + nt
                if full_seq[up+pos] == nt:
                    match = 1
                else:
                    match = 0
                feat_dict[feat_name] = match
        for pos in range(4,down+4):
            for nt in 'ATCG': # Loop downstream area
                feat_name = 'K_mers_position_plus_' + str(pos) + '_is_' + nt
                if full_seq[up-1+pos] == nt:
                    match = 1
                else:
                    match = 0
                feat_dict[feat_name] = match
        
        # frequency features:
        ## codon:
        ### codon_all: K_mers_codon
        n_sites = full_len // 3
        all_tri = n_sites*3 - 2
        for codon in all_type_of_codon_list:
            feat_name = 'K_mers_' + codon
            freq = full_seq.count(codon)
            feat_dict[feat_name] = freq/all_tri
        ### codon_upstream: K_mers_upstream_codon
        n_sites = up // 3
        all_tri = n_sites*3 - 2
        for codon in all_type_of_codon_list:
            feat_name = 'K_mers_upstream_' + codon
            freq = full_seq[:up].count(codon)
            feat_dict[feat_name] = freq/all_tri
        ### codon_downstream: K_mers_downstream_codon
        n_sites = down // 3
        all_tri = n_sites*3 - 2
        for codon in all_type_of_codon_list:
            feat_name = 'K_mers_downstream_' + codon
            freq = full_seq[up+3:].count(codon)
            feat_dict[feat_name] = freq/all_tri
        ### in-frame_codon_downstream: K_mers_in_frame_downstream_codon
        n_sites = down // 3
        for codon in all_type_of_codon_list:
            feat_name = 'K_mers_in_frame_downstream_' + codon
            seq = full_seq[up+3:]
            site = 0
            count = 0
            for g in range(len(seq)):
                site = seq.find(codon)
                det = site % 3
                if site == -1:
                    break
                elif det == 0:
                    count += 1
                    seq = seq[site+3:]
                    #print(seq)
                elif det != 0:
                    seq = seq[site+(3-det):]
                    #print(seq)
            feat_dict[feat_name] = count/n_sites
        ### in-frame_codon_upstream: K_mers_in_frame_upstream_codon
        n_sites = up // 3
        for codon in all_type_of_codon_list:
            feat_name = 'K_mers_in_frame_downstream_' + codon
            seq = full_seq[:up]
            seq = seq[::-1] # Reverse the sequence
            site = 0
            count = 0
            for g in range(len(seq)):
                site = seq.find(codon)
                det = site % 3
                if site == -1:
                    break
                elif det == 0:
                    count += 1
                    seq = seq[site+3:]
                    #print(seq)
                elif det != 0:
                    seq = seq[site+(3-det):]
                    #print(seq)
            feat_dict[feat_name] = count/n_sites

        ## counts of NT:
        ### counts of NT in whole sequence: K_mers_nt
        for nt in 'ATCG':
            feat_name = 'K_mers_' + nt
            freq = full_seq.count(nt)
            feat_dict[feat_name] = freq/full_len
        ### counts of upstream NT: K_mers_upstream_nt
        for nt in 'ATCG':
            feat_name = 'K_mers_upstream_' + nt
            freq = full_seq[:up].count(nt)
            feat_dict[feat_name] = freq/up
        ### counts of downstream NT: K_mers_downstream_nt
        for nt in 'ATCG':
            feat_name = 'K_mers_downstream_' + nt
            freq = full_seq[up+3:].count(nt)
            feat_dict[feat_name] = freq/down
        
        ## Counts of amino acid
        ### Counts of amino acid in whole sequence: K_mers_aa
        n_sites = full_len // 3
        all_tri = n_sites*3 - 2
        for aa in amino_acid_and_codon_dict.keys():
            feat_name = 'K_mers_' + aa
            freq = 0
            for cd in amino_acid_and_codon_dict[aa]:
                freq += full_seq.count(cd)
            feat_dict[feat_name] = freq/all_tri
        ### AA_upstream: K_mers_upstream_aa
        n_sites = up // 3
        all_tri = n_sites*3 - 2
        for aa in amino_acid_and_codon_dict.keys():
            feat_name = 'K_mers_upstream_' + aa
            freq = 0
            for cd in amino_acid_and_codon_dict[aa]:
                freq += full_seq[:up].count(cd)
            feat_dict[feat_name] = freq/all_tri
        ### AA_downstream: K_mers_downstream_aa
        n_sites = down // 3
        all_tri = n_sites*3 - 2
        for aa in amino_acid_and_codon_dict.keys():
            feat_name = 'K_mers_downstream_' + aa
            freq = 0
            for cd in amino_acid_and_codon_dict[aa]:
                freq += full_seq[up+3:].count(cd)
            feat_dict[feat_name] = freq/all_tri
        ### in-frame_AA_downstream: K_mers_in_frame_downstream_aa
        n_sites = down // 3
        for aa in amino_acid_and_codon_dict.keys():
            feat_name = 'K_mers_in_frame_downstream_' + aa
            count = 0
            for cd in amino_acid_and_codon_dict[aa]:
                seq = full_seq[up+3:]
                site = 0
                for g in range(len(seq)):
                    site = seq.find(cd)
                    det = site % 3
                    if site == -1:
                        break
                    elif det == 0:
                        count += 1
                        seq = seq[site+3:]
                        #print(seq)
                    elif det != 0:
                        seq = seq[site+(3-det):]
                        #print(seq)
            feat_dict[feat_name] = count/n_sites
        ### in-frame_AA_upstream: K_mers_in_frame_upstream_aa
        n_sites = up // 3
        for aa in amino_acid_and_codon_count_dict.keys():
            feat_name = 'K_mers_in_frame_upstream_' + aa
            count = 0
            for cd in amino_acid_and_codon_dict[aa]:
                seq = full_seq[:up]
                seq = seq[::-1] # Reverse the sequence
                site = 0
                for g in range(len(seq)):
                    site = seq.find(cd)
                    det = site % 3
                    if site == -1:
                        break
                    elif det == 0:
                        count += 1
                        seq = seq[site+3:]
                        #print(seq)
                    elif det != 0:
                        seq = seq[site+(3-det):]
                        #print(seq)
            feat_dict[feat_name] = count/n_sites
        
        ## Counts of Shine-Dalgarno sequence
        ### Counts of SD in whole sequence: K_mers_Shine-Dalgarno
        n_sites = full_len // 6
        all_hex = n_sites*6 - 5
        feat_name = 'K_mers_Shine-Dalgarno'
        freq = full_seq.count('AGGAGG')
        feat_dict[feat_name] = freq/all_hex
        ### SD counts upstream: K_mers_upstream_Shine-Dalgarno
        n_sites = up // 6
        all_hex = n_sites*6 - 5
        feat_name = 'K_mers_upstream_Shine-Dalgarno'
        freq = full_seq[:up].count('AGGAGG')
        feat_dict[feat_name] = freq/all_hex
        ### SD counts downstream: K_mers_downstream_Shine-Dalgarno
        n_sites = down // 6
        all_hex = n_sites*6 - 5
        feat_name = 'K_mers_downstream_Shine-Dalgarno'
        freq = full_seq[up+3:].count('AGGAGG')
        feat_dict[feat_name] = freq/all_hex
        ### SD counts in-frame downstream: K_mers_in_frame_downstream_Shine-Dalgarno
        n_sites = down // 6
        feat_name = 'K_mers_in_frame_downstream_Shine-Dalgarno'
        seq = full_seq[up+3:]
        count = 0
        for SD in range(len(seq)):
            site = seq.find('AGGAGG')
            det = site % 6
            if site == -1:
                break
            elif det == 0:
                count += 1
                seq = seq[site + 6:]
            elif det != 0:
                seq = seq[site+(6-det):]
        feat_dict[feat_name] = count/n_sites
        ### SD counts in-frame upstream: K_mers_in_frame_upstream_Shine-Dalgarno
        n_sites = up // 6
        feat_name = 'K_mers_in_frame_upstream_Shine-Dalgarno'
        seq = full_seq[:up]
        seq = seq[::-1] # Reverse the sequence
        count = 0
        for SD in range(len(seq)):
            site = seq.find('AGGAGG')
            det = site % 6
            if site == -1:
                break
            elif det == 0:
                count += 1
                seq = seq[site + 6:]
            elif det != 0:
                seq = seq[site+(6-det):]
        feat_dict[feat_name] = count/n_sites

        
        # Merge rows of features into a dataframe
        if i == 0:
            feats = pd.DataFrame(data = feat_dict, index = [ID])
        else:
            feat_pd = pd.DataFrame(data = feat_dict, index = [ID])
            feats = [feats, feat_pd]
            feats = pd.concat(feats)
        #print(feats)
        #time.sleep(3)
    Kmers = m.reset_index().rename(columns={'index':'gene_ID'}).join(feats.reset_index(drop=True))
    Kmers.to_csv(path_or_buf= os.getcwd() + '/Kmers_IFtest.txt', sep='\t',index=False)

main()