# -*- coding:utf-8 -*-
import sys
import re
import os
import time
import multiprocessing

#Make a dictionary with all theoretical Kmers.
def theoretical_Kmers(K, Kmer, base):
     if K==0:
         r_Kmer = reverse_complement(Kmer)
         if r_Kmer in theoretical_Kmer_list:
             if Kmer<r_Kmer:
                 theoretical_Kmer_list.add(Kmer)
                 theoretical_Kmer_list.remove(r_Kmer)
             else:
                 pass
         else:
             theoretical_Kmer_list.add(Kmer)
         return
     for i in range(4):
         new_Kmer = Kmer+base[i]
         theoretical_Kmers(K-1, new_Kmer, base)

#Create directory at present workplace if this directory does not exist
def crtdir(work_station, filename):
    if not os.path.exists(work_station + '/' + filename):
        cmd = 'mkdir ' + work_station + '/' + filename
        os.system(cmd)

#Obtain Inverted complement string of a DNA sequence
def reverse_complement(s):
    basecomplement = {"A": "T", "T": "A", "G": "C", "C": "G"}
    letters = list(s)
    letters = [basecomplement[base] for base in letters]
    return ''.join(letters)[::-1]

#Inverted complement Kmer count as default
def Kmercount(freq, Nmer):
    if Nmer in freq:
        freq[Nmer] += 1
    elif reverse_complement(Nmer) in freq:
        freq[reverse_complement(Nmer)] += 1
    else:
        freq[Nmer] = 1
    return freq

#Calculate all Kmer's number
def length(freq):
    Kmer_num = 0
    for i in freq:
        Kmer_num += freq[i]
    return Kmer_num

#Designed for Kmer counting without subtraction of background
def Kmer_statistics0(seq_dict):
    freq = {}
    for i in seq_dict:
        my_seq = seq_dict[i].upper()
        for t in range(len(my_seq)-K+1):
            Kmer = my_seq[t:t+K]
            deviation = K-Kmer.count('A')-Kmer.count('T')\
            -Kmer.count('C')-Kmer.count('G')
            if deviation == 0:
                freq = Kmercount(freq, Kmer)
    return freq

#Designed for Kmer counting with subtraction of background
def Kmer_statistics1(seq_dict):
    freq, freq1, freq2 = {}, {}, {}
    for i in seq_dict:
        my_seq = seq_dict[i].upper()
        for t in range(len(my_seq)-K+1):
            Kmer = my_seq[t:t+K]
            deviation = K-Kmer.count('A')-Kmer.count('T')\
            -Kmer.count('C')-Kmer.count('G')
            if deviation == 0:
                freq = Kmercount(freq, Kmer)
            if t == 0:
                Kmer1 = Kmer[:K-1]
                if Kmer1.count('A')+Kmer1.count('T')+Kmer1.count('C')+Kmer1.count('G') == K-1:
                    freq1 = Kmercount(freq1, Kmer1)
                Kmer1 = Kmer[1:K]
                if Kmer1.count('A')+Kmer1.count('T')+Kmer1.count('C')+Kmer1.count('G') == K-1:
                    freq1 = Kmercount(freq1, Kmer1)
                Kmer2 = Kmer[:K-2]
                if Kmer2.count('A')+Kmer2.count('T')+Kmer2.count('C')+Kmer2.count('G') == K-2:
                    freq2 = Kmercount(freq2, Kmer2)
                Kmer2 = Kmer[1:K-1]
                if Kmer2.count('A')+Kmer2.count('T')+Kmer2.count('C')+Kmer2.count('G') == K-2:
                    freq2 = Kmercount(freq2, Kmer2)
                Kmer2 = Kmer[2:K]
                if Kmer2.count('A')+Kmer2.count('T')+Kmer2.count('C')+Kmer2.count('G') == K-2:
                    freq2 = Kmercount(freq2, Kmer2)
            else:
                Kmer1 = Kmer[1:K]
                Kmer2 = Kmer[2:K]
                if Kmer1.count('A')+Kmer1.count('T')+Kmer1.count('C')+Kmer1.count('G') == K-1:
                    freq1 = Kmercount(freq1, Kmer1)
                if Kmer2.count('A')+Kmer2.count('T')+Kmer2.count('C')+Kmer2.count('G') == K-2:
                    freq2 = Kmercount(freq2, Kmer2)
    return freq, freq1, freq2

def Kmer_statistics2(seq_dict):
    freq = {}
    for i in seq_dict:
        my_seq = seq_dict[i].upper()
        for t in range(len(my_seq)-K+1):
            Kmer = my_seq[t:t+K]
            Nmer = ''.join([Kmer[x] for x in loc])
            deviation = K-Kmer.count('A')-Kmer.count('T')\
            -Kmer.count('C')-Kmer.count('G')
            if deviation == 0:
                freq = Kmercount(freq, Nmer)
    return freq

#Designed for spaced word with substraction of background
def Kmer_statistics3(seq_dict):
    freq, freq1, freq2 = {}, {}, {}
    for i in seq_dict:
        my_seq = seq_dict[i].upper()
        for t in range(len(my_seq)-K+1):
            Kmer = my_seq[t:t+K]
            Nmer = ''.join([Kmer[x] for x in loc])
            deviation = K-Kmer.count('A')-Kmer.count('T')\
                        -Kmer.count('C')-Kmer.count('G')
            if deviation == 0:
                freq = Kmercount(freq, Nmer)
            if t == 0:
                Nmer1 = Nmer[:N-1]
                if Nmer1.count('A')+Nmer1.count('T')+Nmer1.count('C')+Nmer1.count('G') == N-1:
                    freq1 = Kmercount(freq1, Nmer1)
                Nmer1 = Nmer[1:N]
                if Nmer1.count('A')+Nmer1.count('T')+Nmer1.count('C')+Nmer1.count('G') == N-1:
                    freq1 = Kmercount(freq1, Nmer1)
                Nmer2 = Nmer[:N-2]
                if Nmer2.count('A')+Nmer2.count('T')+Nmer2.count('C')+Nmer2.count('G') == N-2:
                    freq2 = Kmercount(freq2, Nmer2)
                Nmer2 = Nmer[1:N-1]
                if Nmer2.count('A')+Nmer2.count('T')+Nmer2.count('C')+Nmer2.count('G') == N-2:
                    freq2 = Kmercount(freq2, Nmer2)
                Nmer2 = Nmer[2:N]
                if Nmer2.count('A')+Nmer2.count('T')+Nmer2.count('C')+Nmer2.count('G') == N-2:
                    freq2 = Kmercount(freq2, Nmer2)
            else:
                Nmer1 = Nmer[1:N]
                Nmer2 = Nmer[2:N]
                if Nmer1.count('A')+Nmer1.count('T')+Nmer1.count('C')+Nmer1.count('G') == N-1:
                    freq1 = Kmercount(freq1, Nmer1)
                if Nmer2.count('A')+Nmer2.count('T')+Nmer2.count('C')+Nmer2.count('G') == N-2:
                    freq2 = Kmercount(freq2, Nmer2)
    return freq, freq1, freq2

#Calculate the frequency of Kmers
def frequency_vector_calculator0(freq):
    prob, index_lst = {}, []
    sum = length(freq)
    for Kmer in theoretical_Kmer_list:
        ICKmer = reverse_complement(Kmer)
        index_forward = str2int(Kmer)
        index_lst.append(index_forward)
        if Kmer not in freq and ICKmer not in freq:
            prob[index_forward] = 0
        elif Kmer in freq:
            prob[index_forward] = freq[Kmer]/sum
        else:
            prob[index_forward] = freq[ICKmer]/sum
    index_lst.sort()
    return prob, index_lst

#Calculate the frequency of Kmers which is substracted by expected degenerated Kmer frequency
def frequency_vector_calculator1(freq_lst):
    freq, freq1, freq2 = freq_lst[0], freq_lst[1], freq_lst[2]
    prob, prob1, prob2, bgd, del_bgd = {}, {}, {}, {}, {}
    sum0, sum1, sum2 = length(freq), length(freq1), length(freq2)
    for key in freq:
        prob[key] = freq[key]/sum0
    for key in freq1:
        prob1[key] = freq1[key]/sum1
    for key in freq2:
        prob2[key] = freq2[key]/sum2
    index_lst = []
    for Nmer in theoretical_Kmer_list:
        Nmer1, Nmer2, Nmer3 = Nmer[1:], Nmer[:-1], Nmer[1:-1]
        ICNmer, ICNmer1, ICNmer2, ICNmer3 = reverse_complement(Nmer)\
            , reverse_complement(Nmer1), reverse_complement(Nmer2), reverse_complement(Nmer3)
        index_forward, index_backward = str2int(Nmer), str2int(ICNmer)
        if (Nmer1 not in prob1 and ICNmer1 not in prob1) or\
            (Nmer2 not in prob1 and ICNmer2 not in prob1) or\
            (Nmer3 not in prob2 and ICNmer3 not in prob2):
            continue
        elif Nmer not in prob and ICNmer not in prob:
            if index_forward < index_backward:
                index_lst.append(index_forward)
                del_bgd[index_forward] = -1
            else:
                index_lst.append(index_backward)
                del_bgd[index_backward] = -1
            continue
        else:
            Nmer1 = ICNmer1 if Nmer1 not in prob1 else Nmer1
            Nmer2 = ICNmer2 if Nmer2 not in prob1 else Nmer2
            Nmer3 = ICNmer3 if Nmer3 not in prob2 else Nmer3
            bgd[Nmer] = prob1[Nmer1]*prob1[Nmer2]/prob2[Nmer3]
        if index_forward < index_backward:
            index_lst.append(index_forward)
            if Nmer in prob:
                del_bgd[index_forward] = prob[Nmer]/bgd[Nmer]-1
            else:
                del_bgd[index_forward] = prob[ICNmer]/bgd[Nmer]-1
        else:
            index_lst.append(index_backward)
            if Nmer in prob:
                del_bgd[index_backward] = prob[Nmer]/bgd[Nmer]-1
            else:
                del_bgd[index_backward] = prob[ICNmer]/bgd[Nmer]-1
    index_lst.sort()
    return del_bgd, index_lst

#General statistics of Kmer number from all scaffolds
def integrate(sub_freq, integrated_freq):
    for Nmer in sub_freq:
        if Nmer in integrated_freq:
            integrated_freq[Nmer] += sub_freq[Nmer]
        elif reverse_complement(Nmer) in integrated_freq:
            integrated_freq[reverse_complement(Nmer)] += sub_freq[Nmer]
        else:
            integrated_freq[Nmer] = sub_freq[Nmer]
    return integrated_freq

#General statistics of Kmer number from all scaffolds
def summary0(sub_freq):
    freq = {}
    for item in sub_freq:
        freq = integrate(item, freq)
    return freq

#General statistics of Kmer number from all scaffolds,for degenerated Kmer background
def summary1(sub_freq_lst):
    freq, freq1, freq2 = {}, {}, {}
    for item in sub_freq_lst:
        freq = integrate(item[0], freq)
        freq1 = integrate(item[1], freq1)
        freq2 = integrate(item[2], freq2)
    return freq, freq1, freq2

def str2int(Kmer):
    nucleotide = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    Kmer_list = list(Kmer)
    location = len(Kmer)-1
    Kmer_value = 0
    for i in Kmer_list:
        Kmer_value += nucleotide[i]*4**location
        location -= 1
    return Kmer_value

def int2str(numeric_kmer):
    nucleotide = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}
    rst = []
    a = len(numeric_kmer)
    while a >= 0:
        rst.append(nucleotide[int(numeric_kmer/4**a)])
        numeric_kmer = numeric_kmer%(4**a)
        a = a-1
    Kmer = ''.join(rst)
    return Kmer


if __name__ == '__main__':
    help = "python3 Kmer_frequence_vector.py -d -N -fna.path -speciesname -CPU.NUM -bgd\n\
-d\t-num\tthe interval length of degenerated Kmer\td>0\n\
-N\t-num\tthe length of spaced word\tN>1\n\
-fna.path\t-str\tthe path of dna sequence as an input,whose suffix is fna,fa,ffn,or fasta\n\
-speciesname\t-str\tthe Latin name of species\n\
-CPU.NUM\t-num\tthe number of CPU be used\n\
-bgd\t-str\tthe type of background substracting\n\
    \t\t0 :\tno substraction of background (CV)\n\
    \t\t1 :\twith substraction of background (CV)\n\
    \t\t2 :\tno substraction of background (spaced word)\n\
    \t\t3 :\twith substraction of background (spaced word)\n"

    # Input Parameter:
    if len(sys.argv) != 7:
        sys.exit(help)
    d, N, filement, speciesname, coreNum, bgdType = \
    int(sys.argv[1]), int(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6]
    K = N+(N-1)*d
    loc = range(0, K, d+1)
    startTime = time.time()

    print('Pretreatment of genomes sequences in the form of scaffolds or chromosomes!!!')
    seq_dict = {}
    for row in open(filement, 'r'):
        if row[0] == '>':
            key = row[1:-1]
            seq_dict[key] = []
        else:
            seq_dict[key].append(row.strip())
    for key, value in seq_dict.items():
        seq_dict[key] = ''.join(value)
    #print(time.time()-startTime)
    #startTime = time.time()

    print("Scanning each subsequence with a k-string(contiguous length k word or\
degenerated word with a length k but less than k characters is count),\
count each k-string's number in each substring")
    queue = [{} for i in range(coreNum)]
    if bgdType == '0':
        i = 0
        for seqid in seq_dict:
            queue[i][seqid] = seq_dict[seqid]
            i += 1
            if i == coreNum:
                i = 0
        p = multiprocessing.Pool(coreNum)
        sub_freq = p.map(Kmer_statistics0, queue)
        p.close()
        p.join()
    elif d == 0:
        i = 0
        for seqid in seq_dict:
            queue[i][seqid] = seq_dict[seqid]
            i += 1
            if i == coreNum:
                i = 0
        p = multiprocessing.Pool(coreNum)
        freq_lst = p.map(Kmer_statistics1, queue)
        p.close()
        p.join()
    elif bgdType == '2':   
        i = 0
        for seqid in seq_dict:
            queue[i][seqid] = seq_dict[seqid]
            i += 1
            if i == coreNum:
                i = 0
        p = multiprocessing.Pool(coreNum)
        sub_freq = p.map(Kmer_statistics2, queue)
        p.close()
        p.join()
    elif bgdType == '3':   
        i = 0
        for seqid in seq_dict:
            queue[i][seqid] = seq_dict[seqid]
            i += 1
            if i == coreNum:
                i = 0
        p = multiprocessing.Pool(coreNum)
        freq_lst = p.map(Kmer_statistics3, queue)
        p.close()
        p.join()
    else:
        p.close()
        exit(help)
    #print(time.time()-startTime)
    #startTime = time.time()
    
    #Make a dictionary with all theoretical Kmers.
    base = ['A','C','G','T']
    theoretical_Kmer_list = {''}
    theoretical_Kmers(N, "", base)
    theoretical_Kmer_list.remove('')

    print("Make a summary statistic for all k-strings' number in all substrings")
    if bgdType == '0' or bgdType == '2':
        grass_freq = summary0(sub_freq)
    else:
        grass_freq_lst = summary1(freq_lst)
    #print(time.time()-startTime)
    #startTime = time.time()
    print("Calculate the frequence of k-strings")
    if bgdType == '0' or bgdType == '2':
        prob, index_lst = frequency_vector_calculator0(grass_freq)
    else:
        del_bgd, index_lst = frequency_vector_calculator1(grass_freq_lst)
    #print(time.time()-startTime)
    #startTime = time.time()

    # Output as a format of cvtree
    output = str(N)+'\n'
    if bgdType == '0' or bgdType == '2':
        countNoneZero, Inner_product = len(index_lst), 0
        for i in index_lst:
            Inner_product += prob[i]**2
    else:
        countNoneZero, Inner_product = len(index_lst), 0
        for i in index_lst:
            Inner_product += del_bgd[i]**2
    output += str(Inner_product)+'\n'+str(countNoneZero)+'\n'
    if bgdType == '0' or bgdType == '2':
        for i in index_lst:
            output += str(i)+' '+str(prob[i])+'\n'
    else:
        for i in index_lst:
            output += str(i)+' '+str(del_bgd[i])+'\n'
    bgdType = bgdType.split('-')[0]
    crtdir('./', 'Mode'+bgdType+'_K_'+str(K))
    open('Mode'+bgdType+'_K_'+str(K)+'/'+speciesname+'.cv.txt', 'w').write(output)
    print(time.time()-startTime)
