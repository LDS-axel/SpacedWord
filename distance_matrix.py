# encoding:utf-8
import sys
import os
import time
import itertools
import math
import multiprocessing

#Create a directory if not exist at present workplace
def crtdir(work_station, filename):
    if not os.path.exists(work_station + '/' + filename):
        cmd = 'mkdir ' + work_station + '/' + filename
        os.system(cmd)

#Calculate the cosine similarity between two species' Kmer frequency vector
def cosine_dissimilariry(pairsSpecies):
    sp1, sp2 = pairsSpecies[0], pairsSpecies[1]
    dict1, dict2, index1, index2 = {}, {}, set(), set()
    sp_cv1 = open(directory+sp1+'.cv.txt', 'r').readlines()
    sp_cv2 = open(directory+sp2+'.cv.txt', 'r').readlines()
    inner_product_for_sp1 = float(sp_cv1[1].strip())
    inner_product_for_sp2 = float(sp_cv2[1].strip())
    for i in sp_cv1[3:]:
        a = i.strip().split(' ')
        index1.add(a[0])
        dict1[a[0]] = float(a[1])
    for i in sp_cv2[3:]:
        a = i.strip().split(' ')
        index2.add(a[0])
        dict2[a[0]] = float(a[1])
    common_sets = index1&index2
    dist = 0
    for i in common_sets:
        dist += (dict1[i]*dict2[i])
    return (1-dist/(math.sqrt(inner_product_for_sp1*inner_product_for_sp2)))/2

#Calculate the euclidean distance between two species' Kmer frequency vector
def euclidean_distance(pairsSpecies):
    sp1, sp2 = pairsSpecies[0], pairsSpecies[1]
    dict1, dict2, index1, index2 = {}, {}, set(), set()
    for i in open(directory+sp1+'.cv.txt', 'r').readlines()[3:]:
        a = i.strip().split(' ')
        index1.add(a[0])
        dict1[a[0]] = float(a[1])
    for i in open(directory+sp2+'.cv.txt', 'r').readlines()[3:]:
        a = i.strip().split(' ')
        index2.add(a[0])
        dict2[a[0]] = float(a[1])
    common_sets = index1&index2
    different_sets = (index1|index2)-common_sets
    dist = 0
    for i in common_sets:
        dist += (dict1[i]-dict2[i])**2
    for i in different_sets:
        if i in dict1:
            dist += dict1[i]**2
        else:
            dist += dict2[i]**2
    return math.sqrt(dist)


if __name__ == '__main__':
    help = "python3 distance_matrix.py cvtxt.directory CPU.NUM distance-type(COS or EUC)\n\
cvtxt.directory\t-str\tThe path of cv.txt of species stored\tother suffix files must nonexistent\n\
CPU.NUM\t-num\tThe number of CPU be used\n\
distance-type\t-str\tEuclidean distance or Cosine disimilarity\n"
    if len(sys.argv) < 4:
        exit(help)
    now = time.time()
    directory, coreNum, distance_type = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    specieslist = \
    [i[:-7] for i in os.listdir(directory) if i[-7:] == '.cv.txt']
    biospecies = list(itertools.combinations(specieslist,2))
    rst = []
    if distance_type == 'COS':
        p = multiprocessing.Pool(coreNum)
        rst = p.map(cosine_dissimilariry, biospecies)
        p.close()
        p.join()
    elif distance_type == 'EUC':
        p = multiprocessing.Pool(coreNum)
        rst = p.map(euclidean_distance, biospecies)
        p.close()
        p.join()
    else:
        exit(help)
    rst = list(rst)[::-1]

    #Output a triangle distance matrix of species
    output = \
    '#mega\n!Title megadist.meg' + ';\n' \
    + '!Format DataType=distance DataFormat=LowerLeft NTaxa=' \
    + str(len(specieslist))+';\n'+'\n'
    for i in specieslist[::-1]:
        output += '#'+i+'\n'
    output += '\n'
    count = 0
    for i in range(1,len(specieslist)+1):
        for j in rst[count:count+i]:
            output += str("%.18f"%j)+'\t'
            count += 1
        output += '\n'
    crtdir('./', 'MEG')
    crtdir('MEG/', directory[:-1])
    if distance_type == 'COS':
        open('MEG/'+directory+'/'+directory[:-1]+'cosdist.meg', 'w').write(output)
    else:
        open('MEG/'+directory+'/'+directory[:-1]+'eucdist.meg', 'w').write(output)
    print(time.time()-now)
