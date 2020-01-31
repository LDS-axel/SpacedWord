Spaced word
==================================
DATE:2020-1-31
version:0.2.0


Author Information:
==================
Name:  Liu Shuai
Department:  Anhui University and the Institute of Zoology, Chinese Academy of Sciences jointly cultivation
E-mail:  ls2659614061@126.com


Supervisor Information:
======================
Name:  Wu Qi
Department:  Institute of Microbiology, Chinese Academy of Sciences
E-mail:  wuq@im.ac.cn


General Procedure
=================
* kmer frequency vector
    * 1.Pretreatment and distribution for scaffolds-like or chromosomes-like genomes sequences.
    * 2.Scanning each sequence with a k-string(contiguous length k word or 
    degenerated word with a length k but less than k characters is count),
    count each k-string's number in each sequence.
    * 3.Make a summary statistic for all k-strings' number in all scaffolds or chromosomes.
    * 4.Calculate the frequence of k-strings subtracted by 
    contiguous k-strings background frequence or 
    degenerated k-strings background frequence or
    Calculate the frequence of k-strings straightly without any subtraction.
    * 5.Make a cv.txt format output.
        #example.cv.txt
        N
        Inner Product
        Significant kmer number
        1 0.044418402253529
        2 -0.13694127815637536
        3 0.02302543931397616
        4 0.08209917560394331
        ......
* distance matrix
    * 1.Employee Cosine dissimilarity or Euclidean distance to calculate the difference among species
    * 2.Make a .meg format output opened up by MEGA.


Manual
========
kmer frequency vector
* python3 kmer_frequence_vector.py -d -m -fna.path -speciesname -CPU.NUM -bgd
-d      -num    the interval length of degenerated kmer     d>0
-N      -num    the length of degenerated kmer      N>1
-fna.path       -str    the path of dna sequence as an input,whose suffix is fna,fa,ffn,or fasta
-speciesname    -str    the Latin name of species
-CPU.NUM        -num    the number of CPU be used
-bgd            -str    whether subtracted background

distance matrix
* python3 distance_matrix.py cvtxt.directory CPU.NUM distance-type(COS or EUC)
cvtxt.directory     -str    The path of cv.txt of species stored      other suffix files must nonexistent
CPU.NUM     -num            The number of CPU be used
distance-type       -str    Euclidean distance or Cosine disimilarity