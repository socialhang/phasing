import sys
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import kendalltau
import random

queryfile = sys.argv[1]
refile = sys.argv[2]
cutoff = 0.5

kmer2lst = {}
with open (queryfile, 'r') as f:
    for i in f:
        i = i.strip()
        kmerid, *value = i.split()
        if kmerid == "kmerid":
            continue
        value = list(map(int, value))
        kmer2lst[kmerid] = value
n = 0
with open (refile, 'r') as f:
    for i in f:
        i = i.strip()
        n += 1
        if n == 1:
            continue
        kmerid, *value = i.split()
        value = list(map(int, value))
        for i in kmer2lst:
            cor, pv = kendalltau(np.array(kmer2lst[i]), np.array(value))
            if cor < cutoff:
                continue
            #过滤和保留
            #cor = (cor + 1)/2  #全部数据，包括负相关的
            cor = (cor-cutoff)/(1-cutoff)  #去掉小于cutoff的值
            rand_num = random.random()
            m = 0
            if cor > rand_num:
                m += 1
                chrid1, pos1, strand1 = i.split("_")
                chrid2, pos2, strand2 = kmerid.split("_")
                if chrid1 == chrid2 and pos1 == pos2:
                    continue
                frag1 = random.randint(0, 1)
                frag2 = 1 - frag1
                #print("{:d}\t{}\t{}\t{}\t{}\t{}\t{}\t60".format(m, chrid1, pos1, chrid2, pos2, strand1, strand2)) # cphasing format
                print(16, chrid1, pos1, frag1, 16, chrid2, pos2, frag2, "60 - - 60 - - -")# 3ddna format

#cphasing format
'''
0	utg001511l	52687	utg002211l	462632	-	+	0
1	utg001511l	52687	utg002480l	4989	-	-	0
2	utg001511l	52687	utg002480l	5843	-	+	0
3	utg001511l	52687	utg004459l	49529	-	+	0
4	utg002211l	462632	utg002480l	4989	+	-	1
5	utg002211l	462632	utg002480l	5843	+	+	1
6	utg002211l	462632	utg004459l	49529	+	+	1
7	utg002480l	4989	utg002480l	5843	-	+	1
8	utg002480l	4989	utg004459l	49529	-	+	1
9	utg002480l	5843	utg004459l	49529	+	+	1
'''





