import sys

kmerfa = sys.argv[1]
kmernum = sys.argv[2]

#保存kmer的坐标
kmer2pos = {}
with open (kmerfa, 'r') as f:
    for i in f:
        i = i.strip()
        if i.startswith(">"):
            seqid = i[1:]
        else:
            kmer2pos[i] = seqid

#统计
n = 0
with open (kmernum, 'r') as f:
    for i in f:
        n += 1
        if n == 1:
            continue
        i = i.strip()
        kmerid, *numlist = i.split()
        kmerid = kmerid.upper()
        numlist = list(map(int, numlist))
        Tnum, Dnum, Anum = sum(numlist[0:50]), sum(numlist[50:100]), sum(numlist)
        if kmerid in kmer2pos:
            print(kmerid, kmer2pos[kmerid].replace("_","\t"), Anum, Tnum, Dnum, sep="\t")
        else:
            print(kmerid)
