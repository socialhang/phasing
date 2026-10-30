import sys
import random

bed_file = sys.argv[1]
out_file = sys.argv[2]
window = int(sys.argv[3])

samplenum = 1
lowcutoff = 6
higcutoff = 180

bed2kmer = {}
with open (bed_file, 'r') as f:
    for i in f:
        i = i.strip()
        chrid, st, end = i.split()
        st, end = int(st), int(end)
        bed2kmer[(chrid, st)] = []

with open (out_file, 'r') as f:
    for i in f:
        i = i.strip()
        if len(i.split()) < 2:
            continue
        if len(i.split()) == 8:
            kmerseq, chrid1, chrid2, pos, strand, Anum, Tnum, Dnum = i.split()
            chrid = chrid1 + "_" + chrid2
        else:
            kmerseq, chrid, pos, strand, Anum, Tnum, Dnum = i.split()
        pos, Anum, Tnum, Dnum = int(pos), int(Anum), int(Tnum), int(Dnum)
        if Anum < lowcutoff or Anum > higcutoff:
            continue
        pos = int(pos/window)*window
        if (chrid, pos) in bed2kmer:
            bed2kmer[(chrid,pos)].append(kmerseq)

for j in bed2kmer:
    if len(bed2kmer[j]) >= samplenum:
        for k in random.sample(bed2kmer[j], samplenum):
            print(k)
