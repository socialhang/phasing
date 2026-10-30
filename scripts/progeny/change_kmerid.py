import sys

fafile = sys.argv[1]
kmerfile = sys.argv[2]

kmer2pos = {}
with open (fafile, 'r') as f:
    for i in f:
        i = i.strip()
        if i.startswith(">"):
            seqid = i[1:]
        else:
            kmer2pos[i] = seqid

n = 0
with open (kmerfile, 'r') as f:
    for i in f:
        i = i.strip()
        n += 1
        if n == 1:
            print(i)
            continue
        kmerid, other = i.split("\t", 1)
        kmerid = kmerid.upper()
        print(kmer2pos[kmerid], other, sep="\t")
