import sys
import gzip

reads_file = sys.argv[1]
fq_file = sys.argv[2]
kmer_file = sys.argv[3]

readsok = {}
with open (reads_file, 'r') as f:
    for i in f:
        i = i.strip()
        readsok[i] = 1

kmerok = {}
with open (kmer_file, 'r') as f:
    for i in f:
        i = i.strip()
        kmerok[i] = 1

def read_fq_file_and_print(fq_file, readsok):
    n = 0
    with gzip.open (fq_file, "r") as f:
        for i in f:
            i = i.decode().strip()
            n += 1
            m = n%4
            if m == 1:
                seqid = i.split(" ")[0][1:]
            elif m == 2:
                seq = i
            elif m == 3:
                plus = i
            elif m == 0:
                qulity = i
                if seqid in readsok:
                    for j in kmerok:
                        if j in seq:
                            print("@{}\n{}\n{}\n{}".format(seqid, seq, plus, qulity))
                            break

read_fq_file_and_print(fq_file, readsok)
