import sys

key_file = sys.argv[1]
all_file = sys.argv[2]

key = {}
with open (key_file, 'r') as f:
    for i in f:
        i = i.strip()
        key[i] = 1

n = 0
with open (all_file, 'r') as f:
    for i in f:
        n += 1
        i = i.strip()
        if n == 1:
            print(i)
            continue
        kmer, *others = i.split()
        if kmer.upper() in key:
            print(i)
