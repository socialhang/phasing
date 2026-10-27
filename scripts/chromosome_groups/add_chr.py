import sys

file_i = sys.argv[1]
ids = file_i.split(".")[0]

utgid2chrid = {}
with open (file_i, 'r') as f:
    for i in f:
        i = i.strip()
        utgid2chrid[i] = 1
        print(i)

with open ("utg2chr.txt", 'r') as f:
    for i in f:
        i = i.strip()
        utgid, chrid = i.split()
        if chrid != ids:
            continue
        if utgid not in utgid2chrid:
            print(utgid)
