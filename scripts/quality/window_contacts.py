import sys

file_i = sys.argv[1]
windows = 10000

with open (file_i, 'r') as f:
    for i in f:
        i = i.strip()
        _, c1, p1, _, _, c2, p2, *others = i.split()
        p1, p2 = int(int(p1)/windows), int(int(p2)/windows)
        if c1+"_"+str(p1) != c2+"_"+str(p2):
            print(c1+"_"+str(p1), c2+"_"+str(p2), sep="\t")
