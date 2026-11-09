import sys

file_i = sys.argv[1]
winlen = 100000

def max_index(lst):
    if not lst:
        return None
    max_val = lst[0]
    max_idx = 0
    for i, val in enumerate(lst):
        if val > max_val:
            max_val = val
            max_idx = i
    return max_idx

pp2cor = {}
with open (sys.argv[1], 'r') as f:
    for i in f:
        i = i.strip()
        rawpos1, rawpos2, cor = i.split()
        if not rawpos1.startswith("chr") or not rawpos2.startswith("chr"):
            continue
        c1, p1 = rawpos1.split("_")
        c2, p2 = rawpos2.split("_")
        g1, h1 = c1.split(".")
        g2, h2 = c2.split(".")
        if g1 != g2:
            continue
        #p1, p2 = int(int(p1)/winlen)*winlen, int(int(p2)/winlen)*winlen
        pos1, pos2 = "{}_{}".format(c1, str(p1)), "{}_{}".format(c2, str(p2))
        cor = float(cor)
        if pos1 not in pp2cor:
            pp2cor[pos1] = [0,0,0,0]
        pp2cor[pos1][int(h2)-1] += cor
        if pos2 not in pp2cor:
            pp2cor[pos2] = [0,0,0,0]
        pp2cor[pos2][int(h1)-1] += cor

for pos in pp2cor:
    hapid = pos.split("_")[0][-1]
    print(pos, hapid, max_index(pp2cor[pos])+1, pp2cor[pos][0], pp2cor[pos][1], pp2cor[pos][2], pp2cor[pos][3], sep="\t")
