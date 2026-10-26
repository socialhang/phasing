#!/usr/bin/env python3
#输出文件第一列是utgid，KDE，average
import sys
import gzip
import numpy as np
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

ctg2depth_file = sys.argv[1]
windows = 100

# 4 min
win2depth = {}
with gzip.open (ctg2depth_file, 'r') as f:
    for i in f:
        i = i.decode().strip()
        utgid, st, end, depth = i.split()
        st, end, depth = int(st), int(end), int(depth)
        if depth == 0:
            continue
        for pos in range(st, end):
            win = int(pos/windows)
            if utgid not in win2depth:
                win2depth[utgid] = {win:[depth,1]}
            elif win not in win2depth[utgid]:
                win2depth[utgid][win] = [depth,1]
            else:
                win2depth[utgid][win][0] += depth
                win2depth[utgid][win][1] += 1

def average(lst):
    return lst[0]/lst[1]

# 1 h
utg2depth = {}
for utgid in win2depth:
    rawlist = list(win2depth[utgid].values())
    utg2depth[utgid] = list(map(average, rawlist))

for k in utg2depth:
    if len(utg2depth[k]) < 5:
        print(k, sum(utg2depth[k])/len(utg2depth[k]), sum(utg2depth[k])/len(utg2depth[k]), sep="\t")
        continue
    data = np.array(utg2depth[k])
    try:
        kde = gaussian_kde(data, bw_method='scott')
        x = np.linspace(data.min()-1, data.max()+1, 1000)
        density = kde.evaluate(x)
        max_index = np.argmax(density)
        peak_value = x[max_index]
        print(k, peak_value, sum(utg2depth[k])/len(utg2depth[k]), sep="\t")
    except np.linalg.LinAlgError as e:
        print(k, sum(utg2depth[k])/len(utg2depth[k]), sum(utg2depth[k])/len(utg2depth[k]), sep="\t")
