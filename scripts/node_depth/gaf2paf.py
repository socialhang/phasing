#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import pybedtools
import re

rawgfa = sys.argv[1]
gaf_file = sys.argv[2]
pattern = r'([><]utg[0-9]+[lc])'

def give_ax(graph, graphlen, graphst, graphend, utg2length, link2olap):
    rawutglist = re.findall(pattern, graph)
    regionlst = []
    if len(rawutglist) == 1:
        utgid = rawutglist[0][1:]
        utglen = utg2length[utgid]
        if rawutglist[0].startswith("<"):
            graphst = utglen - graphst + 1
            graphend = utglen - graphend + 1
            graphst, graphend = graphend, graphst
            strand = "-"
        else:
            strand = "+"
        regionlst.append((utgid, graphst, graphend, strand))
    else:
        st = 1
        bed = []
        cut = [] #用于保存矫正坐标
        cut.append(st-1)
        for j in range(len(rawutglist)-1):
            utg1, utg2 = rawutglist[j], rawutglist[j+1]
            utgid1, utgid2 = utg1[1:], utg2[1:]
            utg1len, utg2len = utg2length[utgid1], utg2length[utgid2]
            overlap = link2olap[(utg1, utg2)]
            if j == 0:
                bedinfo = "graph\t{:d}\t{:d}\t{}".format(st, utg1len, utg1)
                bed.append(bedinfo)
            st = st + utg1len - 1 - overlap + 1
            cut.append(st)
            bedinfo = "graph\t{:d}\t{:d}\t{}".format(st, st+utg2len, utg2)
            bed.append(bedinfo)
        bed1 = pybedtools.BedTool(bed)
        bed2 = pybedtools.BedTool(["graph\t{:d}\t{:d}".format(graphst, graphend)])
        overlap_result = bed1.intersect(bed2)
        k = -1 #矫正一下坐标，将graph坐标转成utg坐标
        for interval in overlap_result:
            k += 1
            utg, st, end = interval.name, interval.start-cut[k], interval.end-cut[k]
            if utg.startswith("<"):
                utglen = utg2length[utg[1:]]
                st = utglen - st + 1
                end = utglen - end + 1
                st, end = end, st
                strand = "-"
            else:
                strand = "+"
            regionlst.append((utg[1:], st, end, strand))
    return regionlst

# read gfa
utg2length = {}
link2olap = {}
with open (rawgfa, 'r') as f:
    for i in f:
        i = i.strip()
        if i.startswith("L"):
            ltype, ctg1, strand1, ctg2, strand2, overlaps, *others = i.split()
            if strand1 == "+":
                ctg1 = ">" + ctg1
            else:
                ctg1 = "<" + ctg1
            if strand2 == "+":
                ctg2 = ">" + ctg2
            else:
                ctg2 = "<" + ctg2
            overlaps = int(overlaps[0:-1])
            link2olap[(ctg1, ctg2)] = overlaps
        elif i.startswith("S"):
            ltype, utg, _, utglen, _ = i.split()
            utglen = int(utglen[5:])
            utg2length[utg] = utglen

readid_f = ""
with open (gaf_file, 'r') as f:
    for i in f:
        i = i.strip()
        readid, readlen, readst, readend, _, graph, graphlen, graphst, graphend, *others = i.split()
        #print(readlen, readst, readend)
        readlen, readst, readend = int(readlen), int(readst), int(readend)
        graphlen, graphst, graphend = int(graphlen), int(graphst), int(graphend)
        if readlen < 1000:
            continue
        #print(readid, readlen, readst, readend, readlen, graph, graphlen, graphst, graphend)
        #节奏控制
        if readid != readid_f:
            tmp = []
            escape = 0
            readid_f = readid
        elif escape == 1:
            continue
        region_lst = ""
        if tmp == [] and (readend - readst)/readlen > 0.9:
            region_lst = give_ax(graph, graphlen, graphst, graphend, utg2length, link2olap)
            escape = 1
        elif tmp == [] and (readend - readst)/readlen <= 0.9:
            tmp = list(range(int(readst/100), int(readend/100)))
            region_lst = give_ax(graph, graphlen, graphst, graphend, utg2length, link2olap)
        else:
            newlst = list(range(int(readst/100), int(readend/100)))
            overlap = set(newlst) & set(tmp)
            if len(overlap)/readlen < 0.05 and readend-readst+1>100:
                region_lst = give_ax(graph, graphlen, graphst, graphend, utg2length, link2olap)
            tmp = list(set(newlst) | set(tmp))
        #对于region_lst处理和输出
        if region_lst != "":
            for j in region_lst:
                print(readid, readlen, j[0], j[1], j[2], j[3], sep="\t")
