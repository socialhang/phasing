import sys

hap2seq = {}
with open (sys.argv[1], 'r') as f:
    for i in f:
        i = i.strip()
        if i.startswith(">"):
            hapid = i.split(":")[0][1:]
        else:
            if hapid not in hap2seq:
                hap2seq[hapid] = i
            else:
                hap2seq[hapid] += i

for seqid in hap2seq:
    print(">{}\n{}".format(seqid, hap2seq[seqid]))
