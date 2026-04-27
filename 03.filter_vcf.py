import sys
import gzip

file_i = sys.argv[1]
geno = 1-float(sys.argv[2])
n = 0
l = open(file_i+"."+sys.argv[2]+".info", 'w')
l.write("chr\tpos\tqual\tdpcount\tvarnum\tAAAA\tAnBn\tABC\tABCD\n")

def count_DP(lst):
    count = 0
    varnum = 0
    h1 = 0
    h2 = 0
    h3 = 0
    h4 = 0
    for i in lst:
        if i.startswith("."):
            continue
        gt, ad, dp, *others = i.split(":")
        if dp == ".":
            continue
        dp = int(dp)
        if dp > 0:
            count += 1
            gt_list = gt.split("/")
            if set(gt_list) != {"0"}:
                varnum += 1
            if len(set(gt_list)) == 1:
                h1 += 1
            elif len(set(gt_list)) == 2:
                h2 += 1
            elif len(set(gt_list)) == 3:
                h3 += 1
            elif len(set(gt_list)) == 4:
                h4 += 1
    return count, varnum, h1, h2, h3, h4

def clean_vcf(line):
    i_list = line.split("\t")
    data = i_list[9:]
    for j in range(len(data)):
        if data[j].startswith("."):
            data[j] = "./././."
        else:
            data[j] = data[j].split(":")[0]
    okline = "{}\t.\t.\tPR\tGT\t{}".format("\t".join(i_list[0:5]), "\t".join(data))
    return okline

       
with gzip.open (file_i, 'r') as f:
    for i in f:
        n += 1
        i = i.decode().strip()
        if i.startswith("#"):
            print(i)
            continue
        i_list = i.split("\t")
        CHROM, POS, ID, REF, ALT, QUAL = i_list[0:6]
        data = i_list[9:]
        count, varnum, h1, h2, h3, h4 = count_DP(data)
        l.write("{}\t{}\t{}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\n".format(CHROM, POS, QUAL, count, varnum, h1, h2, h3, h4))
        if count/len(data) > geno and varnum > 4 and float(QUAL) > 50:
            print(clean_vcf(i))

