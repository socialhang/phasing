import sys
import gzip
#根据readsid提取gz格式fq文件list的reads

fq_file_list = sys.argv[1] #gzip格式fq文件list
readsid_list = sys.argv[2]

def read_ok_reads_id(readsid_list):
    readsok = dict()
    with open (readsid_list, 'r') as f:
        for i in f:
            i = i.strip()
            readsid, *other = i.split("\t")
            readsok[readsid+"/1"] = 1
            readsok[readsid+"/2"] = 1
            readsok[readsid] = 1
    return readsok

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
                    print("@{}\n{}\n{}\n{}".format(seqid, seq, plus, qulity))

if __name__ == "__main__":

    readsok = read_ok_reads_id(readsid_list)
    with open (fq_file_list, 'r') as f:
        for i in f:
            i = i.strip()
            read_fq_file_and_print(i, readsok)
