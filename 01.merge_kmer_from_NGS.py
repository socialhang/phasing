import sys
import gzip
import resource

def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <file_list> <btype>")
        sys.exit(1)
    
    NGS_kmer_count_list = sys.argv[1]
    btype = sys.argv[2]
    
    # 读取文件列表并预处理
    files = []
    groups = []
    with open(NGS_kmer_count_list, 'r') as f:
        for idx, line in enumerate(f):
            filename = line.strip()
            group = filename.replace("_uG.list.gz", "")
            files.append(filename)
            groups.append(group)
    
    total_files = len(files)
    # 打印表头
    print("kmerid", "\t".join(groups), sep="\t")
    
    # 初始化kmer计数字典
    kmer_bitmasks = {}
    log_file = f"{NGS_kmer_count_list}.{btype}.log"
    
    # 处理每个文件
    with open(log_file, 'w') as log_handle:
        for file_index, filename in enumerate(files):
            groupid = filename.replace("_uG.list.gz", "")
            file_bit = 1 << file_index  # 位掩码
            
            try:
                with gzip.open(filename, 'rt') as f:
                    next(f)  # 跳过标题行
                    for line in f:
                        line = line.strip()
                        kmerseq = line.split(":")[1].split("=")[0].strip()
                        a = 0
                        if (a == 0) and (not kmerseq.startswith(btype)):
                            continue
                        if kmerseq.startswith(btype):
                            a = 1
                            if kmerseq not in kmer_bitmasks:
                                kmer_bitmasks[kmerseq] = 0
                            kmer_bitmasks[kmerseq] |= file_bit
                        else:
                            break
                        
                    
            except Exception as e:
                print(f"Error processing {filename}: {str(e)}", file=sys.stderr)
            
            # 记录内存使用
            mem_usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            log_handle.write(f"{groupid}\tok\t{mem_usage/1024/1024:.2f}Gb\n")
            log_handle.flush()
    
    # 输出结果
    for kmerseq, bitmask in kmer_bitmasks.items():
        bits = []
        for i in range(total_files):
            bits.append('1' if bitmask & (1 << i) else '0')
        print(kmerseq, "\t".join(bits), sep="\t")

if __name__ == '__main__':
    main()
