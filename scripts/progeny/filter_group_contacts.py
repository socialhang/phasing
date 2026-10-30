import sys

key_file = sys.argv[1]
merged_file = sys.argv[2]

ok = {}
with open (key_file, 'r') as f:
    for i in f:
        i = i.strip()
        ok[i] = 1

'''
16 utg000830l 1320066 1 16 utg000053l 655284 0 60 - - 60 - - -
16 utg000402l 2799253 0 16 utg000286l 3585308 1 60 - - 60 - - -
16 utg000719l 2737980 1 16 utg000719l 3461827 0 60 - - 60 - - -
16 utg000446l 6522562 1 16 utg000446l 1743132 0 60 - - 60 - - -
16 utg000053l 655284 1 16 utg000830l 1320066 0 60 - - 60 - - -
'''
with open (merged_file, 'r') as f:
    for i in f:
        i = i.strip()
        i_list = i.split()
        if i_list[1] in ok and i_list[5] in ok:
            print(i)
