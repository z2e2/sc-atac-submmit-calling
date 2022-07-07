import pandas as pd

submmit_dir = "results/pbmc_archR_output/"
peak_cells = pd.read_csv(submmit_dir + "PBMC_cellnames_info.txt", sep = "\t")
cell_set = set([i.split('#')[1] for i in list(peak_cells.index)])

##### start-end adjusted #####
# import gzip
# fragment_file = '/Genomics/pritykinlab/zzhao/sc-atac-submmit-calling/data/pbmc10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
# STEP = 2363931
# frag_file = gzip.open(fragment_file, "rb")
# res = {}
# keys = []
# for idx, line in enumerate(frag_file):
#     if idx % STEP == 0:
#         print(f'{idx/STEP}%')
#     line = line.decode().strip('\n')
#     if line[0] == '#':
#         continue
#     chrm, start, end, cell_id, count = line.split('\t')
#     if cell_id not in cell_set:
#         continue
#     start_key = f'{chrm}:{int(start)+5}'
#     end_key = f'{chrm}:{int(end)-4}'
#     if start_key in res:
#         res[start_key] += 1
#     else:
#         res[start_key] = 1
#         keys.append(start_key)
#     if end_key in res:
#         res[end_key] += 1
#     else:
#         res[end_key] = 1
#         keys.append(end_key)  
# frag_file.close()   

import gzip

fragment_file = '/Genomics/pritykinlab/zzhao/sc-atac-submmit-calling/data/pbmc10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz'
STEP = 2363931
frag_file = gzip.open(fragment_file, "rb")

res = {}
keys = []
for idx, line in enumerate(frag_file):
    if idx % STEP == 0:
        print(f'{idx/STEP}%')
    line = line.decode().strip('\n')
    if line[0] == '#':
        continue
    chrm, start, end, cell_id, count = line.split('\t')
    if cell_id not in cell_set:
        continue
    if chrm[:3] != 'chr':
        continue
    start_key = f'{chrm}:{start}'
    end_key = f'{chrm}:{end}'
    if start_key in res:
        res[start_key] += 1
    else:
        res[start_key] = 1
        keys.append(start_key)
    if end_key in res:
        res[end_key] += 1
    else:
        res[end_key] = 1
        keys.append(end_key)
        
        
frag_file.close()   

write_file = open('/Genomics/pritykinlab/zzhao/sc-atac-submmit-calling/fragments_tn5.bed', 'w+')
STEP = len(keys)//10
for idx, key in enumerate(keys):
    if idx % STEP == 0:
        print(f'{10*idx/STEP}%')
    chrm, start = key.split(':')
    count = res[key]
    line = f'{chrm}\t{int(start)}\t{int(start)+1}\t{count}\n'
    write_file.write(line)
write_file.close()
        
