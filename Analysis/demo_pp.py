# %pylab inline
from __future__ import division
import time
# import data_processing
# import PP
from data_processing.data_processing_u2_os_zhuang import DataLoader
from PP.proximal_pairs import ProximalPairs
import numpy as np
import pickle
import os


#####Data processing ###
dir_zhuang = '/Users/anurendrakumar/Desktop/1Nature_paper/Data/Zhuang_u2os/'
rep = 'rep2'
data_path = dir_zhuang + rep +  '/data.csv'
codebook_path = dir_zhuang + 'codebook.csv'
min_intensity = 10**0.75
min_area = 3
print('min intensity and area', min_intensity, min_area)
dataset = DataLoader(data_path,codebook_path, min_intensity = min_intensity, min_area = min_area)


dist=4
cell_id_list = dataset.df.uID.unique()
num_cells = len(cell_id_list)
print('num cells', num_cells)

save_dir = '/Users/anurendrakumar/Desktop/1Nature_paper/Results_final/Zhuang_u2os/' +  'PP_results/' 


if not os.path.exists(save_dir):
    os.makedirs(save_dir)


all_p_val = np.ones((num_cells, len(dataset.geneList), len(dataset.geneList)))
all_gene_count = np.zeros((num_cells, len(dataset.geneList)))
start_time = time.time()
min_genecount = 20

for i,cell_id in enumerate(cell_id_list[:num_cells]):
    df = dataset.df[dataset.df.uID == cell_id].copy()
    
    if df.shape[0]>min_genecount:
        pp_model = ProximalPairs(dataset.geneList, df,  dist_thresh = dist)
        all_p_val[i] = pp_model.p_vals
        all_gene_count[i] = pp_model.genecount.values.reshape(len(dataset.geneList))
    else:
        print('min genecount less than', min_genecount)
    
    if (i%50==0 or i == num_cells - 1):
        print("Cell  {0}, Time taken {1}.".format(i, time.time() - start_time))
        output = open(save_dir +   'pvalues_' + str(dist) + '.pkl', 'wb')
        pickle.dump(all_p_val, output)
        output.close()

        output2 = open(save_dir + 'gene_count' +'.pkl', 'wb')
        pickle.dump(all_gene_count, output2)
        output2.close()
        start_time = time.time()
        # pass
