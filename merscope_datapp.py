from __future__ import division
import pandas as pd
import numpy as np
import time
import scipy.io
import pickle
import time
# import data_processing
# import PP
from data_processing.data_processing import DataLoader
from PP.proximal_pairs_test import ProximalPairs
import numpy as np
import pickle
import os
from InSTAnT.InSTAnT import Instant
import multiprocessing as mp
import timeit

def find_cell_id(args):
    n, x, y, fov = args[0], args[1], args[2], args[3]
    if fov in cell_meta:
        for possible_cell in cell_meta[fov]:
            if (possible_cell[1] <= x and x <= possible_cell[2]):
                if (possible_cell[3] <= y and y <= possible_cell[4]):
                    return n, possible_cell[0]

start = timeit.default_timer()
df_dt = pd.read_csv('data/merscope/ovc2_detected_transcripts.csv', index_col=0)
df_cm = pd.read_csv('data/merscope/ovc2_cell_metadata.csv', index_col=0)
print(len(df_dt.fov.unique()), len(df_cm.fov.unique()))

df_dt = df_dt.loc[df_dt.global_z == 0]
cell_meta = {}
for n, i in df_cm.iterrows():
    if i.fov not in cell_meta:
       cell_meta[int(i.fov)] = [[n, i.min_x, i.max_x, i.min_y, i.max_y]]
    else:
       cell_meta[int(i.fov)].append([n, i.min_x, i.max_x, i.min_y, i.max_y])

#df_dt = df_dt.iloc[:5]
#print(df_dt)
searches = []
for n, i in df_dt.iterrows():
   searches.append([n, i.global_x, i.global_y, i.fov])
   #df_fov = df_cm.loc[df_cm.fov == i.fov]
   #for cell_id, possible_cell in df_fov.iterrows():
   #    if (possible_cell.min_x < i.global_x and i.global_x < possible_cell.max_x):
   #        if (possible_cell.min_y < i.global_y and i.global_y < possible_cell.max_y):
   #            df_dt.at[n, 'uID'] = cell_id
   #            break
pool = mp.Pool(8)
results = pool.map(find_cell_id, searches)
with open("p_list.pkl", 'wb') as fp:
   pickle.dump(results, fp)

print("Map Time: ", timeit.default_timer() - start)
for result in results:
   if result is not None:
       df_dt.at[result[0], 'uID'] = result[1]
#print(results)
df_dt.to_csv('data/merscope/data.csv')

df_dt = pd.read_csv('data/merscope/data.csv', index_col=0)
print(len(df_dt),df_dt.uID.isna().sum(), len(df_dt) - df_dt.uID.isna().sum())

# 9975.304   5625.6196