from __future__ import division
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import binom_test
from collections import Counter


class ProximalPairs():
    def __init__(self, geneList, df,  dist_thresh, mode='normal'):
        self.geneList = geneList
        self.orig_df = df.copy()
        self.curr_cell_df = df[['absX', 'absY']]

        self.dist_thresh = dist_thresh
        if mode=='normal':
            self.genecount, self.num_trial = self.num_trial_pairs()
            self.prob_null = self.prob_null_hypothesis()  
            self.obs = self.obs_trial_pairs()
            self.p_vals = self.compute_p_val()
        else:
            self.obs_all = self.obs_spatial_cat()
        

    def compute_p_val(self):
        p_vals = np.ones((len(self.geneList), len(self.geneList)))
        for i in range(self.obs.shape[0]):
            for j in range(i,self.obs.shape[1]):
                p_vals[i,j] = binom_test(self.obs[i,j], self.num_trial[i,j], \
                self.prob_null, alternative = 'greater')
                p_vals[j,i] = p_vals[i,j]
        
        return p_vals

    def prob_null_hypothesis(self):
        point_tree = cKDTree(self.curr_cell_df)
        pairs = point_tree.query_pairs(self.dist_thresh)
        if len(pairs):   #Later change the condition to min gene count and other heuristic
            prob_null = len(pairs)*2/(self.curr_cell_df.shape[0]*(self.curr_cell_df.shape[0]-1))
        else:
            prob_null = 1

        return prob_null




    def num_trial_pairs(self):  #Reindexed with geneList
        genecount = pd.DataFrame.from_dict(Counter(self.curr_cell_df.index) , orient='index').reindex(self.geneList).fillna(0)
        num_trial_pairs = genecount.dot(genecount.T).values   #n1*n2
        return genecount, num_trial_pairs

    def obs_trial_pairs(self): #debug for large d
        point_tree = cKDTree(self.curr_cell_df)
        pairs = point_tree.query_pairs(self.dist_thresh)

        # prob_null = len(pairs)*2/(self.curr_cell_df.shape[0]*(self.curr_cell_df.shape[0]-1))

        pairs = [(self.curr_cell_df.index[i], self.curr_cell_df.index[j]) for (i,j) in pairs]
        pairs = Counter(pairs)
        pairs = pd.Series(pairs).reset_index()

        if len(pairs):
            obs_df = pd.pivot_table(pairs, index = 'level_0', columns ='level_1', values = 0).fillna(0)
            col_na_genes = [i for i in self.geneList if i not in obs_df.columns]
            row_na_genes = [i for i in self.geneList if i not in obs_df.index]
            col_na_genes = dict.fromkeys(col_na_genes, 0)  #nan genes
            obs_df = obs_df.assign(**col_na_genes)
            for row_na in row_na_genes:
                obs_df.loc[row_na] = 0
            obs_df = obs_df.reindex(index = self.geneList, columns = self.geneList)
        else:                                       #if no entry less than dist thresh
            print('no entry less than dist thresh, total rna', self.curr_cell_df.shape[0])
            obs_df = pd.DataFrame(0, self.geneList, self.geneList)

        arr2 = np.triu(obs_df) + np.triu(obs_df,1).T   #Making matrix symmetric

        return arr2

    def obs_spatial_cat(self): #debug for large d
        spatial_cat_all = spatialCat(self.orig_df)
        # spatial_cat_all = spatialCatPercentile(self.orig_df)
        self.orig_df.loc[:,'InnerNuc'] = spatial_cat_all[0]
        self.orig_df.loc[:,'Perinuc'] = spatial_cat_all[1]
        self.orig_df.loc[:,'Cyto'] = spatial_cat_all[2]
        self.orig_df.loc[:,'CellPeri'] = spatial_cat_all[3]

        point_tree = cKDTree(self.curr_cell_df)
        queried_pairs = point_tree.query_pairs(self.dist_thresh)
        pairs = [(self.curr_cell_df.index[i], self.curr_cell_df.index[j]) for (i,j) in queried_pairs]

        pairs_isinnernuc = [(self.orig_df.iloc[i,:]['InnerNuc'] + self.orig_df.iloc[j,:]['InnerNuc'])/2 for (i,j) in  queried_pairs]  #update later
        pairs_isperinuc = [(self.orig_df.iloc[i,:]['Perinuc'] + self.orig_df.iloc[j,:]['Perinuc'])/2 for (i,j) in  queried_pairs] 
        pairs_iscyto = [(self.orig_df.iloc[i,:]['Cyto'] + self.orig_df.iloc[j,:]['Cyto'])/2 for (i,j) in  queried_pairs] 
        pairs_iscellperi = [(self.orig_df.iloc[i,:]['CellPeri'] + self.orig_df.iloc[j,:]['CellPeri'])/2 for (i,j) in  queried_pairs] 
        pairs,spatial_cat_0,spatial_cat_1,spatial_cat_2,spatial_cat_3  = AnuCounter(pairs,[pairs_isinnernuc, pairs_isperinuc, pairs_iscyto, pairs_iscellperi]).counter
        pairs = pd.Series(pairs).reset_index()

        cat_pairs_0 = pd.Series(spatial_cat_0).reset_index()
        cat_pairs_1 = pd.Series(spatial_cat_1).reset_index()
        cat_pairs_2 = pd.Series(spatial_cat_2).reset_index()
        cat_pairs_3 = pd.Series(spatial_cat_3).reset_index()

        if len(pairs):
            obs_df = pd.pivot_table(pairs, index = 'level_0', columns ='level_1', values = 0).fillna(0)
    
            obs_cat_0 = pd.pivot_table(cat_pairs_0 , index = 'level_0', columns ='level_1', values = 0).fillna(0)
            obs_cat_1 = pd.pivot_table(cat_pairs_1, index = 'level_0', columns ='level_1', values = 0).fillna(0)
            obs_cat_2 = pd.pivot_table(cat_pairs_2, index = 'level_0', columns ='level_1', values = 0).fillna(0)
            obs_cat_3 = pd.pivot_table(cat_pairs_3, index = 'level_0', columns ='level_1', values = 0).fillna(0)
            col_na_genes = [i for i in self.geneList if i not in obs_df.columns]
            row_na_genes = [i for i in self.geneList if i not in obs_df.index]
            col_na_genes = dict.fromkeys(col_na_genes, 0)  #nan genes
            obs_df = obs_df.assign(**col_na_genes)
            obs_cat_0 = obs_cat_0.assign(**col_na_genes)
            obs_cat_1 = obs_cat_1.assign(**col_na_genes)
            obs_cat_2 = obs_cat_2.assign(**col_na_genes)
            obs_cat_3 = obs_cat_3.assign(**col_na_genes)
            

            for row_na in row_na_genes:
                obs_df.loc[row_na] = 0
                obs_cat_0.loc[row_na] = 0
                obs_cat_1.loc[row_na] = 0
                obs_cat_2.loc[row_na] = 0
                obs_cat_3.loc[row_na] = 0
                

            obs_df = obs_df.reindex(index = self.geneList, columns = self.geneList)
            obs_cat_0 = obs_cat_0.reindex(index = self.geneList, columns = self.geneList)
            obs_cat_1 = obs_cat_1.reindex(index = self.geneList, columns = self.geneList)
            obs_cat_2 = obs_cat_2.reindex(index = self.geneList, columns = self.geneList)
            obs_cat_3 = obs_cat_3.reindex(index = self.geneList, columns = self.geneList)
        else:                                       #if no entry less than dist thresh
            print('no entry less than dist thresh, total rna', self.curr_cell_df.shape[0])
            obs_df = pd.DataFrame(0, self.geneList, self.geneList)
            obs_cat_0,obs_cat_1,obs_cat_2,obs_cat_3  = obs_df.copy(), obs_df.copy(), obs_df.copy(),obs_df.copy()

        arr_cat_0 = np.triu(obs_cat_0) + np.triu(obs_cat_0,1).T 
        arr_cat_1 = np.triu(obs_cat_1) + np.triu(obs_cat_1,1).T 
        arr_cat_2 = np.triu(obs_cat_2) + np.triu(obs_cat_2,1).T 
        arr_cat_3 = np.triu(obs_cat_3) + np.triu(obs_cat_3,1).T 

        return arr_cat_0, arr_cat_1, arr_cat_2, arr_cat_3



    
def spatialCat(df,dist_thresh_nucleus = 2.5,dist_thresh_cyto_nuc = 2.5,dist_thresh_cyto_peri = 4):
    
    nuc_distance = df.distNucleus.values
    cyt_distance = df.distPeriphery.values
    inNucleus = df.inNucleus.values
    
    
    
    temp1 = np.zeros(nuc_distance.shape)
    temp2 = np.zeros(nuc_distance.shape)
    temp3 = np.zeros(nuc_distance.shape)
    temp4 = np.zeros(nuc_distance.shape)

    temp1[(inNucleus==1) & (nuc_distance>dist_thresh_nucleus )] = 1
    temp2[(inNucleus==1) & (nuc_distance<=dist_thresh_nucleus )] = 1
    temp2[(inNucleus==0) & (nuc_distance<=dist_thresh_cyto_nuc )] = 1
    temp3[(inNucleus==0) & (nuc_distance>dist_thresh_cyto_nuc) & (cyt_distance>dist_thresh_cyto_peri)] = 1
    temp4[(inNucleus==0) & (cyt_distance<=dist_thresh_cyto_peri)] = 1


    return [temp1,temp2,temp3, temp4]

def spatialCatPercentile(df, nuc_percentile=25, cyt_percentile_nuc = 75, cyt_percentile_peri= 30): #40,75,30-previous
    
    nuc_distance = df.distNucleus.values
    cyt_distance = df.distPeriphery.values

    temp_nucleus = df[df.inNucleus==1]
    temp_cyt = df[df.inNucleus ==0]
    dist_thresh_nucleus = np.percentile(temp_nucleus.distNucleus,nuc_percentile)
    dist_thresh_cyt_nuc = np.percentile(temp_cyt.distPeriphery,cyt_percentile_nuc)
    dist_thresh_cyt_peri = np.percentile(temp_cyt.distPeriphery,cyt_percentile_peri)
    
    inNucleus = df.inNucleus.values
    temp1 = np.zeros(nuc_distance.shape)
    temp2 = np.zeros(nuc_distance.shape)
    temp3 = np.zeros(nuc_distance.shape)
    temp4 = np.zeros(nuc_distance.shape)

    temp1[(inNucleus==1) & (nuc_distance>dist_thresh_nucleus)] = 1
    temp2[(inNucleus==1) & (nuc_distance<=dist_thresh_nucleus)] = 1
    temp2[(inNucleus==0) & (cyt_distance>=dist_thresh_cyt_nuc)] = 1
    temp3[(inNucleus==0) & (cyt_distance<dist_thresh_cyt_nuc ) & (cyt_distance>dist_thresh_cyt_peri )] = 1
    temp4[(inNucleus==0) & ( (cyt_distance<dist_thresh_cyt_peri ))] = 1
    
    return [temp1,temp2,temp3, temp4]



def anu_count_elements(pairs, pairs_all_cat):
    counter = {}
    spatial_cat_0,spatial_cat_1,spatial_cat_2,spatial_cat_3 = {},{},{},{}
    pairs_isinnernuc, pairs_isperinuc, pairs_iscyto, pairs_iscellperi = pairs_all_cat[0],pairs_all_cat[1],pairs_all_cat[2],pairs_all_cat[3]
    
    for i,elem in enumerate(pairs):
        
        counter[elem] = counter.get(elem, 0) + 1
        spatial_cat_0[elem] = spatial_cat_0.get(elem, 0) + pairs_isinnernuc[i]
        spatial_cat_1[elem] = spatial_cat_1.get(elem, 0) + pairs_isperinuc[i]
        spatial_cat_2[elem] = spatial_cat_2.get(elem, 0) + pairs_iscyto[i]
        spatial_cat_3[elem] = spatial_cat_3.get(elem, 0) + pairs_iscellperi[i]

    return counter, spatial_cat_0,spatial_cat_1,spatial_cat_2,spatial_cat_3#,spatial_cat_4
        
class AnuCounter:
    def __init__(self, pairs,pairs_all_cat):
        
        self.counter = anu_count_elements(pairs, pairs_all_cat)

