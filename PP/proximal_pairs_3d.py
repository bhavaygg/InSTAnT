from __future__ import division
import pandas as pd
import time
import numpy as np
from scipy.spatial import cKDTree, distance_matrix
from scipy.stats import binom_test
from collections import Counter

class ColocalizationModel():
    def __init__(self, geneList, df,  dist_thresh):  #Later remove dist_thresh from here?
        start_time = time.time()
        self.geneList = geneList
        self.curr_cell_df = df
        self.curr_cell_df = self.curr_cell_df[['absX', 'absY','geneName']].set_index('geneName') #Later change this implementation
        self.dist_thresh = dist_thresh
        self.genecount, self.num_trial_df = self.num_trial_pairs()
        self.prob_null, _ = self.prob_null_hypothesis()
        self.obs_df = self.obs_trial_pairs()
        # self.p_vals_df = self.compute_p_val()
        

    def compute_p_val(self):
        p_vals = np.ones((len(self.geneList), len(self.geneList)))
        for i in range(self.obs_df.shape[0]):
            for j in range(i,self.obs_df.shape[1]):
                p_vals[i,j] = binom_test(self.obs_df.iloc[i,j], self.num_trial_df.iloc[i,j], \
                self.prob_null, alternative = 'greater')
                p_vals[j,i] = p_vals[i,j]
        p_vals_df = pd.DataFrame(p_vals, index = self.geneList, columns = self.geneList)
        return p_vals_df

    def prob_null_hypothesis(self):
        point_tree = cKDTree(self.curr_cell_df)
        pairs = point_tree.query_pairs(self.dist_thresh)
        if len(pairs):   #Later change the condition to min gene count and other heuristic
            prob_null = len(pairs)*2/(self.curr_cell_df.shape[0]*(self.curr_cell_df.shape[0]-1))
        else:
            prob_null = 1

        return prob_null, len(pairs)

    def num_trial_pairs(self):  #Reindexed with geneList
        genecount = pd.DataFrame.from_dict(Counter(self.curr_cell_df.index) , orient='index').reindex(self.geneList).fillna(0)
        num_trial_pairs = genecount.dot(genecount.T)   #n1*n2
        return genecount, num_trial_pairs

    def obs_trial_pairs(self): #debug for large d
        point_tree = cKDTree(self.curr_cell_df)
        pairs = point_tree.query_pairs(self.dist_thresh)
        pairs = [(self.curr_cell_df.index[i], self.curr_cell_df.index[j]) for (i,j) in pairs]
        pairs = Counter(pairs)
        pairs = pd.Series(pairs).reset_index()
        if len(pairs):
            obs_df = pd.pivot_table(pairs, index = 'level_0', columns ='level_1', values = 0).fillna(0)
            col_na_genes = [i for i in self.geneList if i not in obs_df.columns]
            row_na_genes = [i for i in self.geneList if i not in obs_df.index]
            col_na_genes = dict.fromkeys(col_na_genes, 0)
            obs_df = obs_df.assign(**col_na_genes)
            for row_na in row_na_genes:
                obs_df.loc[row_na] = 0
            obs_df = obs_df.reindex(index = self.geneList, columns = self.geneList)
        else:                                       #if no entry less than dist thresh
            print('no entry less than dist thresh, total rna', self.curr_cell_df.shape[0])
            obs_df = pd.DataFrame(0, self.geneList, self.geneList)

        arr2 = np.triu(obs_df) + np.triu(obs_df,1).T   #Making matrix symmetric
        return pd.DataFrame(arr2, self.geneList, self.geneList)



class ColocalizationModelSingleaxis(ColocalizationModel):
    def __init__(self, geneList, df, dist_thresh, z_axis = 0):
        df = df[['uID', 'geneName', 'absX', 'absY', 'absZ']]
        df = df[df.absZ == z_axis][['uID', 'geneName', 'absX', 'absY']]
        super().__init__(geneList, df,  dist_thresh)


class ColocalizationModelAllaxis():
    def __init__(self, geneList, df, dist_thresh, min_gene_counts = 10):
        self.df = df
        self.dist_thresh = dist_thresh
        self.geneList = geneList
        self.min_gene_counts = min_gene_counts
        self.z_list = np.sort(df.absZ.unique())
        self.obs, self.num_trial, self.genecount_all, self.prob_null = self.binom_params_all_axis()
        # print('prob null is ', self.prob_null)
        self.p_vals = self.compute_p_val()
        

    def binom_params_all_axis(self):
        obs_all = []
        num_trial_all = []
        prob_null_all = []
        len_prob_null_all = []
        genecount_all = []
        
        for z_axis in self.z_list:
            
            df = self.df[['uID', 'geneName', 'absX', 'absY', 'absZ']]
            df = df[df.absZ == z_axis][['uID', 'geneName', 'absX', 'absY']]
            if df.shape[0] > self.min_gene_counts:
                model_single_axis = ColocalizationModelSingleaxis(self.geneList, self.df,  self.dist_thresh, z_axis)
                obs_all.append(model_single_axis.obs_df.values)
                num_trial_all.append(model_single_axis.num_trial_df.values)
                genecount_all.append(model_single_axis.genecount.iloc[:,0].values)  #df series doesn't return values with iloc directly
                prob_null, len_prob_null = model_single_axis.prob_null_hypothesis()
                prob_null_all.append(prob_null)
                len_prob_null_all.append(len_prob_null)
                # print('all prob', prob_null_all)
                # print('total count slices', len_prob_null_all)

        
        if len(prob_null_all):
            obs_all = np.array(obs_all).sum(axis = 0)
            num_trial_all = np.array(num_trial_all).sum(axis = 0)
            genecount_all = np.array(genecount_all).sum(axis = 0)
            prob_null_all = sum(x * y for x, y in zip(prob_null_all, len_prob_null_all))/sum(len_prob_null_all)
        else:       #useless cells, edge case
            obs_all = np.zeros((len(self.geneList), len(self.geneList)))
            num_trial_all = np.zeros((len(self.geneList), len(self.geneList)))
            genecount_all = np.zeros(len(self.geneList))
            prob_null_all = 1      
        return obs_all, num_trial_all, genecount_all, prob_null_all

    def compute_p_val(self):
        p_vals = np.ones((len(self.geneList), len(self.geneList)))
        for i in range(self.obs.shape[0]):
            for j in range(i,self.obs.shape[1]):
                p_vals[i,j] = binom_test(self.obs[i,j], self.num_trial[i,j], \
                self.prob_null, alternative = 'greater')
                p_vals[j,i] = p_vals[i,j]

        # p_vals_df = pd.DataFrame(p_vals, index = self.geneList, columns = self.geneList)
        return p_vals#_df

    






