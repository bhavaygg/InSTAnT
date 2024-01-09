import gc
import pickle
import timeit
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import multiprocessing as mp
import matplotlib.pyplot as plt
from InSTAnT.poibin import PoiBin
from scipy.spatial import cKDTree
from collections import Counter
from scipy.stats import hypergeom
from scipy.stats import binom_test
from scipy.optimize import minimize
from sklearn.preprocessing import LabelEncoder
from scipy.stats import multivariate_hypergeom
from InSTAnT.poisson_binomial import PoissonBinomial
import anndata as ad

class ConditionalGlobalColocalization():
    '''
        Performs conditional global colocalization
        Requires output from ProximalPairs()
        Arguments: 
            - all_pvals: (Array) Gene-Gene pairwise pvalues for all cells calcuated using ProximalPairs().
            - transcript_count: (Array) Expression count of each gene across each cell.
            - alpha_cellwise: (Float) pvalue signifcance threshold (>alpha_cellwise are converted to 1).
            - min_transcript: (Float) Gene expression lower threshold.
            - show_det_pairs: (Not used)
            - high_precision: (Boolean) High precision pvalue. Expect longer computer.
            - threads: (Integer) Number of threads to use.
    '''
    def __init__(self, all_pvals, transcript_count, alpha_cellwise = 0.01, min_transcript = 0, show_det_pairs = 0, high_precision = False, threads = 1, precision_mode = 'high'):
        global all_cell_pval
        if precision_mode == 'high':
            self.p_mode = 1e-64
        else:
            self.p_mode = 1e-16
        all_cell_pval = all_pvals
        self.alpha, self.show_det_pairs = alpha_cellwise, show_det_pairs
        self.transcript_count = transcript_count   #num_cells x num_genes
        self.min_transcript = min_transcript
        self.high_precision = high_precision
        self.threads = threads
        self.binarize_adj_matrix()
        del all_cell_pval
        self.global_weight_pair()
        print("Global Colocalization initialized ..")

    def binarize_adj_matrix(self):
        '''
        Convert all pvals > alpha to 1
        '''
        edge_all  = np.zeros(all_cell_pval.shape)
        edge_all[all_cell_pval < self.alpha] = 1
        global edge, obs
        edge = edge_all
        obs = edge_all.sum(axis=0) #pairwise colocalisation
        self.num_cells = edge_all.shape[0]
        self.num_genes = edge_all.shape[1]
        print('Number of cells: %d, Number of genes: %d' %(self.num_cells, self.num_genes))

    def global_weight_pair(self):
        global_genecount = edge.sum(axis = (0,1)).reshape([-1,1])
        weight_pair = np.matmul(global_genecount, np.transpose(global_genecount))
        self.gene_pair_weight = weight_pair

    def prob_cells_coloc(self):
        #global prob_pairwise_all_cells
        prob_pairwise_all_cells = np.ones((self.num_cells, self.num_genes, self.num_genes))
        for i in range(self.num_cells):
            curr_edge = edge[i].copy()
            curr_edge = curr_edge[np.triu_indices(curr_edge.shape[0])]
            #zero count genes weight
            curr_weight_pair = self.gene_pair_weight.copy() 
            curr_weight_pair[self.transcript_count[i] <= self.min_transcript] = 0 #removing lower count genes
            temp = curr_weight_pair[np.triu_indices(curr_weight_pair.shape[0])]
            # assert temp.sum() != 0
            curr_weight_pair = curr_weight_pair/(temp.sum() + self.p_mode)   #adding small number in denominator to avoid inf
            prob_pairwise_all_cells[i] = 1 - np.power((1 - curr_weight_pair), curr_edge.sum())
        return prob_pairwise_all_cells
        
    def _compute_pval(self, args):
        i, j = args[0], args[1]
        pp_all_cells_local = np.frombuffer(pp_all_cells).reshape(pp_all_cells_shape)
        p_ij = []
        for cell_id in range(self.num_cells):
            if self.transcript_count[cell_id, i] > self.min_transcript and self.transcript_count[cell_id, j] > self.min_transcript:
                p_ij.append(pp_all_cells_local[cell_id,i,j])
        pb = PoiBin(p_ij)
        pval_curr = pb.pval(obs[i,j].astype(int))
        if self.high_precision and pval_curr < 1e-13:  #Poibin has numerical error in low p val region
            pb_highres = PoissonBinomial(p_ij)
            pval_curr = pb_highres.pval(obs[i,j].astype(int))
        return i, j, pval_curr, np.sum(p_ij)
    
    def _initializer_func(self, X, X_shape):
        global pp_all_cells, pp_all_cells_shape
        pp_all_cells = X
        pp_all_cells_shape = X_shape

    def global_colocalization(self):
        start = timeit.default_timer()
        prob_pairwise_all_cells = self.prob_cells_coloc()
        share = mp.RawArray('d', self.num_cells*self.num_genes*self.num_genes)
        share_np = np.frombuffer(share).reshape(prob_pairwise_all_cells.shape)
        # Copy data to our shared array.
        np.copyto(share_np, prob_pairwise_all_cells)
        with mp.Pool(processes=self.threads, initializer=self._initializer_func, initargs=(prob_pairwise_all_cells, prob_pairwise_all_cells.shape), maxtasksperchild = 1) as pool:
            results = pool.map(self._compute_pval, [[i, j] for i in range(self.num_genes) for j in range(i, self.num_genes)])
        coloc_matrix = np.ones((self.num_genes, self.num_genes))
        expected_coloc = np.zeros((self.num_genes, self.num_genes))
        for genepair_ij_result in results:
            i, j = genepair_ij_result[0], genepair_ij_result[1]
            coloc_matrix[i,j], expected_coloc[i,j] = genepair_ij_result[2], genepair_ij_result[3]
            coloc_matrix[j,i], expected_coloc[j,i] = coloc_matrix[i,j], expected_coloc[i,j]
        if self.high_precision:
            mode = "High"
        else:
            mode = "Low"
        print(f"{mode} Precision Global Colocalization Time: {round(timeit.default_timer() - start, 2)} seconds")
        return coloc_matrix, expected_coloc
    
    def _compute_pval_serial(self, args):
        i, j = args[0], args[1]
        p_ij = []
        for cell_id in range(self.num_cells):
            if self.transcript_count[cell_id, i] > self.min_transcript and self.transcript_count[cell_id, j] > self.min_transcript:
                p_ij.append(self.pp_all_cells_local[cell_id,i,j])
        pb = PoiBin(p_ij)
        pval_curr = pb.pval(obs[i,j].astype(int))
        if self.high_precision and pval_curr < 1e-13:  #Poibin has numerical error in low p val region
            pb_highres = PoissonBinomial(p_ij)
            pval_curr = pb_highres.pval(obs[i,j].astype(int))
        return i, j, pval_curr, np.sum(p_ij)
    
    def global_colocalization_serial(self):
        start = timeit.default_timer()
        self.pp_all_cells_local = self.prob_cells_coloc()
        coloc_matrix = np.ones((self.num_genes, self.num_genes))
        expected_coloc = np.zeros((self.num_genes, self.num_genes))
        for i in range(self.num_genes):
            for j in range(i, self.num_genes):
                _, _, pv, exp_pv = self._compute_pval_serial([i, j])
                coloc_matrix[i,j], expected_coloc[i,j] = pv, exp_pv
                coloc_matrix[j,i], expected_coloc[j,i] = coloc_matrix[i,j], expected_coloc[i,j]
        if self.high_precision:
            mode = "High"
        else:
            mode = "Low"
        print(f"{mode} Precision Global Colocalization Time: {round(timeit.default_timer() - start, 2)} seconds")
        return coloc_matrix, expected_coloc

class ProximalPairs():
    '''
        Calculates proximal pairs for all gene-gene pairs for a given cell.
        Arguments: 
            - geneList: (Array) List of genes.
            - df_loc: (DataFrame) Coordinates of the gene transcript.
            - distance_threshold: (Integer) distance threshold at which to consider 2 genes proximal.
            - mode: (Not used).
    '''
    def __init__(self, geneList, df_loc, distance_threshold, mode="global", genes = None):
        self.geneList = geneList

        self.distance_threshold = distance_threshold
        if mode=='global':
            position_matrix_local = np.frombuffer(position_all_cells).reshape(position_all_cells_shape)
            cell_positions = np.where(position_matrix_local[:,2] == df_loc)
            self.curr_cell_df = position_matrix_local[cell_positions, 0:2][0]
            self.genes = genes[cell_positions]
        elif mode == "global_3d":
            position_matrix_local = np.frombuffer(position_all_cells).reshape(position_all_cells_shape)
            cell_positions = np.where(position_matrix_local[:,3] == df_loc)
            self.curr_cell_df = position_matrix_local[cell_positions, 0:3][0]
            self.genes = genes[cell_positions]
        else:
            self.curr_cell_df = df_loc
            self.genes = genes
        self.genecount, self.num_trial = self.num_trial_pairs()
        self.prob_null = self.prob_null_hypothesis()  
        self.obs = self.obs_trial_pairs()
        self.p_vals = self.compute_p_val()
    
    def compute_p_val(self):
        p_vals = np.ones((len(self.geneList), len(self.geneList)))
        for i in range(self.obs.shape[0]):
            for j in range(i,self.obs.shape[1]):
                p_vals[i,j] = binom_test(self.obs[i,j], self.num_trial[i,j], \
                self.prob_null, alternative = 'greater')
                p_vals[j,i] = p_vals[i,j]
        return p_vals

    def prob_null_hypothesis(self):
        '''
        Using a KDTree, find all pairs of genes in a cell that are under distance `distance_threshold`
        '''
        self.point_tree = cKDTree(self.curr_cell_df)
        self.pairs = self.point_tree.query_pairs(self.distance_threshold)
        if len(self.pairs):   #Later change the condition to min gene count and other heuristic
            prob_null = len(self.pairs)*2/(self.curr_cell_df.shape[0]*(self.curr_cell_df.shape[0]-1))
        else:
            prob_null = 1
        return prob_null

    def num_trial_pairs(self):  #Reindexed with geneList
        '''
        Count for each gene pair, the number of times they are proximal given the `distance_threshold`
        '''
        genecount = pd.DataFrame.from_dict(Counter(self.genes) , orient='index').reindex(self.geneList).fillna(0)
        num_trial_pairs = np.dot(genecount, genecount.T)   #n1*n2  #using numpy here now
        return genecount, num_trial_pairs

    def obs_trial_pairs(self): #debug for large d
        '''
        After counting all gene pairs within the distance threshold, create a pivot table to create a 2d represtation
        and fill in the missing genes.
        '''
        pairs = [(self.genes[i], self.genes[j]) for (i,j) in self.pairs]
        pairs = Counter(pairs)
        pairs = pd.Series(pairs).reset_index()
        if len(pairs):
            obs_df = pd.pivot_table(pairs, index = 'level_0', columns ='level_1', values = 0).fillna(0)
            col_na_genes = list(set(self.geneList).difference(obs_df.columns))
            row_na_genes = list(set(self.geneList).difference(obs_df.index))
            obs_df.loc[ : , col_na_genes] = 0
            for row_na in row_na_genes:
                obs_df.loc[row_na] = 0
            obs_df = obs_df.reindex(index = self.geneList, columns = self.geneList)
        else:                                       #if no entry less than dist thresh
            obs_df = pd.DataFrame(0, self.geneList, self.geneList)
        temp= obs_df.values
        arr2 = temp + temp.T - temp*np.identity(len(temp))
        return arr2
    
class ProximalPairs3D():
    '''
        Calculates proximal pairs for all gene-gene pairs for a given cell.
        Arguments: 
            - geneList: (Array) List of genes.
            - df_loc: (DataFrame) Coordinates of the gene transcript.
            - distance_threshold: (Integer) distance threshold at which to consider 2 genes proximal.
            - mode: (Not used).
    '''
    def __init__(self, geneList, df_loc, distance_threshold, min_genecount = 10, cell_id = None, genes = None):
        global df_cell, obs, p_vals
        position_matrix_local = np.frombuffer(position_all_cells).reshape(position_all_cells_shape)
        cell_positions = np.where(position_matrix_local[:,3] == df_loc)
        self.df_cell = position_matrix_local[cell_positions, 0:3][0]
        self.genes = genes[cell_positions]
        self.distance_threshold = distance_threshold
        self.geneList = geneList
        self.min_genecount = min_genecount
        self.cell_id = cell_id
        self.z_list = np.unique(self.df_cell[:, 2])
        obs, self.num_trial, self.genecount, self.prob_null = self.binom_params_all_axis()
        self.p_vals = self.compute_p_val()

    def binom_params_all_axis(self):
        obs_all = []
        num_trial_all = []
        prob_null_all = []
        len_prob_null_all = []
        genecount_all = []
        for z_axis in self.z_list:
            z_positions = np.where(self.df_cell[:,2] == z_axis)
            df_z = self.df_cell[z_positions, 0:2][0]
            genes_z = self.genes[z_positions]
            if df_z.shape[0] > self.min_genecount:
                model_single_axis = ProximalPairs(self.geneList, df_z, self.distance_threshold, mode = "3D", genes = genes_z)
                obs_all.append(model_single_axis.obs)
                num_trial_all.append(model_single_axis.num_trial)
                genecount_all.append(model_single_axis.genecount.iloc[:,0].values)  #df series doesn't return values with iloc directly
                prob_null = model_single_axis.prob_null
                prob_null_all.append(prob_null)
                len_prob_null_all.append(len(model_single_axis.pairs))

        if len(prob_null_all) > 0 and sum(len_prob_null_all) > 0:
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
        for i in range(obs.shape[0]):
            for j in range(i, obs.shape[1]):
                p_vals[i,j] = binom_test(obs[i,j], self.num_trial[i,j], \
                self.prob_null, alternative = 'greater')
                p_vals[j,i] = p_vals[i,j]
        return p_vals


class SpatialModulation:
    def __init__(self, edge_all_cells, pos_cells, dist, geneList, file_name, min_num_neighbor = 1, threads = 1, adata_filename = None):
        self.edge_all_cells = edge_all_cells
        self.num_genes = edge_all_cells.shape[1]
        self.pos_cells = pos_cells
        self.dist, self.min_num_neighbor = dist, min_num_neighbor
        self.threads = threads
        self.geneList = geneList
        self.filename = file_name
        self.adata_filename = adata_filename
        self.estimate_prob()
        self.optimize_all_pair_h1()
        self.unstack_df_llr()
        if self.adata_filename == None:
            self.adata_filename = "null"

    def estimate_prob(self, show_neighbor_hist = 0):
        print('Estimating Global and Local Probabilities')
        num_cells = self.pos_cells.shape[0]
        prob_global = self.edge_all_cells.sum(axis = 0)/num_cells
        prob_neighbors = np.zeros((num_cells, self.num_genes, self.num_genes))
        point_tree = cKDTree(self.pos_cells)

        num_neighbors = []
        for cell_index in range(num_cells):
            curr_cell_loc = self.pos_cells.iloc[cell_index,:]
            neighbors_index = point_tree.query_ball_point(curr_cell_loc, self.dist)
            neighbors_index.remove(cell_index)
            num_neighbors.append(len(neighbors_index))
            if len(neighbors_index) < self.min_num_neighbor :
                curr_prob_neighbors = prob_global
            else:
                curr_prob_neighbors = (self.edge_all_cells[neighbors_index , :, :].sum(axis=0))/(len(neighbors_index)) 
            prob_neighbors[cell_index] = curr_prob_neighbors

        if show_neighbor_hist:
            plt.hist(num_neighbors)
        self.prob_global, self.prob_neighbors = prob_global, prob_neighbors

    def find_log_likelihood_pair(self, i,j, obs_edge, with_neighbor = 0):
        if with_neighbor : 
            neigh_prob = self.prob_neighbors[:,i,j]
            neigh_prob_0 = neigh_prob[obs_edge == 0]
            neigh_prob_1 = neigh_prob[obs_edge == 1]
            ll = np.log(1 - neigh_prob_0).sum() + np.log(neigh_prob_1).sum()
        else:
            p1 = self.prob_global[i,j]
            if p1 ==0:
                p1 = 1e-64
            ll = obs_edge.sum() * np.log(p1) + (obs_edge.shape[0] - obs_edge.sum()) * np.log(1 - p1)
        return ll

    def objective_function(self, x, obs_edge, neigh_prob):
        w = x[0]
        p1 = x[1]
        t = p1*(1-w) + (w)* neigh_prob  #x is weight
        t_0 = t[obs_edge == 0]
        t_1 = t[obs_edge == 1]
        log_ll = (np.log(1 - t_0)).sum() + np.log(t_1).sum()
        return -log_ll
    
    def _parallelize_optimization(self, args):
        i, j, obs_edge, neigh_prob = args[0], args[1], args[2], args[3]
        return i, j, self.find_log_likelihood_pair(i,j, obs_edge, with_neighbor = 0), minimize(self.objective_function, x0 = [0.1, self.prob_global[i,j]], args=(obs_edge, neigh_prob), bounds= ((0.00001,0.99999),(0.00001,0.99999)))
    
    def optimize_all_pair_h1(self):
        print(f"Running Spatial Modulation now on {self.threads} threads for, {self.pos_cells.shape[0]} cells, {self.num_genes} genes ..")
        pool = mp.Pool(processes = self.threads)
        parallel_arguments = []
        for i in range(self.num_genes):
            for j in range(i, self.num_genes):
                parallel_arguments.append([i, j, self.edge_all_cells[:,i,j], self.prob_neighbors[:,i,j]])
        del self.edge_all_cells, self.prob_neighbors
        results = pool.map(self._parallelize_optimization, parallel_arguments)
        ll_h0_all = np.zeros((self.num_genes, self.num_genes))
        ll_h1_all = np.zeros(( self.num_genes, self.num_genes))
        w_h1_all = np.zeros((self.num_genes, self.num_genes))
        p_h1_all = np.zeros((self.num_genes, self.num_genes))
        for result in results:
            i, j, ll_ij, res = result[0], result[1], result[2], result[3]
            [w, p_h1] = res.x
            w_h1_all[i,j] = w
            p_h1_all[i,j] = p_h1
            ll_h1_all[i,j]  = -res.fun
            ll_h1_all[j,i] = ll_h1_all[i,j]
            ll_h0_all[i,j] = ll_ij
        self.w_h1_all = w_h1_all
        self.p_h1_all = p_h1_all
        self.ll_h0_all, self.ll_h1_all = ll_h0_all, ll_h1_all
    
    def unstack_df_llr(self):
        print("Saving ..")
        llr = self.ll_h1_all - self.ll_h0_all
        llr = llr[np.triu_indices(self.num_genes)]
        self.w_h1_all = self.w_h1_all[np.triu_indices(self.num_genes)]
        self.p_h1_all = self.p_h1_all[np.triu_indices(self.num_genes)]
        self.prob_global = self.prob_global[np.triu_indices(self.num_genes)]
        gene_id1 = [self.geneList[i] for i in np.triu_indices(self.num_genes)[0]]
        gene_id2 = [self.geneList[i] for i in np.triu_indices(self.num_genes)[1]]
        data = {'gene_id1': gene_id1, 'gene_id2': gene_id2, 'llr': llr, 'w_h1': self.w_h1_all, 'p_g_h1': self.p_h1_all, 'p_g_h0' : self.prob_global}
        pairwise_p_val_df = pd.DataFrame(data)
        pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'].astype(str) + ', ' + pairwise_p_val_df['gene_id2'].astype(str)
        pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
        pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['llr'],ascending=False)
        if Path(self.adata_filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.adata_filename)
            adata.uns["spatial_modulation"] = pairwise_p_val_df
            adata.write(self.adata_filename, compression="gzip")
        if self.filename != None:
            pairwise_p_val_df.to_excel(self.filename)

class DifferentialColocalization():
    def __init__(self, all_p_vals, genecount, geneList, ct_spec_indices, rest_indices, cell_type, file_location, alpha = 0.01, threads = 1, filename = None):
        self.threads = threads
        self.cell_type = cell_type
        self.ct_spec_indices, self.rest_indices = ct_spec_indices, rest_indices
        self.all_p_vals = all_p_vals[ct_spec_indices+rest_indices, :, :]
        global gc
        gc = genecount[ct_spec_indices+rest_indices, :]
        self.geneList, self.alpha = geneList, alpha
        self.binarize_adj_matrix()
        self.all_pairs_cond_df, self.all_pairs_uncond_df, self.all_pairs_g1_df = self.compute_all_pairs()
        if file_location != None:
            self.all_pairs_cond_df.to_csv(Path(file_location) / f"{cell_type}_conditional.csv")
            self.all_pairs_uncond_df.to_csv(Path(file_location) / f"{cell_type}_unconditional.csv")
            self.all_pairs_g1_df.to_csv(Path(file_location) / f"{cell_type}_genemarkers.csv")
        self.obs = obs
        self.filename = filename
        if self.filename == None:
            self.filename = "null"
        self._save_unstacked_pvals(file_location, cell_type)

    def _save_unstacked_pvals(self, file_location, cell_type, min_transcript = 0):
        present_cells = self._num_present_cells(min_transcript)
        unstacked_global_pvals = self._unstack_df_both(present_cells)
        if Path(self.filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.filename)
            adata.uns['differential_colocalization'][f"{cell_type}"] = unstacked_global_pvals
            adata.write(self.filename, compression="gzip")
        if file_location != None:
            unstacked_global_pvals.to_excel(Path(file_location) / f"{cell_type}_unstacked.xlsx") 

    def _unstack_df_both(self, present_cells, second_df_name="Unconditional"):
        num_genes = self.all_pairs_uncond_df.values.shape[0]
        agg_coloc = self.all_pairs_uncond_df.values[np.triu_indices(num_genes)]
        cond_agg_coloc = self.all_pairs_cond_df.values[:,:][np.triu_indices(num_genes)]
        obs = self.obs[np.triu_indices(num_genes)]
        present_cells = present_cells[np.triu_indices(num_genes)]
        gene_id1 = [self.geneList[i] for i in np.triu_indices(num_genes)[0]]
        gene_id2 = [self.geneList[i] for i in np.triu_indices(num_genes)[1]]
        data = {'gene_id1': gene_id1, 'gene_id2': gene_id2, 'ct_cond': cond_agg_coloc, second_df_name: agg_coloc, f'Conditional cells(Threshold<{str(0.01)})': obs, 'Present cells' : present_cells}
        pairwise_p_val_df = pd.DataFrame(data)
        pairwise_p_val_df['frac_cells'] = pairwise_p_val_df[f'Conditional cells(Threshold<{str(0.01)})']/pairwise_p_val_df['Present cells']
        pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'].astype(str) + ', ' + pairwise_p_val_df['gene_id2'].astype(str)
        pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
        pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['ct_cond'])
        return pairwise_p_val_df

    def _num_present_cells(self, min_transcript):
        indicator = np.zeros(gc.shape)
        indicator[gc > min_transcript] = 1
        ind_n1n2 = np.matmul(indicator.reshape(indicator.shape[0],indicator.shape[1],1),indicator.reshape(indicator.shape[0], 1,indicator.shape[1]))
        present_cells = ind_n1n2.sum(axis=0)
        return present_cells

    def binarize_adj_matrix(self):
        edge_all  = np.zeros(self.all_p_vals.shape)
        edge_all[self.all_p_vals < self.alpha] = 1
        global edge, obs
        edge = edge_all
        obs = edge_all.sum(axis=0)
        self.num_cells = edge_all.shape[0]
        self.num_genes = edge_all.shape[1]
        del self.all_p_vals
        
    def all_set_intersection(self, i, j):
        M = range(len(self.ct_spec_indices)) #Cell type index
        O = self.colocalized_set(i,j)
        E = self.high_exp_set(i,len(M))
        gamma = len([k for k in M if k in E])
        lambd = len([k for k in M if k in O])
        alpha = len([k for k in O if k in E])
        return gamma,lambd, alpha, len(M), len(E), len(O)
    
    def compute_cond_p_val(self,i,j):
        gamma,lambd, alpha,m, n1, n2 = self.all_set_intersection(i,j)
        N = self.num_cells 
        p_val_num = 0
        for k in range(lambd, min(m,n2)+1):
            for beta in range(k+1):
                x = [beta, k-beta, alpha-beta, n2-alpha-k+beta]
                r = [gamma, m-gamma, n1-gamma, N-m-n1+gamma]
                curr_pmf = multivariate_hypergeom.pmf(x=x, m=r, n=n2)
                p_val_num = p_val_num + curr_pmf
        p_val_denom = hypergeom(N, n1, n2).pmf(alpha)
        return p_val_num/p_val_denom
        
    def colocalized_set(self,i,j):
        all_obs_ij = edge[:,i,j]
        O = [k for k,item in enumerate(all_obs_ij) if item == 1]
        return O
    
    def high_exp_set(self,i,n1):
        genecount_i = gc[:,i]
        E = np.argpartition(genecount_i, -n1)[-n1:] #top n1
        return E
    
    def compute_uncond_p_val(self, i, j, isPair = True):
        gamma, lambd, alpha, m, n1, n2 = self.all_set_intersection(i,j)
        if isPair:
            p_val = hypergeom.sf(lambd-1, self.num_cells, n2, m)
        else:
            p_val = hypergeom.sf(gamma-1, self.num_cells, n1, m)
        return p_val
    
    def _compute_genemarker(self, args):
        i = args
        return i, self.compute_uncond_p_val(i, 0, isPair = False)
    
    def _compute_conditional(self, args):
        i, j = args[0], args[1]
        return i, j, self.compute_cond_p_val(i, j)
    
    def _compute_unconditional(self, args):
        i, j = args[0], args[1]
        return i, j, self.compute_uncond_p_val(i, j, isPair = True)
    
    def compute_all_pairs(self):
        print(f'Computing Cell Type Characterization for Hypergeometric for cell type - {self.cell_type}')
        pool = mp.Pool(processes = self.threads)
        print("Unconditional Started ....")
        unconditional_data = [(i, j) for i in range(self.num_genes) for j in range(i, self.num_genes)]
        unconditional_data_results = pool.map(self._compute_unconditional, unconditional_data)
        all_pairs_uncond = np.ones((self.num_genes, self.num_genes))
        for result in unconditional_data_results:
            all_pairs_uncond[result[0]][result[1]] = result[2]
            all_pairs_uncond[result[1]][result[0]] = result[2]
        del unconditional_data_results
        print("Unconditional Completed ....")

        print("Gene Marker Started ....")
        gene_marker_data = list(range(self.num_genes))
        gene_marker_data_results = pool.map(self._compute_genemarker, gene_marker_data)
        all_pairs_g1 = np.ones((self.num_genes,))
        for result in gene_marker_data_results:
            all_pairs_g1[result[0]] = result[1]
        del gene_marker_data_results
        print("Gene Marker Completed ....")

        print("Conditional Started ....")
        conditional_data = [(i, j) for i in range(self.num_genes) for j in range(self.num_genes) if all_pairs_uncond[i][j] <= 1e-3]
        conditional_data_results = pool.map(self._compute_conditional, conditional_data)
        all_pairs_cond = np.ones((self.num_genes, self.num_genes))
        for result in conditional_data_results:
            all_pairs_cond[result[0]][result[1]] = result[2]
        del conditional_data_results
        print("Conditional Completed ....")


        all_pairs_cond_df = pd.DataFrame(all_pairs_cond, index = self.geneList, columns = self.geneList)
        all_pairs_uncond_df = pd.DataFrame(all_pairs_uncond, index = self.geneList, columns = self.geneList)
        all_pairs_g1_df = pd.DataFrame(all_pairs_g1, index = self.geneList)
        return all_pairs_cond_df, all_pairs_uncond_df, all_pairs_g1_df
    
    def compute_all_pairs_serial(self):
        print(f'Computing Cell Type Characterization for Hypergeometric for cell type - {self.cell_type}')
        print("Gene Marker Started ....")
        gene_marker_data = list(range(self.num_genes))
        all_pairs_g1 = np.ones((self.num_genes,))
        for i in gene_marker_data:
            result = self._compute_genemarker(i)
            all_pairs_g1[result[0]] = result[1]
        print("Gene Marker Completed ....")

        print("Conditional Started ....")
        conditional_data = [(i, j) for i in range(self.num_genes) for j in range(self.num_genes)]
        all_pairs_cond = np.ones((self.num_genes, self.num_genes))
        for i in conditional_data:
            result = self._compute_conditional(i)
            all_pairs_cond[result[0]][result[1]] = result[2]
        print("Conditional Completed ....")

        print("Unconditional Started ....")
        unconditional_data = [(i, j) for i in range(self.num_genes) for j in range(i, self.num_genes)]
        all_pairs_uncond = np.ones((self.num_genes, self.num_genes))
        for i in unconditional_data:
            result = self._compute_unconditional(i)
            all_pairs_uncond[result[0]][result[1]] = result[2]
            all_pairs_uncond[result[1]][result[0]] = result[2]
        print("Unconditional Completed ....")

        all_pairs_cond_df = pd.DataFrame(all_pairs_cond, index = self.geneList, columns = self.geneList)
        all_pairs_uncond_df = pd.DataFrame(all_pairs_uncond, index = self.geneList, columns = self.geneList)
        all_pairs_g1_df = pd.DataFrame(all_pairs_g1, index = self.geneList)
        return all_pairs_cond_df, all_pairs_uncond_df, all_pairs_g1_df

class GeneModuleDiscovery():
    '''
    TODO #Globalcolo MAP , Frequent subgraph mining(using gspan)
    '''
    def __init__(self) -> None:
        pass

class Instant():
    '''
    Intracellular Spatial Transcriptomic Analysis Toolkit (InSTAnT).
    This class is used to calculate the proximal gene pairs in each cell and 
    used to find which gene pairs are d-cololocalised across all the cells.
        Arguments: 
            - threads: (Integer) Number of threads to use.
            - min_intensity: (Optional) (Integer) Minimum intensity for Merfish.
            - min_area: (Optional) (Integer) Minimum Area for Merfish.
    '''
    def __init__(self, threads = 1, precision_mode = "high", min_intensity = 0, min_area = 0):
        self.min_intensity, self.min_area = min_intensity, min_area
        self.threads = threads
        if precision_mode:
            self.precision_mode = np.float64
        else:
            self.precision_mode = np.float16

    def load_preprocessed_data(self, data, inNucleus = False):
        '''
        Load preprocessed data. Data can either be a '.csv' or a '.h5ad' file.
        If AnnData file is passed, the csv file is expected in `adata.uns['transcripts']`
        The csv should have should have the following columns - 
            ['gene', 'absX', 'absY', 'uID'] with 'gene' being the 1st column
            Arguments: 
                - data: (String) Path to dataframe in the required format.
        '''
        self.filename = data
        if Path(self.filename).suffix.lower() == ".csv":
            self.df = pd.read_csv(data, index_col=0).sort_index()
        elif Path(self.filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.filename)
            self.df = adata.uns['transcripts']
            self.df = self.df.set_index('gene').sort_index()
        else:
            raise("Input File Format Incorrect")
        if inNucleus:
            self.df = self.df[self.df.inNucleus == 1]
        self.geneList = self.df.index.unique()
        print("Loaded Data. Number of Transcripts: ", len(self.df), ", Number of Genes: ", len(self.geneList), 
              ", Number of Cells: ", len(self.df.uID.unique()))
    
    def load_preprocessed_data_randomize(self, data):
        '''
        Load preprocessed data and randomize the genes (to establish baselines). Data should have the following columns - 
            ['gene', 'absX', 'absY', 'uID'] with 'gene' being the 1st column
            Arguments: 
                - data: (String) Path to dataframe in the required format.
        '''
        self.filename = data
        self.df = pd.read_csv(data, index_col=0)
        self.df.index = np.random.permutation(self.df.index.values)
        self.geneList = self.df.index.unique()
    
    def load_preprocessed_data_xenium(self, data, cluster_file, clusters, randomize = False):
        '''
        Load preprocessed data. Data should have the following columns - 
            ['gene', 'absX', 'absY', 'uID'] with 'gene' being the 1st column
            Arguments: 
                - data: (String) Path to dataframe in the required format.
        '''
        self.filename = data
        self.df = pd.read_csv(data, index_col=0).sort_index()
        self.geneList = self.df.index.unique()
        print("Loaded Data. Number of Transcripts: ", len(self.df), ", Number of Genes: ", len(self.geneList), 
              ", Number of Cells: ", len(self.df.uID.unique()))
        cluster_df = pd.read_csv(cluster_file)
        cluster_df.rename(columns={'Barcode': 'uID'}, inplace=True)
        cluster_df.uID = cluster_df.uID - 1
        cluster_df.set_index('uID', inplace=True)
        self.df = self.df.join(cluster_df, on='uID')
        self.df = self.df[self.df.Cluster.isin(clusters)]
        if randomize:
            self.df.index = np.random.permutation(self.df.index.values)
            self.geneList = self.df.index.unique()
    
    def load_preprocessed_data_merfish(self, data, min_intensity = 0, min_area = 0):
        '''
        Load preprocessed data. Data should have the following columns - 
            ['gene', 'absX', 'absY', 'uID'] with 'gene' being the 1st column
            Arguments: 
                - data: (String) Path to dataframe in the required format.
        '''
        self.filename = data
        self.df = pd.read_csv(data, index_col=0)
        self.df.dropna(subset=['uID'], inplace=True)
        self.df = self.df[self.df.Total_brightness > min_intensity] #moved up
        self.df = self.df[self.df.Area >= min_area] #moved up
        self.geneList = self.df.index.unique()
        print("Loaded Data. Number of Transcripts: ", len(self.df), ", Number of Genes: ", len(self.geneList), 
              ", Number of Cells: ", len(self.df.uID.unique()))
        
    def load_preprocessed_data_filter(self, data, threshold = 100):
        '''
        Load preprocessed data and filter cells to ensure that minimum threshold transcripts are present.
        Data should have the following columns - 
        ['gene', 'absX', 'absY', 'uID'] with 'gene' being the 1st column
            Arguments: 
                - data: (String) Path to dataframe in the required format.
                - threshold: (Integer) Minimum number of transcripts in each cell.
        '''
        self.df = pd.read_csv(data, index_col=0)
        threshold_cells = self.df.groupby('uID').size()
        threshold_cells = threshold_cells[threshold_cells > threshold].index.values
        self.df = self.df.loc[self.df.uID.isin(threshold_cells)]
        cell_ids = self.df.uID.unique()
        self.geneList = self.df.index.unique()
    
    def load_preprocessed_data_filter_genes(self, data, threshold = 100):
        self.df = pd.read_csv(data, index_col=0)
        threshold_cells = self.df.groupby('uID').size()
        threshold_cells = threshold_cells[threshold_cells > threshold].index.values
        self.df = self.df.loc[self.df.uID.isin(threshold_cells)]
        threshold_genes = self.df.groupby('gene').size().sort_values(ascending=False).index.values[:100]
        self.df = self.df.loc[self.df.gene.isin(threshold_genes)]
        cell_ids = self.df.uID.unique()
        self.geneList = self.df.index.unique()
        print("Unique Z", len(self.df.absZ.unique()))
        self.df.hist(column="absZ").get_figure().savefig(f'/storage/coda1/p-ssinha338/0/shared/InSTAnT/z.png')
    
    def load_preprocessed_data_random(self, data, num_cells = 5000):
        self.df = pd.read_csv(data, index_col=0)
        cell_ids = self.df.uID.unique()
        random_cell_ids = np.random.choice(cell_ids, num_cells, replace=False)
        self.df = self.df.loc[self.df.uID.isin(random_cell_ids)]
        self.geneList = self.df.index.unique()
    
    def load_preprocessed_data_randomize(self, data):
        '''
        Load preprocessed data and randomize the genes (to establish baselines). Data should have the following columns - 
            ['gene', 'absX', 'absY', 'uID'] with 'gene' being the 1st column
            Arguments: 
                - data: (String) Path to dataframe in the required format.
        '''
        self.filename = data
        self.df = pd.read_csv(data, index_col=0)
        self.df.index = np.random.permutation(self.df.index.values)
        self.df.index.name = 'gene'
        self.geneList = self.df.index.unique()
    
    def load_preprocessed_data_filter(self, data, threshold = 100):
        '''
        Load preprocessed data and filter cells to ensure that minimum threshold transcripts are present.
        Data should have the following columns - 
        ['gene', 'absX', 'absY', 'uID'] with 'gene' being the 1st column
            Arguments: 
                - data: (String) Path to dataframe in the required format.
                - threshold: (Integer) Minimum number of transcripts in each cell.
        '''
        self.df = pd.read_csv(data, index_col=0)
        threshold_cells = self.df.groupby('uID').size()
        threshold_cells = threshold_cells[threshold_cells > threshold].index.values
        self.df = self.df.loc[self.df.uID.isin(threshold_cells)]
        cell_ids = self.df.uID.unique()
        self.geneList = self.df.index.unique()
    
    def load_preprocessed_data_filter_genes(self, data, threshold = 100, random=False):
        self.df = pd.read_csv(data, index_col=0)
        threshold_cells = self.df.groupby('uID').size()
        threshold_cells = threshold_cells[threshold_cells > threshold].index.values
        self.df = self.df.loc[self.df.uID.isin(threshold_cells)]
        threshold_genes = self.df.reset_index().groupby('gene').size().sort_values(ascending=False).index.values[:100]
        self.df = self.df.loc[self.df.index.isin(threshold_genes)]
        cell_ids = self.df.uID.unique()
        self.geneList = self.df.index.unique()
    
    def load_preprocessed_data_random(self, data, num_cells = 5000):
        self.df = pd.read_csv(data, index_col=0)
        cell_ids = self.df.uID.unique()
        random_cell_ids = np.random.choice(cell_ids, num_cells, replace=False)
        self.df = self.df.loc[self.df.uID.isin(random_cell_ids)]
        self.geneList = self.df.index.unique()

    def preprocess_and_load_data(self, expression_data, barcode_data):
        '''
        Data Loading and preprocessing for Merfish data.
            Arguments: 
                - expression_data: (String) Path to expression data csv file.
                - barcode_data: (String) Path to expression data csv file.
        '''
        self.expression_data = expression_data
        self.barcode_data = barcode_data
        df = pd.read_csv(self.expression_data)
        codebook = self._load_codebook()
        df = df[df.normalized_intensity > self.min_intensity] #moved up
        df = df[df.area >= self.min_area] #moved up
        df['geneName'] = df['barcode'].apply(lambda x: codebook.loc[x,'name'])
        self.geneList = df.geneName.unique()
        df = df.rename(columns={'cell_id': 'uID', 'abs_x': 'absX', 'abs_y': 'absY', 'geneName': 'gene'})
        df = df.drop(['barcode', 'area','is_exact', 'normalized_intensity'],axis=1)  
        self.df = df.set_index('gene')
    
    def _load_codebook(self):
        codebook = pd.read_csv(self.barcode_data, converters={'bit_barcode': lambda x: int(str(x)[::-1],2)})
        codebook = codebook[['name', 'bit_barcode']]
        codebook.rename(columns={'bit_barcode': 'barcode'}, inplace=True)
        codebook = codebook.set_index('barcode')
        return codebook

    def save_cell_id_list(self, f=None):
        if f is not None:
            with open(f, 'wb') as f:
                pickle.dump(np.unique(self.df.uID), f)
        return np.unique(self.df.uID)
    
    def _save_pval_matrix(self, filename):
        with open(filename, 'wb') as fp:
            pickle.dump(self.all_pval, fp)
    
    def load_pval_matrix(self, filename):
        with open(filename, 'rb') as fp:
            self.all_pval = pickle.load(fp)
    
    def _save_gene_count(self, filename):
        with open(filename, 'wb') as fp:
            pickle.dump(self.all_gene_count, fp)
    
    def load_gene_count(self, filename):
        with open(filename, 'rb') as fp:
            self.all_gene_count = pickle.load(fp)
    
    def save_gene_list(self, filename):
        with open(filename, 'wb') as fp:
            pickle.dump(self.geneList, fp)
    
    def load_gene_list(self, filename):
        with open(filename, 'rb') as fp:
            self.geneList = pickle.load(fp)
    
    def _calculate_ProximalPairs(self, args):
        cell_num, cell_id, genes = args[0], args[1], args[2]
        pp_model = ProximalPairs(geneList = self.geneList, df_loc = cell_id, distance_threshold = self.distance_threshold, genes = genes) #df_batch.loc[df_batch.uID == cell_id][['absX', 'absY']]
        return cell_num, pp_model.p_vals, pp_model.genecount.values.reshape(len(self.geneList))
    
    def _initializer_func_pp(self, X, X_shape):
        global position_all_cells, position_all_cells_shape
        position_all_cells = X
        position_all_cells_shape = X_shape

    def run_ProximalPairs(self, distance_threshold, min_genecount, pval_matrix_name = None, gene_count_name = None, randomize = False):
        '''
        Calculates Proximal pairs for each gene pair for each input cell. Runs the ProximalPairs() class for each cell
        and generates a 2d (num_genes, num_genes) matrix whose index (i,j) represents the p-value significance
        of whether gene i and gene j are proximal gene pairs in that cell.
            Arguments: 
                - distance_threshold: (Integer) distance threshold at which to consider 2 genes proximal.
                - min_genecount: (integer) Minimum number of transcripts in each cell.
                - pval_matrix_name: (String) (Optional) if provided saves pvalue matrix using pickle at the input path.
                - gene_count_name: (String) (Optional) if provided saves gene expression count matrix using pickle at the input path.
        '''
        self.distance_threshold = distance_threshold
        cell_ids = np.unique(self.df.uID)
        num_cells = len(cell_ids)
        num_genes = len(self.geneList)
        print(f"Initialised PP now on {self.threads} threads")
        print("Number of cells: ", num_cells, ", Number of Genes: ", num_genes)
        start = timeit.default_timer()
        valid_cell_data = []
        num_transcripts = 0
        celllabel_encoder = LabelEncoder()
        cell_labels = celllabel_encoder.fit_transform(self.df['uID'])
        self.df['uID_encoded'] = cell_labels
        genes = self.df.index
        ids, counts = np.unique(self.df.uID_encoded.values, return_counts=True)
        for n, cell_id in enumerate(ids):
            num_transcripts += counts[n]
            if counts[n] > min_genecount:
                valid_cell_data.append([n, cell_id, genes])
            else:
                print(f"min genecount less than {min_genecount} for cell id {cell_ids[n]} across all Z, Skipping ...")
        df_copy = self.df.copy()
        del self.df
        position_matrix = df_copy[['absX', 'absY', 'uID_encoded']].to_numpy().copy(order='C')
        share_pp = mp.RawArray('d', len(position_matrix)*3)
        share_pp_np = np.frombuffer(share_pp).reshape(position_matrix.shape)
        np.copyto(share_pp_np, position_matrix)
        check = timeit.default_timer()
        pool = mp.Pool(processes = self.threads)
        with mp.Pool(processes=self.threads, initializer=self._initializer_func_pp, initargs=(position_matrix, position_matrix.shape), maxtasksperchild = 1) as pool:
            if randomize:
                print(f"Running randomized PP now on {self.threads} threads for, {len(valid_cell_data)} cells, {num_transcripts} transcripts")
                results = pool.map(self._calculate_ProximalPairs_random, valid_cell_data)
            else:
                print(f"Running PP now on {self.threads} threads for, {len(valid_cell_data)} cells, {num_transcripts} transcripts")
                results = pool.map(self._calculate_ProximalPairs, valid_cell_data)
        self.all_pval = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = self.precision_mode)
        self.all_gene_count = np.zeros((num_cells, len(self.geneList)), dtype = self.precision_mode) 
        for cell_i_result in results:
            self.all_pval[cell_i_result[0]] = cell_i_result[1]
            self.all_gene_count[cell_i_result[0]] = cell_i_result[2]
        print(f"Time to complete PP : ", timeit.default_timer() - check)
        if Path(self.filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.filename)
            adata.obsm['genecount'] = self.all_gene_count
            adata.uns['pp_test_pvalues'] = self.all_pval
            adata.write(self.filename, compression="gzip")
        if pval_matrix_name:
            self._save_pval_matrix(pval_matrix_name)
        if gene_count_name:
            self._save_gene_count(gene_count_name)
        self.df = df_copy
        del df_copy
        print(f"Cell-wise Proximal Pairs Time : {round(timeit.default_timer() - start, 2)} seconds")
    
    def _calculate_ProximalPairs3D(self, args):
        cell_num, cell_id, genes = args[0], args[1], args[2]
        pp_model = ProximalPairs(geneList = self.geneList, df_loc = cell_id, distance_threshold = self.distance_threshold, 
                                 genes = genes, mode = "global_3d") #df_batch.loc[df_batch.uID == cell_id][['absX', 'absY']]
        return cell_num, pp_model.p_vals, pp_model.genecount.values.reshape(len(self.geneList))
    
    def run_ProximalPairs3D(self, distance_threshold, min_genecount, pval_matrix_name = None, gene_count_name = None, randomize = False):
        '''
        Calculates Proximal pairs for each gene pair for each input cell. Runs the ProximalPairs() class for each cell
        and generates a 2d (num_genes, num_genes) matrix whose index (i,j) represents the p-value significance
        of whether gene i and gene j are proximal gene pairs in that cell.
            Arguments: 
                - distance_threshold: (Integer) distance threshold at which to consider 2 genes proximal.
                - min_genecount: (integer) Minimum number of transcripts in each cell.
                - pval_matrix_name: (String) (Optional) if provided saves pvalue matrix using pickle at the input path.
                - gene_count_name: (String) (Optional) if provided saves gene expression count matrix using pickle at the input path.
        '''
        self.distance_threshold = distance_threshold
        cell_ids = np.unique(self.df.uID)
        num_cells = len(cell_ids)
        num_genes = len(self.geneList)
        print(f"Initialised PP now on {self.threads} threads")
        print("Number of cells: ", num_cells, ", Number of Genes: ", num_genes)
        start = timeit.default_timer()
        valid_cell_data = []
        num_transcripts = 0
        celllabel_encoder = LabelEncoder()
        cell_labels = celllabel_encoder.fit_transform(self.df['uID'])
        self.df['uID_encoded'] = cell_labels
        genes = self.df.index
        ids, counts = np.unique(self.df.uID_encoded.values, return_counts=True)
        for n, cell_id in enumerate(ids):
            num_transcripts += counts[n]
            if counts[n] > min_genecount:
                valid_cell_data.append([n, cell_id, genes])
            else:
                print(f"min genecount less than {min_genecount} for cell id {cell_ids[n]} across all Z, Skipping ...")
        df_copy = self.df.copy()
        del self.df
        position_matrix = df_copy[['absX', 'absY', 'absZ', 'uID_encoded']].to_numpy().copy(order='C')
        share_pp = mp.RawArray('d', len(position_matrix)*4)
        share_pp_np = np.frombuffer(share_pp).reshape(position_matrix.shape)
        np.copyto(share_pp_np, position_matrix)
        check = timeit.default_timer()
        pool = mp.Pool(processes = self.threads)
        with mp.Pool(processes=self.threads, initializer=self._initializer_func_pp, initargs=(position_matrix, position_matrix.shape), maxtasksperchild = 1) as pool:
            print(f"Running PP-3D now on {self.threads} threads for, {len(valid_cell_data)} cells, {num_transcripts} transcripts")
            results = pool.map(self._calculate_ProximalPairs3D, valid_cell_data)
        self.all_pval = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = self.precision_mode)
        self.all_gene_count = np.zeros((num_cells, len(self.geneList)), dtype = self.precision_mode) 
        for cell_i_result in results:
            self.all_pval[cell_i_result[0]] = cell_i_result[1]
            self.all_gene_count[cell_i_result[0]] = cell_i_result[2]
        print(f"Time to complete PP : ", timeit.default_timer() - check)
        if Path(self.filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.filename)
            adata.obsm['genecount'] = self.all_gene_count
            adata.uns['pp_test_pvalues'] = self.all_pval
            adata.write(self.filename, compression="gzip")
        if pval_matrix_name:
            self._save_pval_matrix(pval_matrix_name)
        if gene_count_name:
            self._save_gene_count(gene_count_name)
        self.df = df_copy
        del df_copy
        print(f"Cell-wise Proximal Pairs Time : {round(timeit.default_timer() - start, 2)} seconds")
    
    def _spatial_category(self, distance_thresold_nucleus = 2.5, distance_threshold_cyto_nuclear = 2.5, distance_threshold_cyto_peri = 4):
        nuc_distance = self.curr_cell_df.distNucleus.values
        cyt_distance = self.curr_cell_df.distPeriphery.values
        inNucleus = self.curr_cell_df.inNucleus.values
        inner_nuclear = np.zeros(nuc_distance.shape)
        peri_nuclear = np.zeros(nuc_distance.shape)
        cytosolic = np.zeros(nuc_distance.shape)
        peri_membrane = np.zeros(nuc_distance.shape)

        inner_nuclear[(inNucleus == 1) & (nuc_distance > distance_thresold_nucleus )] = 1
        peri_nuclear[(inNucleus == 1) & (nuc_distance <= distance_thresold_nucleus )] = 1
        peri_nuclear[(inNucleus == 0) & (nuc_distance <= distance_threshold_cyto_nuclear )] = 1
        cytosolic[(inNucleus == 0) & (nuc_distance > distance_threshold_cyto_nuclear) & (cyt_distance > distance_threshold_cyto_peri)] = 1
        peri_membrane[(inNucleus == 0) & (cyt_distance <= distance_threshold_cyto_peri)] = 1
        return [inner_nuclear, peri_nuclear, cytosolic, peri_membrane]
    
    def _category_counter(self, all_pairs, category_pairs):
        counter = {}
        spatial_cat_0,spatial_cat_1,spatial_cat_2,spatial_cat_3 = {},{},{},{}
        pairs_isinnernuc, pairs_isperinuc, pairs_iscyto, pairs_iscellperi = category_pairs[0], category_pairs[1], category_pairs[2], category_pairs[3]
        for i,elem in enumerate(all_pairs):
            counter[elem] = counter.get(elem, 0) + 1
            spatial_cat_0[elem] = spatial_cat_0.get(elem, 0) + pairs_isinnernuc[i]
            spatial_cat_1[elem] = spatial_cat_1.get(elem, 0) + pairs_isperinuc[i]
            spatial_cat_2[elem] = spatial_cat_2.get(elem, 0) + pairs_iscyto[i]
            spatial_cat_3[elem] = spatial_cat_3.get(elem, 0) + pairs_iscellperi[i]
        return counter, spatial_cat_0, spatial_cat_1, spatial_cat_2, spatial_cat_3
    
    def annotate_ProximalPairs(self, distance_threshold, distance_thresold_nucleus = 2.5, distance_threshold_cyto_nuclear = 2.5, distance_threshold_cyto_peri = 4):
        self.distance_threshold = distance_threshold
        cell_ids = np.unique(self.df.uID)
        num_cells = len(cell_ids)
        self.inner_nuc = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = self.precision_mode)
        self.peri_nuc = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = self.precision_mode)
        self.cytosolic = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = self.precision_mode)
        self.perimem = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = self.precision_mode)
        for n, cell_id in enumerate(cell_ids):
            self.curr_cell_df = self.df[self.df.uID == cell_id].copy()
            spatial_cat_all = self._spatial_category(distance_thresold_nucleus, distance_threshold_cyto_nuclear, distance_threshold_cyto_peri)
            self.curr_cell_df.loc[:,'InnerNuc'] = spatial_cat_all[0]
            self.curr_cell_df.loc[:,'PeriNuc'] = spatial_cat_all[1]
            self.curr_cell_df.loc[:,'Cytosolic'] = spatial_cat_all[2]
            self.curr_cell_df.loc[:,'PeriMem'] = spatial_cat_all[3]
            point_tree = cKDTree(self.curr_cell_df)
            queried_pairs = point_tree.query_pairs(self.distance_threshold)

            pairs = [(self.curr_cell_df.index[i], self.curr_cell_df.index[j]) for (i,j) in queried_pairs]
            pairs_isinnernuc = [(self.curr_cell_df.iloc[i,:]['InnerNuc'] + self.curr_cell_df.iloc[j,:]['InnerNuc'])/2 for (i,j) in queried_pairs]  #update later
            pairs_isperinuc = [(self.curr_cell_df.iloc[i,:]['PeriNuc'] + self.curr_cell_df.iloc[j,:]['PeriNuc'])/2 for (i,j) in  queried_pairs] 
            pairs_iscyto = [(self.curr_cell_df.iloc[i,:]['Cytosolic'] + self.curr_cell_df.iloc[j,:]['Cytosolic'])/2 for (i,j) in  queried_pairs] 
            pairs_iscellperi = [(self.curr_cell_df.iloc[i,:]['PeriMem'] + self.curr_cell_df.iloc[j,:]['PeriMem'])/2 for (i,j) in queried_pairs] 
            pairs, spatial_cat_0, spatial_cat_1, spatial_cat_2, spatial_cat_3  = self._category_counter(pairs, [pairs_isinnernuc, pairs_isperinuc, pairs_iscyto, pairs_iscellperi])
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
                col_na_genes = list(set(self.geneList).difference(obs_df.columns))
                row_na_genes = list(set(self.geneList).difference(obs_df.index))
                obs_df.loc[ : , col_na_genes] = 0
                obs_cat_0.loc[ : , col_na_genes] = 0
                obs_cat_1.loc[ : , col_na_genes] = 0
                obs_cat_2.loc[ : , col_na_genes] = 0
                obs_cat_3.loc[ : , col_na_genes] = 0

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
                obs_cat_0, obs_cat_1, obs_cat_2, obs_cat_3  = obs_df.copy(), obs_df.copy(), obs_df.copy(),obs_df.copy()
            temp = obs_cat_0.values
            arr_cat_0 = temp + temp.T - temp*np.identity(len(temp))
            temp = obs_cat_1.values
            arr_cat_1 = temp + temp.T - temp*np.identity(len(temp))
            temp = obs_cat_2.values
            arr_cat_2 = temp + temp.T - temp*np.identity(len(temp))
            temp = obs_cat_3.values
            arr_cat_3 = temp + temp.T - temp*np.identity(len(temp))
            self.inner_nuc[n] = arr_cat_0
            self.peri_nuc[n] = arr_cat_1
            self.cytosolic[n] = arr_cat_2
            self.perimem[n] = arr_cat_3
        return self.inner_nuc, self.peri_nuc, self.cytosolic, self.perimem

    def _calculate_ProximalPairs3D_slice(self, args):
        cell_num, cell_id, genes = args[0], args[1], args[2]
        pp_model = ProximalPairs3D(geneList = self.geneList, df_loc = cell_id, distance_threshold = self.distance_threshold, 
                                   min_genecount = self.min_genecount, cell_id = cell_id, genes = genes)
        return cell_num, pp_model.p_vals, pp_model.genecount.reshape(len(self.geneList))
    
    def run_ProximalPairs3D_slice(self, distance_threshold, min_genecount, pval_matrix_name = None, gene_count_name = None):
        '''
        Calculates Proximal pairs for each gene pair for each input cell. Runs the ProximalPairs() class for each cell
        and generates a 2d (num_genes, num_genes) matrix whose index (i,j) represents the p-value significance
        of whether gene i and gene j are proximal gene pairs in that cell.
            Arguments: 
                - distance_threshold: (Integer) distance threshold at which to consider 2 genes proximal.
                - min_genecount: (integer) Minimum number of transcripts in each cell.
                - pval_matrix_name: (String) (Optional) if provided saves pvalue matrix using pickle at the input path.
                - gene_count_name: (String) (Optional) if provided saves gene expression count matrix using pickle at the input path.
        '''
        self.distance_threshold = distance_threshold
        self.min_genecount = min_genecount
        cell_ids = np.unique(self.df.uID)
        num_cells = len(cell_ids)
        num_genes = len(self.geneList)
        print(f"Initialised PP-3D now on {self.threads} threads")
        print("Number of cells: ", num_cells, ", Number of Genes: ", num_genes)
        start = timeit.default_timer()
        valid_cell_data = []
        num_transcripts = 0
        celllabel_encoder = LabelEncoder()
        cell_labels = celllabel_encoder.fit_transform(self.df['uID'])
        self.df['uID_encoded'] = cell_labels
        genes = self.df.index
        ids, counts = np.unique(self.df.uID_encoded.values, return_counts=True)
        for n, cell_id in enumerate(ids):
            num_transcripts += counts[n]
            if counts[n] > min_genecount:
                valid_cell_data.append([n, cell_id, genes])
            else:
                print(f"min genecount less than {min_genecount} for cell id {cell_id}, Skipping ...")
        check = timeit.default_timer()
        global df
        df = self.df.copy()
        del self.df
        position_matrix_3d = df[['absX', 'absY', 'absZ', 'uID_encoded']].to_numpy().copy(order='C')
        share_pp3d = mp.RawArray('d', len(position_matrix_3d)*4)
        share_pp3d_np = np.frombuffer(share_pp3d).reshape(position_matrix_3d.shape)
        np.copyto(share_pp3d_np, position_matrix_3d)
        print(f"Running PP-3D slice now on {self.threads} threads for, {len(valid_cell_data)} cells, {num_transcripts} transcripts")
        with mp.Pool(processes=self.threads, initializer=self._initializer_func_pp, initargs=(position_matrix_3d, position_matrix_3d.shape), maxtasksperchild = 1) as pool:
            results = pool.map(self._calculate_ProximalPairs3D_slice, valid_cell_data)
        self.all_pval = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = self.precision_mode)
        self.all_gene_count = np.zeros((num_cells, len(self.geneList)), dtype = self.precision_mode) 
        for cell_i_result in results:
            self.all_pval[cell_i_result[0]] = cell_i_result[1]
            self.all_gene_count[cell_i_result[0]] = cell_i_result[2]
        print(f"Time to complete PP : ", timeit.default_timer() - check)
        if Path(self.filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.filename)
            adata.obsm['genecount'] = self.all_gene_count
            adata.uns['pp_test_pvalues'] = self.all_pval
            adata.write(self.filename, compression="gzip")
        if pval_matrix_name:
            self._save_pval_matrix(pval_matrix_name)
        if gene_count_name:
            self._save_gene_count(gene_count_name)
        self.df = df
        del df
        print(f"Cell-wise Proximal Pairs Time : {round(timeit.default_timer() - start, 2)} seconds")
    
    def run_spatial_modulation(self, inter_cell_distance, cell_locations = None, spatial_modulation_name = None, alpha = 0.01, randomize = False):
        '''
        Probabilistic graphical model to detect spatially modulated gene pairs.
        Requires `run_ProximalPairs()` be run first to generate the p-value matrix for all cells. 
        Arguments: 
            - inter_cell_distance: (Float) Maximum distance between cells at which they are considered proximal.
            - cell_locations: (String) (Optional) Path to file contains locations for each cell. Should be in sorted order. If not provided, cell locations are expected to be provided in `adata.uns['cell_locations']` in the AnnData file specified during initialization.
            - spatial_modulation_name: (String) (Optional) Path and name of the output excel file.
            - alpha: (Float) (Optional) pvalue signifcance threshold (>alpha_cellwise are converted to 1). Default = 0.01.
            - randomize: (Boolean) (Optional) Shuffle cell locations. Default = False".
        '''
        if cell_locations != None:
            if Path(cell_locations).suffix.lower() == ".csv":
                self.cell_locations = pd.read_csv(cell_locations)
        elif Path(self.filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.filename)
            self.cell_locations = adata.uns['cell_locations']
            self.all_pval = adata.uns['pp_test_pvalues']
            self.all_gene_count = adata.obsm['genecount']
        else:
            raise("Cell Locations file format incorrect")  
        print(self.cell_locations) 
        binary_pp_pval = np.zeros(self.all_pval.shape)
        binary_pp_pval[self.all_pval < alpha] = 1
        if randomize:
            self.cell_locations.uID = np.random.permutation(self.cell_locations.uID.values)
        SpatialModulation(binary_pp_pval, self.cell_locations, inter_cell_distance, file_name = spatial_modulation_name, geneList = self.geneList, threads = self.threads, adata_filename = self.filename)

    def run_differentialcolocalization(self, cell_type, cell_labels = None, file_location = None, cell_type_2 = None, mode = "1va", alpha = 0.01, folder_name = "differential_colocalization"):
        '''
        Calculates differentially colocaliztion between cell types across all the cells.
        There are 3 different modes in which it can be run - 
            - "1va" : Compares colocalization for genes in the input cell type vs all other cell types.
            - "1v1" : Compares colocalization for genes in the input cell type 1 vs input cell type 2.
            - "ava" : Compares colocalization for genes for all cell types vs all other cell types.
        Requires `run_ProximalPairs()` be run first to generate the p-value matrix for all cells. 
        Arguments: 
            - cell_type: (String) Cell type to calculate differential colocalization for. Is ignored if mode == "ava".
            - cell_labels: (String) (Optional) Path to file contains cell type for each cell. If not provided, cell labels are expected to be provided in `adata.obs` in the AnnData file specified during initialization. 
            - file_location: (String) (Optional) Directory in which to store output files. if mode == "ava", creates a new directory in this path to store all results.
            - cell_type_2: (String) (Optional) Cell type 2 to calculate differential colocalization for. Required if mode == "1v1".
            - mode: (String) Either "1va" (One cell type vs All cell types), "1v1" (One cell type vs one cell type) or "ava" (All cell types vs All cell types).
            - alpha: (Float) (Optional) pvalue signifcance threshold (>alpha_cellwise are converted to 1). Default = 0.01.
            - folder_name: (String) (Optional) if mode == "ava", folder name inside specified path in which to store results for each cell type. Default = "differential_colocalization".
        '''
        celllabel_encoder = LabelEncoder()
        cell_ids = celllabel_encoder.fit_transform(self.df['uID'])
        cell_ids, _ = np.unique(cell_ids, return_counts=True)
        if cell_labels != None:
            if Path(cell_labels).suffix.lower() == ".csv":
                self.cell_labels = pd.read_csv(cell_labels)
        elif Path(self.filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.filename)
            self.cell_labels = adata.obs
            self.cell_labels['uID'] = self.cell_labels.index.values
            self.all_pval = adata.uns['pp_test_pvalues']
            self.all_gene_count = adata.obsm['genecount']
            if 'differential_colocalization' not in adata.uns:
                adata.uns['differential_colocalization'] = {}
            adata.write(self.filename, compression="gzip")
        else:
            raise("Cell Label file format incorrect")        
        start = timeit.default_timer()
        print(f"Running Differential Colocalization now on {self.threads} threads")
        if mode == "1va":
            selected_cell_ids = self.cell_labels[self.cell_labels.cell_type == cell_type].uID.values
            rest_cell_ids = list(set(self.cell_labels.uID.values).difference(set(selected_cell_ids)))
            selected_cell_indices = [np.where(cell_ids == celllabel_encoder.transform([x])[0])[0][0] for x in selected_cell_ids]
            rest_cell_indices = [np.where(cell_ids == celllabel_encoder.transform([x])[0])[0][0] for x in rest_cell_ids]
            DifferentialColocalization(self.all_pval, self.all_gene_count, self.geneList, selected_cell_indices, rest_cell_indices, file_location=file_location, cell_type=cell_type, threads = self.threads, alpha = alpha, filename = self.filename)
            print(f"Differential Colocalization for cell type {cell_type} Time : {round(timeit.default_timer() - start, 2)} seconds")
        elif mode == "1v1":
            selected_cell_ids = self.cell_labels[self.cell_labels.cell_type == cell_type].uID.values
            rest_cell_ids = self.cell_labels[self.cell_labels.cell_type == cell_type_2].uID.values
            selected_cell_indices = [np.where(cell_ids == celllabel_encoder.transform([x])[0])[0][0] for x in selected_cell_ids]
            rest_cell_indices = [np.where(cell_ids == celllabel_encoder.transform([x])[0])[0][0] for x in rest_cell_ids]
            DifferentialColocalization(self.all_pval, self.all_gene_count, self.geneList, selected_cell_indices, rest_cell_indices, file_location=file_location, cell_type=f"{cell_type}_vs_{cell_type_2}", threads = self.threads, alpha = alpha, filename = self.filename)
            print(f"Differential Colocalization for cell type {cell_type} Time : {round(timeit.default_timer() - start, 2)} seconds")
        else:
            if file_location != None:
                file_location = Path(file_location) / folder_name
                subprocess.run(["mkdir", file_location])
            for selected_cell_type in self.cell_labels.cell_type.unique():
                file_location_ct = None
                if file_location != None:
                    file_location_ct = Path(file_location) / str(selected_cell_type)
                    subprocess.run(["mkdir", file_location_ct])
                selected_cell_ids = self.cell_labels[self.cell_labels.cell_type == selected_cell_type].uID.values
                rest_cell_ids = list(set(self.cell_labels.uID.values).difference(set(selected_cell_ids)))
                selected_cell_indices = [np.where(cell_ids == celllabel_encoder.transform([x])[0])[0][0] for x in selected_cell_ids]
                rest_cell_indices = [np.where(cell_ids == celllabel_encoder.transform([x])[0])[0][0] for x in rest_cell_ids]
                DifferentialColocalization(self.all_pval, self.all_gene_count, self.geneList, selected_cell_indices, rest_cell_indices, file_location=file_location_ct, cell_type=selected_cell_type, threads = self.threads, alpha = alpha, filename = self.filename)
            print(f"All2All Differential Colocalization Time : {round(timeit.default_timer() - start, 2)} seconds")

    def _save_globcolocal_csv(self, filename):
        self.global_coloc_df.to_csv(filename)
    
    def _save_expcolocal_csv(self, filename):
        self.expected_coloc_df.to_csv(filename)

    def _save_unstacked_pvals(self, filename, alpha_cellwise, min_transcript):
        present_cells = self._num_present_cells(min_transcript)
        obs = self._binarize_adj_matrix(alpha_cellwise).sum(axis=0)
        unstacked_global_pvals = self._unstack_df_both(obs, present_cells, alpha_cellwise)
        if unstacked_global_pvals.shape[0] < 1048576 and unstacked_global_pvals.shape[1] < 16384:
            unstacked_global_pvals.to_excel(filename) 
        else:
            unstacked_global_pvals.to_csv(filename[:-5]+".csv") 
    
    def _save_unstacked_pvals_adata(self, alpha_cellwise, min_transcript):
        present_cells = self._num_present_cells(min_transcript)
        obs = self._binarize_adj_matrix(alpha_cellwise).sum(axis=0)
        unstacked_global_pvals = self._unstack_df_both(obs, present_cells, alpha_cellwise)
        return unstacked_global_pvals

    def _binarize_adj_matrix(self, alpha):
        edge_all  = np.zeros(self.all_pval.shape)
        edge_all[self.all_pval < alpha] = 1
        return edge_all

    def _num_present_cells(self, min_transcript):
        indicator = np.zeros(self.all_gene_count.shape)
        indicator[self.all_gene_count > min_transcript] = 1
        ind_n1n2 = np.matmul(indicator.reshape(indicator.shape[0],indicator.shape[1],1),indicator.reshape(indicator.shape[0], 1,indicator.shape[1]))
        present_cells = ind_n1n2.sum(axis=0)
        return present_cells

    def _unstack_df_both(self, obs, present_cells, alpha, second_df_name = 'Expected coloc'):
        num_genes = self.expected_coloc_df.values.shape[0]
        agg_coloc = self.expected_coloc_df.values[np.triu_indices(num_genes)]
        cond_agg_coloc = self.global_coloc_df.values[:,:][np.triu_indices(num_genes)]
        obs = obs[np.triu_indices(num_genes)]
        present_cells = present_cells[np.triu_indices(num_genes)]
        gene_id1 = [self.geneList[i] for i in np.triu_indices(num_genes)[0]]
        gene_id2 = [self.geneList[i] for i in np.triu_indices(num_genes)[1]]
        data = {'gene_id1': gene_id1, 'gene_id2': gene_id2, 'p_val_cond': cond_agg_coloc, second_df_name: agg_coloc, f'Coloc. cells(Threshold<{str(alpha)})': obs, 'Present cells' : present_cells}
        pairwise_p_val_df = pd.DataFrame(data)
        pairwise_p_val_df['frac_cells'] = pairwise_p_val_df[f'Coloc. cells(Threshold<{str(alpha)})']/pairwise_p_val_df['Present cells']
        pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'].astype(str) + ', ' + pairwise_p_val_df['gene_id2'].astype(str)
        pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
        pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['p_val_cond'])
        return pairwise_p_val_df
    
    def run_GlobalColocalization(self, alpha_cellwise = 0.01, min_transcript = 0, high_precision = False, glob_coloc_name = None, exp_coloc_name = None, unstacked_pvals_name = None):
        '''
        Calculates global colocalization significance for each gene pair across all the cells. 
        Requires `run_ProximalPairs()` be run first to generate the p-value matrix for all cells. 
        Generates 3 dataframes, all  of size (num_genes, num_genes)
        Arguments: 
            - alpha_cellwise: (Float) pvalue signifcance threshold (>alpha_cellwise are converted to 1). Default = 0.01.
            - min_transcript: (Float) Gene expression lower threshold. Default = 0.
            - high_precision: (Boolean) High precision pvalue. Expect longer computer. Default = False.
            - glob_coloc_name: (String) (Optional) if provided saves global colocalization matrix as a csv at the input path.
            - exp_coloc_name: (String) (Optional) if provided saves expected colocalization matrix as a csv at the input path.
            - unstacked_pvals_name: (String) (Optional) if provided saves interpretable global colocalization matrix as a csv at the input path.
        '''
        print(f"Running Global Colocalization now on {self.threads} threads")
        if Path(self.filename).suffix.lower() == ".h5ad":
            adata = ad.read_h5ad(self.filename)
            self.all_pval = adata.uns['pp_test_pvalues']
            self.all_gene_count = adata.obsm['genecount']
        print(self.all_pval.shape, self.all_gene_count.shape)
        start = timeit.default_timer()
        global_coloc_model = ConditionalGlobalColocalization(all_pvals = self.all_pval, transcript_count = self.all_gene_count, alpha_cellwise = alpha_cellwise, min_transcript = min_transcript, threads = self.threads, high_precision = high_precision, precision_mode = self.precision_mode)
        global_coloc, expected_coloc = global_coloc_model.global_colocalization()
        self.global_coloc_df = pd.DataFrame(global_coloc, index = self.geneList, columns = self.geneList)
        self.expected_coloc_df = pd.DataFrame(expected_coloc, index = self.geneList, columns = self.geneList)
        if Path(self.filename).suffix.lower() == ".h5ad":
            adata.uns['cpb_global_colocalization'] = self.global_coloc_df
            adata.uns['cpb_expected_colocalization'] = self.expected_coloc_df
            adata.uns['cpb_unstacked_colocalization'] = self._save_unstacked_pvals_adata(alpha_cellwise, min_transcript)
            adata.write(self.filename, compression="gzip")
        if glob_coloc_name:
            self._save_globcolocal_csv(glob_coloc_name)
        if exp_coloc_name:
            self._save_expcolocal_csv(exp_coloc_name)
        if unstacked_pvals_name:
            self._save_unstacked_pvals(unstacked_pvals_name, alpha_cellwise, min_transcript)
        print(f"Cell-wise Global Colocalization Time : {round(timeit.default_timer() - start, 2)} seconds")