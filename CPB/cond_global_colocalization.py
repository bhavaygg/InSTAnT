import matplotlib.pyplot as plt
import time
import numpy as np
import pandas as pd
import pickle
import os
# setting path

from poisson_binomial import PoissonBinomial
from poibin import PoiBin

class ConditionalGlobalColocalization():
    def __init__(self, all_cell_p_val, transcript_count, alpha_cellwise = 0.01, min_transcript = 0, show_det_pairs = 0):
        self.all_cell_p_val = all_cell_p_val
        self.alpha, self.show_det_pairs = alpha_cellwise, show_det_pairs
        self.transcript_count = transcript_count   #num_cells x num_genes
        self.min_transcript = min_transcript
        self.binarize_adj_matrix()
        self.global_weight_pair()
        self.prob_cells_coloc()

    def binarize_adj_matrix(self):
        edge_all  = np.zeros(self.all_cell_p_val.shape)
        edge_all[self.all_cell_p_val < self.alpha] = 1
        self.edge = edge_all
        self.obs = edge_all.sum(axis=0)
        self.num_cells = edge_all.shape[0]
        self.num_genes = edge_all.shape[1]
        print('num cells: %d, num_genes: %d' %(self.num_cells, self.num_genes))


    def global_weight_pair(self):
        global_genecount = self.edge.sum(axis = (0,1)).reshape([-1,1])
        weight_pair = np.matmul(global_genecount, np.transpose(global_genecount))
        self.gene_pair_weight = weight_pair

    def compute_pmf(self,i,j):
        p_ij = []
        for cell_id in range(self.num_cells):
            if self.transcript_count[cell_id, i] > self.min_transcript and self.transcript_count[cell_id, j] > self.min_transcript:
                p_ij.append(self.prob_pairwise_all_cells[cell_id,i,j])
#         pb = PoissonBinomial(p_ij)
        pb = PoiBin(p_ij)
        # print('test', self.obs[i,j].astype(int), len(p_ij), i,j)
        pmf_curr = pb.pmf(self.obs[i,j].astype(int))

        # if pval_curr < 1e-13:  #Poibin has numerical error in low p val region
        #     pb_highres = PoissonBinomial(p_ij)
        #     pval_curr = pb_highres.pval(self.obs[i,j].astype(int))
        
        return pmf_curr

    def find_likelihood(self):
        likelihood = 1e-128*np.ones((self.num_genes, self.num_genes))
        start_time = time.time()
        for i in range(self.num_genes):
            print('gene',i, time.time() - start_time)
            for j in range(i,self.num_genes):
                likelihood[i,j] = -np.log(self.compute_pmf(i,j)+1e-128)
                likelihood[j,i] = likelihood[i,j]
#                 print('gene', i, j, time.time() - start_time)
        return likelihood
        
    def compute_pval(self, i, j):
        p_ij = []
        for cell_id in range(self.num_cells):
            if self.transcript_count[cell_id, i] > self.min_transcript and self.transcript_count[cell_id, j] > self.min_transcript:
                p_ij.append(self.prob_pairwise_all_cells[cell_id,i,j])
#         pb = PoissonBinomial(p_ij)
        # print(p_ij)
        pb = PoiBin(p_ij)
        # print('test', self.obs[i,j].astype(int), len(p_ij), i,j)clear
        
        pval_curr = pb.pval(self.obs[i,j].astype(int))

        if pval_curr < 1e-13:  #Poibin has numerical error in low p val region
            pb_highres = PoissonBinomial(p_ij)
            pval_curr = pb_highres.pval(self.obs[i,j].astype(int))
        
        return pval_curr, np.sum(p_ij)
    
    def global_colocalization(self):
        coloc_matrix = np.ones((self.num_genes, self.num_genes))
        expected_coloc = np.zeros((self.num_genes, self.num_genes))
        start_time = time.time()
        for i in range(self.num_genes):
            print('gene',i, time.time() - start_time)
            for j in range(i,self.num_genes):
                coloc_matrix[i,j], expected_coloc[i,j] = self.compute_pval(i,j)
                coloc_matrix[j,i], expected_coloc[j,i] = coloc_matrix[i,j], expected_coloc[i,j]
#                 print('gene', i, j, time.time() - start_time)
        return coloc_matrix, expected_coloc
    
    def prob_cells_coloc(self):
        prob_pairwise_all_cells = np.ones((self.num_cells, self.num_genes, self.num_genes))
        
        for i in range(self.num_cells):
            curr_edge = self.edge[i].copy()
            curr_edge = curr_edge[np.triu_indices(curr_edge.shape[0])]
            
            #zero count genes weight
            curr_weight_pair = self.gene_pair_weight.copy() 
            curr_weight_pair[self.transcript_count[i] <= self.min_transcript] = 0
            temp = curr_weight_pair[np.triu_indices(curr_weight_pair.shape[0])]
            # assert temp.sum() != 0
            curr_weight_pair = curr_weight_pair/(temp.sum() + 1e-64)   #adding small number in denominator to avoid inf
            
            prob_pairwise_all_cells[i] = 1 - np.power((1 - curr_weight_pair), curr_edge.sum())

        self.prob_pairwise_all_cells = prob_pairwise_all_cells
    
                
        
        
