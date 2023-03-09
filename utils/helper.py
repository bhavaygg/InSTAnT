from __future__ import division
import pickle
import numpy as np
import os
import argparse
import pandas as pd

def pkl_load(f):
    with open(f,'rb') as f:
        data = pickle.load(f)
    return np.asarray(data)

def pkl_save(data,f):
    with open(f, 'wb') as f:
        pickle.dump(data, f)
    
def print_and_log(msg, log_file, write_mode="a"):
    """
    print `msg` (string) on stdout and also append ('a') or write ('w') (default 'a') it to `log_file`
    """
    print(msg)
    with open(log_file, write_mode) as f:
        f.write(msg + "\n")

def agg_cmdline_args_parser():
    """
    Commandline argument parser
    """
    parser = argparse.ArgumentParser("Cmdline Arguments")
    parser.add_argument("-alpha", "--alpha", type=float, default=0.01)
    parser.add_argument("-min_transcript", "--min_transcript", type=int, default=0)
    # parser.add_argument("-filter", "--isFiltered", type=bool, default=False)
    parser.add_argument("-ct", "--cell_type_name", type=str, default= 'exc')
    parser.add_argument("-d", "--dist", type=float, default= 4)
    parser.add_argument("-rep", "--rep", type=str, default= 'rep3')
    parser.add_argument("-inter_dist", "--inter_dist", type=int, default= 500)

    return parser

def num_present_cells(genecount, min_transcript):
    indicator = np.zeros(genecount.shape)
    indicator[genecount>min_transcript] = 1
    ind_n1n2 = np.matmul(indicator.reshape(indicator.shape[0],indicator.shape[1],1),indicator.reshape(indicator.shape[0], 1,indicator.shape[1]))
    present_cells = ind_n1n2.sum(axis=0)
    return present_cells

def unstack_df_both(cond_agg_coloc, agg_coloc, obs, present_cells, geneList, alpha, second_df_name = 'Expected coloc'):
    num_genes = agg_coloc.shape[0]
    agg_coloc = agg_coloc[np.triu_indices(num_genes)]
    cond_agg_coloc = cond_agg_coloc[np.triu_indices(num_genes)]
    obs = obs[np.triu_indices(num_genes)]
    present_cells = present_cells[np.triu_indices(num_genes)]
    gene_id1 = [geneList[i] for i in np.triu_indices(num_genes)[0]]
    gene_id2 = [geneList[i] for i in np.triu_indices(num_genes)[1]]
    data = { 'gene_id1': gene_id1, 'gene_id2': gene_id2, 'p_val_cond': cond_agg_coloc, second_df_name: agg_coloc, 'Coloc. cells(Threshold<'+ str(alpha) +')':obs, 'Present cells' : present_cells}
    
    pairwise_p_val_df = pd.DataFrame(data)
    
    pairwise_p_val_df['frac_cells'] = pairwise_p_val_df['Coloc. cells(Threshold<' + str(alpha)  +')']/pairwise_p_val_df['Present cells']
    pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'] + ', ' + pairwise_p_val_df['gene_id2']
#     pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['frac_cells'], ascending=False)
    pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
    # pairwise_p_val_df = pairwise_p_val_df.sort_values(by=second_df_name)
    pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['p_val_cond'])
    return pairwise_p_val_df

def unstack_df(p_vals_all, obs, present_cells, geneList, alpha):
    num_genes = p_vals_all.shape[0]
    p_val_all = p_vals_all[np.triu_indices(num_genes)]
    obs = obs[np.triu_indices(num_genes)]
    present_cells = present_cells[np.triu_indices(num_genes)]
    gene_id1 = [geneList[i] for i in np.triu_indices(num_genes)[0]]
    gene_id2 = [geneList[i] for i in np.triu_indices(num_genes)[1]]
    data = { 'gene_id1': gene_id1, 'gene_id2': gene_id2, 'p_vals_all': p_val_all, 'Coloc. cells(Threshold<'+ str(alpha) +')':obs, 'Present cells' : present_cells}
    pairwise_p_val_df = pd.DataFrame(data)
    pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['p_vals_all'])
    pairwise_p_val_df['frac_cells'] = pairwise_p_val_df['Coloc. cells(Threshold<' + str(alpha)  +')']/pairwise_p_val_df['Present cells']
    pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'] + ', ' + pairwise_p_val_df['gene_id2']
#     pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['frac_cells'], ascending=False)
    pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
    return pairwise_p_val_df

def unstack_df_corr(corr_matrix,  geneList):
    print(corr_matrix.shape)
    num_genes = corr_matrix.shape[0]
    corr_matrix = corr_matrix[np.triu_indices(num_genes)]
    gene_id1 = [geneList[i] for i in np.triu_indices(num_genes)[0]]
    gene_id2 = [geneList[i] for i in np.triu_indices(num_genes)[1]]
    data = { 'gene_id1': gene_id1, 'gene_id2': gene_id2, 'corr': corr_matrix}
    pairwise_p_val_df = pd.DataFrame(data)
    pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['corr'], ascending=False)
    pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'] + ', ' + pairwise_p_val_df['gene_id2']
#     pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['frac_cells'], ascending=False)
    pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
    return pairwise_p_val_df


def unstack_df_simple(p_vals_all, geneList,column_name = 'p_vals_all'):
    num_genes = p_vals_all.shape[0]
    p_vals_all = p_vals_all[np.triu_indices(num_genes)]
    gene_id1 = [geneList[i] for i in np.triu_indices(num_genes)[0]]
    gene_id2 = [geneList[i] for i in np.triu_indices(num_genes)[1]]
    data = { 'gene_id1': gene_id1, 'gene_id2': gene_id2, 'llr': p_vals_all}
    pairwise_p_val_df = pd.DataFrame(data)
    pairwise_p_val_df = pairwise_p_val_df.sort_values(by=[column_name],ascending=False)
    pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'] + ', ' + pairwise_p_val_df['gene_id2']
    pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
    return pairwise_p_val_df

def create_dir(save_dir):
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

def binarize_adj_matrix(all_cell_p_val, alpha):
    edge_all  = np.zeros(all_cell_p_val.shape)
    edge_all[all_cell_p_val < alpha] = 1
    return edge_all

def load_cells_loc(cells_loc_path,cell_id_list):
    cells_loc = pd.read_csv(cells_loc_path, sep = '\t', names = ['uID','absX', 'absY'], header = None)
    cells_loc = cells_loc.set_index('uID')
    return cells_loc.loc[cell_id_list]

def load_cells_loc_brain(cells_loc_path,cell_id_list):
    cells_loc = pd.read_csv(cells_loc_path)
    cells_loc = cells_loc[['feature_ID','centroid_1','centroid_2']]
    cells_loc = cells_loc.set_index('feature_ID')
    return cells_loc.loc[cell_id_list]

def unstack_df_llr(llr, w_h1_all, p_h1_all, p_h0_all, geneList):
    num_genes = len(geneList)
    llr = llr[np.triu_indices(num_genes)]
    w_h1_all = w_h1_all[np.triu_indices(num_genes)]
    p_h1_all = p_h1_all[np.triu_indices(num_genes)]
    p_h0_all = p_h0_all[np.triu_indices(num_genes)]
    gene_id1 = [geneList[i] for i in np.triu_indices(num_genes)[0]]
    gene_id2 = [geneList[i] for i in np.triu_indices(num_genes)[1]]
    data = { 'gene_id1': gene_id1, 'gene_id2': gene_id2, 'llr': llr, 'w_h1': w_h1_all, 'p_g_h1':p_h1_all, 'p_g_h0' : p_h0_all}
    
    pairwise_p_val_df = pd.DataFrame(data)

    pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'] + ', ' + pairwise_p_val_df['gene_id2']
    pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
    pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['llr'],ascending=False)
    return pairwise_p_val_df


