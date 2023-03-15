from __future__ import division
import os
import numpy as np
import sys
import pandas as pd
sys.path.append('/Users/anu/Desktop/InSTAnT/utils')
sys.path.append('/Users/anu/Desktop/InSTAnT/AggregationModel')
# import utils
# import AggregationModel
from CPB.cond_global_colocalization import ConditionalGlobalColocalization

from utils.helper import pkl_load,  unstack_df_both, agg_cmdline_args_parser, create_dir, num_present_cells, binarize_adj_matrix

parser = agg_cmdline_args_parser()
args = parser.parse_args()

##Hyperparams
alpha = args.alpha
min_transcript_cpb = args.min_transcript  #curr only for min transcript 0
d = int(args.dist)
rep = args.rep



data_dir = '/Users/anurendrakumar/Desktop/1Nature_paper/Results_final/U2OS_FILTERED/' + rep + '/'
f_all_p_val = data_dir + 'all_cells_p_val/pvalues_' + str(d) + '.pkl'
f_genecount = data_dir + 'all_cells_p_val/gene_count.pkl'
f_geneList = '/Users/anurendrakumar/Desktop/1Nature_paper/Results_final/U2OS_FILTERED/geneList.pkl'


all_p_vals, genecount, geneList = pkl_load(f_all_p_val), pkl_load(f_genecount ), pkl_load(f_geneList)

assert  genecount.shape[0] == all_p_vals.shape[0]  #num cells should be same

##Model
def run_aggregate_model(isConditional):
    expected_coloc_df = []
    print('running conditional aggregate model')
    agg_model = ConditionalGlobalColocalization(all_p_vals, genecount, alpha_cellwise = alpha, min_transcript = min_transcript_cpb)

    agg_coloc, expected_coloc = agg_model.global_colocalization()
    agg_coloc_df = pd.DataFrame(agg_coloc, index = geneList, columns = geneList)
    expected_coloc_df = pd.DataFrame(expected_coloc, index = geneList, columns = geneList)
    return agg_coloc_df, expected_coloc_df

save_dir = data_dir + 'agg/'  + str(alpha) + '/' + str(d)  + '/'
create_dir(save_dir)

if not os.path.exists(save_dir + 'cond_agg_p_val.csv'):
    cond_agg_coloc_df, expected_coloc_df = run_aggregate_model(isConditional = True)
    cond_agg_coloc_df.to_csv(save_dir + 'cond_agg_p_val.csv') 
    expected_coloc_df.to_csv(save_dir + 'expected_coloc.csv') 
else:
    cond_agg_coloc_df = pd.read_csv(save_dir + 'cond_agg_p_val.csv', index_col=0)
    expected_coloc_df = pd.read_csv(save_dir + 'expected_coloc.csv', index_col = 0)


present_cells = num_present_cells(genecount, min_transcript_cpb)
obs = binarize_adj_matrix(all_p_vals, alpha).sum(axis=0)
p_vals_agg_unstacked = unstack_df_both(cond_agg_coloc_df.values[:,:], expected_coloc_df.values, obs, present_cells, geneList, alpha)
# p_vals_agg_unstacked.to_csv(save_dir + 'unstacked_agg_p_vals.csv')
p_vals_agg_unstacked.to_excel(save_dir + 'unstacked_agg_p_vals.xlsx') 
