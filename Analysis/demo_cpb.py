from __future__ import division
import os
import numpy as np
import sys
import pandas as pd
sys.path.append('/Users/anu/Desktop/InSTAnT/utils')
sys.path.append('/Users/anu/Desktop/InSTAnT/AggregationModel')
import utils
import AggregationModel
from AggregationModel.cond_global_colocalization import ConditionalGlobalColocalization
from AggregationModel.global_colocalization import GlobalColocalization
from utils.constants import  FilePathU2OSHanAllCells, FileFalsePathU2OSArea3Filtered, FilePathU2OSArea3Filtered,  FilePathU2OSHanAllCells, FilePathU2OSArea3FilteredNucleus, FilePathM1HanAllCells,FileFalsePathU2OSHanAllCells, FileFalsePathM1HanAllCells
from utils.helper import pkl_load,  unstack_df_both, agg_cmdline_args_parser, create_dir, num_present_cells, binarize_adj_matrix

parser = agg_cmdline_args_parser()
args = parser.parse_args()

##Hyperparams
alpha = args.alpha
min_transcript_cpb = args.min_transcript  #curr only for min transcript 0
d = int(args.dist)
rep = args.rep

##filepaths
# filepath = FilePathU2OSArea3FilteredNucleus(rep,d)#Nucleus
filepath = FilePathU2OSArea3Filtered(rep,d)
# filepath = FilePathM1HanAllCells
# filepath = FileFalsePathU2OSHanAllCells(d)
# filepath = FileFalsePathM1HanAllCells
# filepath = FileFalsePathU2OSAllCells('rep4')
f_pvalues,f_genecount,f_geneList = filepath.f_all_p_val, filepath.f_genecount, filepath.f_geneList
all_p_vals, genecount, geneList = pkl_load(f_pvalues), pkl_load(f_genecount), pkl_load(f_geneList)

assert  genecount.shape[0] == all_p_vals.shape[0]  #num cells should be same

##Model
def run_aggregate_model(isConditional):
    expected_coloc_df = []
    if isConditional:
        print('running conditional aggregate model')
        agg_model = ConditionalGlobalColocalization(all_p_vals, genecount, alpha_cellwise = alpha, min_transcript = min_transcript_cpb)
    else:
        print('running aggregate model')
        agg_model = GlobalColocalization(all_p_vals, genecount, alpha_cellwise = alpha, min_transcript = min_transcript_cpb)
    agg_coloc, expected_coloc = agg_model.global_colocalization()
    agg_coloc_df = pd.DataFrame(agg_coloc, index = geneList, columns = geneList)
    expected_coloc_df = pd.DataFrame(expected_coloc, index = geneList, columns = geneList)
    return agg_coloc_df, expected_coloc_df

save_dir = filepath.f_agg + str(alpha) + '/' + str(d)  + '/'

create_dir(save_dir)

if not os.path.exists(save_dir + 'cond_agg_p_val.csv'):
    cond_agg_coloc_df, expected_coloc_df = run_aggregate_model(isConditional = True)
    cond_agg_coloc_df.to_csv(save_dir + 'cond_agg_p_val.csv') 
    expected_coloc_df.to_csv(save_dir + 'expected_coloc.csv') 
else:
    cond_agg_coloc_df = pd.read_csv(save_dir + 'cond_agg_p_val.csv', index_col=0)
    expected_coloc_df = pd.read_csv(save_dir + 'expected_coloc.csv', index_col = 0)

# if not os.path.exists(save_dir + 'agg_p_val.csv'):
#     agg_coloc_df, _ = run_aggregate_model(isConditional = False)
#     agg_coloc_df.to_csv(save_dir + 'agg_p_val.csv')
# else:
#     agg_coloc_df = pd.read_csv(save_dir + 'agg_p_val.csv', index_col = 0)


present_cells = num_present_cells(genecount, min_transcript_cpb)
obs = binarize_adj_matrix(all_p_vals, alpha).sum(axis=0)
num_blank_genes=10
p_vals_agg_unstacked = unstack_df_both(cond_agg_coloc_df.values[num_blank_genes:,num_blank_genes:], expected_coloc_df.values[num_blank_genes:,num_blank_genes:], obs[num_blank_genes:,num_blank_genes:], present_cells[num_blank_genes:,num_blank_genes:], geneList[num_blank_genes:], alpha)
# p_vals_agg_unstacked.to_csv(save_dir + 'unstacked_agg_p_vals.csv')
p_vals_agg_unstacked.to_excel(save_dir + 'unstacked_agg_p_vals.xlsx') 
