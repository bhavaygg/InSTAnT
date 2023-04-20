import pandas as pd
import numpy as np
import pickle
import timeit
from scipy.spatial import cKDTree
from scipy.stats import binom_test
from collections import Counter
from scipy.sparse import coo_matrix
import multiprocessing as mp
from sklearn.neighbors import NearestNeighbors
from InSTAnT.poibin import PoiBin
from InSTAnT.poisson_binomial import PoissonBinomial
from numpy.random import default_rng


class ConditionalGlobalColocalization():
    '''
        Performs conditional global colocalization
        Requires output from ProximalPairs()
        Arguments: 
            - all_cell_pval: (Array) Gene-Gene pairwise pvalues for all cells calcuated using ProximalPairs().
            - transcript_count: (Array) Expression count of each gene across each cell.
            - alpha_cellwise: (Float) pvalue signifcance threshold (>alpha_cellwise are converted to 1).
            - min_transcript: (Float) Gene expression lower threshold.
            - show_det_pairs: (Not used)
            - high_precision: (Boolean) High precision pvalue. Expect longer computer.
            - threads: (Integer) Number of threads to use.
    '''
    def __init__(self, all_cell_pval, transcript_count, alpha_cellwise = 0.01, min_transcript = 0, show_det_pairs = 0, high_precision = False, threads = 1):
        self.all_cell_pval = all_cell_pval
        self.alpha, self.show_det_pairs = alpha_cellwise, show_det_pairs
        self.transcript_count = transcript_count   #num_cells x num_genes
        self.min_transcript = min_transcript
        self.high_precision = high_precision
        self.threads = threads
        self.binarize_adj_matrix()
        self.global_weight_pair()
        self.prob_cells_coloc()

    def binarize_adj_matrix(self):
        '''
        Convert all pvals > alpha to 1
        '''
        edge_all  = np.zeros(self.all_cell_pval.shape)
        edge_all[self.all_cell_pval < self.alpha] = 1
        self.edge = edge_all
        self.obs = edge_all.sum(axis=0) #pairwise colocalisation
        self.num_cells = edge_all.shape[0]
        self.num_genes = edge_all.shape[1]
        print('Number of cells: %d, Number of genes: %d' %(self.num_cells, self.num_genes))

    def global_weight_pair(self):
        global_genecount = self.edge.sum(axis = (0,1)).reshape([-1,1])
        weight_pair = np.matmul(global_genecount, np.transpose(global_genecount))
        self.gene_pair_weight = weight_pair

    def prob_cells_coloc(self):
        prob_pairwise_all_cells = np.ones((self.num_cells, self.num_genes, self.num_genes))
        for i in range(self.num_cells):
            curr_edge = self.edge[i].copy()
            curr_edge = curr_edge[np.triu_indices(curr_edge.shape[0])]
            #zero count genes weight
            curr_weight_pair = self.gene_pair_weight.copy() 
            curr_weight_pair[self.transcript_count[i] <= self.min_transcript] = 0 #removing lower count genes
            temp = curr_weight_pair[np.triu_indices(curr_weight_pair.shape[0])]
            # assert temp.sum() != 0
            curr_weight_pair = curr_weight_pair/(temp.sum() + 1e-64)   #adding small number in denominator to avoid inf
            prob_pairwise_all_cells[i] = 1 - np.power((1 - curr_weight_pair), curr_edge.sum())
        self.prob_pairwise_all_cells = prob_pairwise_all_cells
        
    def _compute_pval(self, args):
        i, j = args[0], args[1]
        p_ij = []
        for cell_id in range(self.num_cells):
            if self.transcript_count[cell_id, i] > self.min_transcript and self.transcript_count[cell_id, j] > self.min_transcript:
                p_ij.append(self.prob_pairwise_all_cells[cell_id,i,j])
        pb = PoiBin(p_ij)
        pval_curr = pb.pval(self.obs[i,j].astype(int))
        if self.high_precision and pval_curr < 1e-13:  #Poibin has numerical error in low p val region
            pb_highres = PoissonBinomial(p_ij)
            pval_curr = pb_highres.pval(self.obs[i,j].astype(int))
        return i, j, pval_curr, np.sum(p_ij)
    
    def global_colocalization(self):
        coloc_matrix = np.ones((self.num_genes, self.num_genes))
        expected_coloc = np.zeros((self.num_genes, self.num_genes))
        start = timeit.default_timer()
        pool = mp.Pool(self.threads)
        results = pool.map(self._compute_pval, [[i, j] for i in range(self.num_genes) for j in range(i, self.num_genes)])
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
    

class ProximalPairs():
    '''
        Calculates proximal pairs for all gene-gene pairs for a given cell.
        Arguments: 
            - geneList: (Array) List of genes.
            - df_loc: (DataFrame) Coordinates of the gene transcript.
            - distance_threshold: (Integer) distance threshold at which to consider 2 genes proximal.
            - mode: (Not used).
    '''
    def __init__(self, geneList, df_loc, distance_threshold, mode="normal"):
        self.geneList = geneList
        #self.orig_df = df_loc.copy()
        self.curr_cell_df = df_loc[['absX', 'absY']]
        self.distance_threshold = distance_threshold
        if mode=='normal':
            self.genecount, self.num_trial = self.num_trial_pairs()
            self.prob_null = self.prob_null_hypothesis()  
            self.obs = self.obs_trial_pairs()
            self.p_vals = self.compute_p_val()
        else:
            self.obs_all = self.obs_spatial_cat()
    
    def compute_p_val(self):
        p_vals = np.ones((len(self.geneList), len(self.geneList)), dtype = np.float16)
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
        genecount = pd.DataFrame.from_dict(Counter(self.curr_cell_df.index) , orient='index').reindex(self.geneList).fillna(0)
        num_trial_pairs = np.dot(genecount, genecount.T)   #n1*n2  #using numpy here now
        return genecount, num_trial_pairs

    def obs_trial_pairs(self): #debug for large d
        '''
        After counting all gene pairs within the distance threshold, create a pivot table to create a 2d represtation
        and fill in the missing genes.
        '''
        pairs = [(self.curr_cell_df.index[i], self.curr_cell_df.index[j]) for (i,j) in self.pairs]
        pairs = Counter(pairs)
        pairs = pd.Series(pairs).reset_index()
        if len(pairs):
            obs_df = pd.pivot_table(pairs, index = 'level_0', columns ='level_1', values = 0).fillna(0)
            col_na_genes = list(set(self.geneList).difference(obs_df.columns))
            row_na_genes = list(set(self.geneList).difference(obs_df.index))
            obs_df.loc[ : , col_na_genes] = 0
            #obs_df.loc[row_na_genes] = 0
            for row_na in row_na_genes:
                obs_df.loc[row_na] = 0
            obs_df = obs_df.reindex(index = self.geneList, columns = self.geneList)
        else:                                       #if no entry less than dist thresh
            #print('no entry less than dist thresh, total rna', self.curr_cell_df.shape[0])
            obs_df = pd.DataFrame(0, self.geneList, self.geneList)
        arr2 = np.triu(obs_df) + np.triu(obs_df,1).T   #Making matrix symmetric  #ASK
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
    def __init__(self, geneList, df, distance_threshold, min_gene_counts = 10):
        self.df = df
        self.distance_threshold = distance_threshold
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
            df = self.df[['uID', 'gene', 'absX', 'absY', 'absZ']]
            df = df[df.absZ == z_axis][['uID', 'gene', 'absX', 'absY']]
            if df.shape[0] > self.min_gene_counts:
                model_single_axis = ProximalPairs(self.geneList, df, self.dist_thresh)
                obs_all.append(model_single_axis.obs_df.values)
                num_trial_all.append(model_single_axis.num_trial_df.values)
                genecount_all.append(model_single_axis.genecount.iloc[:,0].values)  #df series doesn't return values with iloc directly
                prob_null, len_prob_null = model_single_axis.prob_null_hypothesis()
                prob_null_all.append(prob_null)
                len_prob_null_all.append(len_prob_null)

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
        return p_vals

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
    def __init__(self, threads = 1, min_intensity = 0, min_area = 0):
        self.min_intensity, self.min_area = min_intensity, min_area
        self.threads = threads

    def load_preprocessed_data(self, data):
        '''
        Load preprocessed data. Data should have the following columns - 
            ['gene', 'absX', 'absY', 'uID'] with 'gene' being the 1st column
            Arguments: 
                - data: (String) Path to dataframe in the required format.
        '''
        self.filename = data
        self.df = pd.read_csv(data, index_col=0)
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
    
    def load_preprocessed_data_filter_genes(self, data, threshold = 100):
        self.df = pd.read_csv(data, index_col=0)
        threshold_cells = self.df.groupby('uID').size()
        threshold_cells = threshold_cells[threshold_cells > threshold].index.values
        self.df = self.df.loc[self.df.uID.isin(threshold_cells)]
        threshold_genes = self.df.groupby('gene').size().sort_values(ascending=False).index.values[:100]
        self.df = self.df.loc[self.df.gene.isin(threshold_genes)]
        cell_ids = self.df.uID.unique()
        self.geneList = self.df.index.unique()
    
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
    
    def load_preprocessed_data_filter_genes(self, data, threshold = 100):
        self.df = pd.read_csv(data, index_col=0)
        threshold_cells = self.df.groupby('uID').size()
        threshold_cells = threshold_cells[threshold_cells > threshold].index.values
        self.df = self.df.loc[self.df.uID.isin(threshold_cells)]
        threshold_genes = self.df.groupby('gene').size().sort_values(ascending=False).index.values[:100]
        self.df = self.df.loc[self.df.gene.isin(threshold_genes)]
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
                pickle.dump(self.df.uID.unique(), f)
        return self.df.uID.unique()
    
    def _calculate_ProximalPairs(self, args):
        cell_num, cell_id = args[0], args[1]
        pp_model = ProximalPairs(geneList = self.geneList, df_loc = self.df_batch.loc[self.df_batch.uID == cell_id][['absX', 'absY']], distance_threshold = self.distance_threshold)
        return cell_num, pp_model.p_vals, pp_model.genecount.values.reshape(len(self.geneList))

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
    
    def _save_gene_list(self, filename):
        with open(filename, 'wb') as fp:
            pickle.dump(self.geneList, fp)
    
    def load_gene_list(self, filename):
        with open(filename, 'rb') as fp:
            self.geneList = pickle.load(fp)

    def run_ProximalPairs_parallel(self, distance_threshold, min_genecount, pval_matrix_name = None, gene_count_name = None):
        '''
        TODO Implement batched ProximalPairs calculation to decrease memory overhead.
            
            Arguments: 
                - expression_data: (String) Path to expression data csv file.
                - barcode_data: (String) Path to expression data csv file.
        '''
        self.distance_threshold = distance_threshold
        cell_ids = self.df.uID.unique()
        num_cells = len(cell_ids)
        print(f"Running PP now on {self.threads} threads, {mp.cpu_count()}")
        print("Number of cells: ", len(cell_ids), ", Number of Genes: ", len(self.geneList))
        start = timeit.default_timer()
        for n, i in enumerate(range(0, num_cells, 500)):
            print(f"Batc")
            self.df = pd.read_csv(self.filename, index_col=0)
            batch_cell_ids = cell_ids[i : i+500]
            valid_cell_data = []
            self.df_batch = self.df.loc[self.df.uID.isin(batch_cell_ids)]
            del self.df
            num_transcripts = 0
            for i, cell_id in enumerate(batch_cell_ids):
                num_transcripts += self.df_batch.loc[self.df_batch.uID == cell_id].shape[0]
                print(self.df_batch.loc[self.df_batch.uID == cell_id].shape)
                if self.df_batch[self.df_batch.uID == cell_id].shape[0] > min_genecount:
                    valid_cell_data.append([i, cell_id])
                else:
                    print(f"min genecount less than {min_genecount} for cell id {cell_id}, Skipping ...")
            check = timeit.default_timer()
            pool = mp.Pool(processes = self.threads)
            print(f"Running PP now on {self.threads} threads for cell batch {n}, {len(valid_cell_data)} cells, {num_transcripts} transcripts")
            results = pool.map(self._calculate_ProximalPairs, valid_cell_data)
            print(f"Done")  
            with open(pval_matrix_name[:-4]+f"_batch_{n}_results.pkl", 'wb') as fp:
                pickle.dump(results, fp)
            print(f"Time to complete PP for cell batch {n}: ", timeit.default_timer() - check)
        print(f"Cell-wise Proximal Pairs Time : {round(timeit.default_timer() - start, 2)} seconds")
    
    def run_ProximalPairs(self, distance_threshold, min_genecount, pval_matrix_name = None, gene_count_name = None):
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
        cell_ids = self.df.uID.unique()
        num_cells = len(cell_ids)
        print(f"Running PP now on {self.threads} threads, {mp.cpu_count()}")
        print("Number of cells: ", len(cell_ids), ", Number of Genes: ", len(self.geneList))
        start = timeit.default_timer()
        valid_cell_data = []
        self.df_batch = self.df.loc[self.df.uID.isin(cell_ids)]
        del self.df
        num_transcripts = 0
        for i, cell_id in enumerate(cell_ids):
            num_transcripts += self.df_batch.loc[self.df_batch.uID == cell_id].shape[0]
            print(self.df_batch.loc[self.df_batch.uID == cell_id].shape)
            if self.df_batch[self.df_batch.uID == cell_id].shape[0] > min_genecount:
                valid_cell_data.append([i, cell_id])
            else:
                print(f"min genecount less than {min_genecount} for cell id {cell_id}, Skipping ...")
        check = timeit.default_timer()
        pool = mp.Pool(processes = self.threads)
        print(f"Running PP now on {self.threads} threads for, {len(valid_cell_data)} cells, {num_transcripts} transcripts")
        results = pool.map(self._calculate_ProximalPairs, valid_cell_data)
        print(f"Done")  
        #with open(pval_matrix_name[:-4]+f"_batch_{n}_results.pkl", 'wb') as fp:
        #    pickle.dump(results, fp)
        self.all_pval = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = np.float16)
        self.all_gene_count = np.zeros((num_cells, len(self.geneList)), dtype = np.float16) 
        for cell_i_result in results:
            self.all_pval[cell_i_result[0]] = cell_i_result[1]
            self.all_gene_count[cell_i_result[0]] = cell_i_result[2]
        print(f"Time to complete PP : ", timeit.default_timer() - check)
        if pval_matrix_name:
            self._save_pval_matrix(pval_matrix_name)
        if gene_count_name:
            self._save_gene_count(gene_count_name)
        print(f"Cell-wise Proximal Pairs Time : {round(timeit.default_timer() - start, 2)} seconds")

    def _calculate_ProximalPairs3D(self, args):
        cell_num, cell_id = args[0], args[1]
        pp_model = ProximalPairs3D(geneList = self.geneList, df_loc = self.df.loc[self.df.uID == cell_id], distance_threshold = self.distance_threshold)
        return cell_num, pp_model.p_vals, pp_model.genecount.values.reshape(len(self.geneList))
    
    def run_ProximalPairs(self, distance_threshold, min_genecount, pval_matrix_name = None, gene_count_name = None):
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
        cell_ids = self.df.uID.unique()
        num_cells = len(cell_ids)
        print(f"Running PP now on {self.threads} threads, {mp.cpu_count()}")
        print("Number of cells: ", len(cell_ids), ", Number of Genes: ", len(self.geneList))
        start = timeit.default_timer()
        valid_cell_data = []
        self.df = self.df.loc[self.df.uID.isin(cell_ids)]
        num_transcripts = 0
        for i, cell_id in enumerate(cell_ids):
            num_transcripts += self.df.loc[self.df.uID == cell_id].shape[0]
            print(self.df.loc[self.df.uID == cell_id].shape)
            if self.df[self.df.uID == cell_id].shape[0] > min_genecount:
                valid_cell_data.append([i, cell_id])
            else:
                print(f"min genecount less than {min_genecount} for cell id {cell_id}, Skipping ...")
        check = timeit.default_timer()
        pool = mp.Pool(processes = self.threads)
        print(f"Running PP now on {self.threads} threads for, {len(valid_cell_data)} cells, {num_transcripts} transcripts")
        results = pool.map(self._calculate_ProximalPairs, valid_cell_data)
        print(f"Done")  
        self.all_pval = np.ones((num_cells, len(self.geneList), len(self.geneList)), dtype = np.float16)
        self.all_gene_count = np.zeros((num_cells, len(self.geneList)), dtype = np.float16) 
        for cell_i_result in results:
            self.all_pval[cell_i_result[0]] = cell_i_result[1]
            self.all_gene_count[cell_i_result[0]] = cell_i_result[2]
        print(f"Time to complete PP : ", timeit.default_timer() - check)
        if pval_matrix_name:
            self._save_pval_matrix(pval_matrix_name)
        if gene_count_name:
            self._save_gene_count(gene_count_name)
        print(f"Cell-wise Proximal Pairs Time : {round(timeit.default_timer() - start, 2)} seconds")

    def _save_globcolocal_csv(self, filename):
        self.global_coloc_df.to_csv(filename)
    
    def _save_expcolocal_csv(self, filename):
        self.expected_coloc_df.to_csv(filename)

    def _save_unstacked_pvals(self, filename, alpha_cellwise, min_transcript):
        present_cells = self._num_present_cells(min_transcript)
        obs = self._binarize_adj_matrix(alpha_cellwise).sum(axis=0)
        unstacked_global_pvals = self._unstack_df_both(obs, present_cells, alpha_cellwise)
        unstacked_global_pvals.to_csv(filename) 

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
        pairwise_p_val_df['g1g2'] = pairwise_p_val_df['gene_id1'] + ', ' + pairwise_p_val_df['gene_id2']
        pairwise_p_val_df = pairwise_p_val_df.set_index(['g1g2'])
        pairwise_p_val_df = pairwise_p_val_df.sort_values(by=['p_val_cond'])
        return pairwise_p_val_df
    
    def run_GlobalColocalization(self, g = None, alpha_cellwise = 0.01, min_transcript = 0, show_det_pairs = 0, high_precision = False, glob_coloc_name = None, exp_coloc_name = None, unstacked_pvals_name = None):
        '''
        Calculates global colocalization significance for each gene pair across all the cells. 
        Requires `run_ProximalPairs()` be run first to generate the p-value matrix for all cells. 
        Generates 3 dataframes both of size (num_genes, num_genes)
        Arguments: 
            - alpha_cellwise: (Float) pvalue signifcance threshold (>alpha_cellwise are converted to 1). Default = 0.05.
            - min_transcript: (Float) Gene expression lower threshold. Default = 0.
            - show_det_pairs: (Not used)
            - high_precision: (Boolean) High precision pvalue. Expect longer computer. Default = False.
            - glob_coloc_name: (String) (Optional) if provided saves global colocalization matrix as a csv at the input path.
            - exp_coloc_name: (String) (Optional) if provided saves expected colocalization matrix as a csv at the input path.
            - unstacked_pvals_name: (String) (Optional) if provided saves interpretable global colocalization matrix as a csv at the input path.
        '''
        print(f"Running GCL now on {self.threads} threads")
        start = timeit.default_timer()
        global_coloc_model = ConditionalGlobalColocalization(all_cell_pval = self.all_pval, transcript_count = self.all_gene_count, alpha_cellwise = alpha_cellwise, min_transcript = min_transcript, show_det_pairs = show_det_pairs, threads = self.threads, high_precision = high_precision)
        global_coloc, expected_coloc = global_coloc_model.global_colocalization()
        self.global_coloc_df = pd.DataFrame(global_coloc, index = self.geneList, columns = self.geneList)
        self.expected_coloc_df = pd.DataFrame(expected_coloc, index = self.geneList, columns = self.geneList)
        if glob_coloc_name:
            self._save_globcolocal_csv(glob_coloc_name)
        if exp_coloc_name:
            self._save_expcolocal_csv(exp_coloc_name)
        #self.geneList = pd.read_csv(g, index_col=0).index.unique()
        #self.global_coloc_df = pd.read_csv(glob_coloc_name, index_col=0)
        #self.expected_coloc_df = pd.read_csv(exp_coloc_name, index_col=0)
        if unstacked_pvals_name:
            self._save_unstacked_pvals(unstacked_pvals_name, alpha_cellwise, min_transcript)
        print(f"Cell-wise Global Colocalization Time : {round(timeit.default_timer() - start, 2)} seconds")