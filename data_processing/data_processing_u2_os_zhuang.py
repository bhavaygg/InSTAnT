import pandas as pd
import numpy as np
import time
import scipy.io
import pickle



class DataLoader():
    def __init__(self, exp_path, codebook_path, min_intensity = 0, min_area = 0):
        start_time = time.time()
        self.exp_path = exp_path
        self.codebook_path = codebook_path
        self.min_intensity, self.min_area = min_intensity, min_area
        self.df = self.load_data()
        
        
        print('time taken to load data', time.time() - start_time)

        

    def load_data(self):
        df = pd.read_csv(self.exp_path)
        
        codebook = self.load_codebook()
        df['geneName'] = df['barcode'].apply(lambda x: codebook.loc[x,'name'])
        self.geneList = df.geneName.unique()
        df = df[df.normalized_intensity > self.min_intensity]
        df = df[df.area >= self.min_area]
        df = df.rename(columns={'cell_id': 'uID', 'abs_x': 'absX', 'abs_y': 'absY'})
        # df = df.drop(['barcode', 'area','is_exact', 'normalized_intensity',  'distPeriphery', 'distNucleus'],axis=1)
        # df = df[df.is_exact==1]
        df = df.drop(['barcode', 'area','is_exact', 'normalized_intensity'],axis=1)  
        df = df.set_index('geneName')
        return df

    def load_codebook(self):
        codebook = pd.read_csv(self.codebook_path, converters={'bit_barcode': lambda x: str(x)})
        codebook['barcode'] = codebook['bit_barcode'].apply(lambda x: int(x[::-1],2))
        codebook = codebook[['name', 'barcode']]
        codebook = codebook.set_index('barcode')
        return codebook

    def save_cell_id_list(self, f=None):
        if f is not None:
            with open(f, 'wb') as f:
                pickle.dump(self.df.uID.unique(), f)
        return self.df.uID.unique()

        

class BoundaryLoader:
    def __init__(self):
        self.file_cyt_x = '/Users/anu/Desktop/1Nature_paper/Data/all_rep_zhuang_u2_os/rep3/cytoplasm_boundary_x_rep3.mat'
        self.file_cyt_y = '/Users/anu/Desktop/1Nature_paper/Data/all_rep_zhuang_u2_os/rep3/cytoplasm_boundary_y_rep3.mat'
        self.file_nuc_x = '/Users/anu/Desktop/1Nature_paper/Data/all_rep_zhuang_u2_os/rep3/nucleus_boundary_x_rep3.mat'
        self.file_nuc_y = '/Users/anu/Desktop/1Nature_paper/Data/all_rep_zhuang_u2_os/rep3/nucleus_boundary_y_rep3.mat'
        self.file_cellid = '/Users/anu/Desktop/1Nature_paper/Data/all_rep_zhuang_u2_os/rep3/all_cell_id_rep3.mat'
        self.cyt_boundary, self.nuc_boundary = self.all_cell_boundary()

    def cell_boundary(self,cell_id):
        return self.cyt_boundary[cell_id],self.nuc_boundary[cell_id]

    def all_cell_boundary(self):  #Improve this later
        
        cyt_x = scipy.io.loadmat(self.file_cyt_x)['cytoplasm_boundary_x']
        cyt_y = scipy.io.loadmat(self.file_cyt_y)['cytoplasm_boundary_y']
        nuc_x = scipy.io.loadmat(self.file_nuc_x)['nucleus_boundary_x']
        nuc_y = scipy.io.loadmat(self.file_nuc_y)['nucleus_boundary_y']
        mat_cellid = scipy.io.loadmat(self.file_cellid)['all_cell_id']

        cyt_boundary = {}
        nuc_boundary = {}
        for i in range(cyt_x.shape[0]):
            cyt_boundary[mat_cellid[i][0]] = np.stack((cyt_x[i][0].flatten(), cyt_y[i][0].flatten()), axis=1)
            nuc_boundary[mat_cellid[i][0]] = np.stack((nuc_x[i][0].flatten(), nuc_y[i][0].flatten()), axis=1)
        return cyt_boundary,nuc_boundary
