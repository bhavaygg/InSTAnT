import matplotlib.pyplot as plt
import pandas as pd
import time
import scipy
import scipy.io
import numpy as np

class DataLoader():
    def __init__(self, rna_loc_path,  gene_name_path, min_area = 0):
        start_time = time.time()
        self.rna_loc_path = rna_loc_path
        self.gene_name_path = gene_name_path
        self.min_area = min_area

        self.rna_loc_df, self.geneList = self.load_rna_loc_data()
        self.df, self.geneList = self.load_rna_loc_data()
        print('time taken to load data', time.time() - start_time)

    def load_rna_loc_data(self):
        '''loads rna loc for all cells and genenames 
        '''

        rna_loc_df = pd.read_csv(self.rna_loc_path)  #sep = '\t'
        col_name = ['barcode_id', 'abs_position_1', 'abs_position_2', 'abs_position_3', 'feature_id', 'in_feature', 'area']
        rna_loc_df = rna_loc_df.loc[:,col_name]

        rna_loc_df = rna_loc_df[rna_loc_df.area >= self.min_area]
        rna_loc_df = rna_loc_df.drop(['area'],axis=1)  

        rna_loc_df = rna_loc_df.rename(columns={"barcode_id": "geneID", "feature_id": "uID", "in_feature" : "in_nucleus", \
        'abs_position_1' : 'absX', 'abs_position_2' : 'absY', 'abs_position_3' : 'absZ' })  #Not true axes, just convention
        

        geneList = pd.read_csv(self.gene_name_path, header = None)[0].tolist()

        rna_loc_df['geneName'] = rna_loc_df.geneID.apply(lambda x: geneList[x-1])  #Indexing of list starts from 0


        return rna_loc_df, geneList

class BoundaryLoader():
    def __init__(self,  cyt_boundary_path, cells_loc_path):
        self.cells_loc_path = cells_loc_path
        self.cyt_boundary_path = cyt_boundary_path

    @staticmethod
    def str_bdr_to_float(bx):
        l1 = bx.values[0].split(';')[:-1]
        l1 = [float(i) for i in l1]
        return l1

    def cell_boundary(self,cell_id):
        cells_loc = pd.read_csv(self.cells_loc_path)
        bx, by = cells_loc[cells_loc.feature_ID == cell_id].boundaryX, cells_loc[cells_loc.feature_ID == cell_id].boundaryY
        bx, by = BoundaryLoader.str_bdr_to_float(bx), BoundaryLoader.str_bdr_to_float(by)
        return bx,by

    def cell_boundary_slice(self,cell_id,z_slice=3):
        cyt_boundary = pd.read_csv(self.cyt_boundary_path)
        cyt_boundary = cyt_boundary[cyt_boundary.feature_ID==cell_id]
        cyt_x,cyt_y = cyt_boundary.abs_x_boundary_3, cyt_boundary.abs_y_boundary_3
        cyt_x, cyt_y = BoundaryLoader.str_bdr_to_float(cyt_x), BoundaryLoader.str_bdr_to_float(cyt_y)
        return cyt_x, cyt_y

    