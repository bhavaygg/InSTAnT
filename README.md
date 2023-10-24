InSTAnT
=========

**InSTAnT** is a toolkit to idetify gene pairs which are d-colocalized
from single molecule measurement data e.g. MERFISH or SeqFISH. A gene
pair is d-colocalized when their transcripts are within distance d
across many cells.

This repository contains implementation of PP Test and CPB test and demo
on a U2OS dataset. The dataset can be downloaded from here (Moffit et
al., 2016, PNAS ) -
http://zhuang.harvard.edu/MERFISHData/data_for_release.zip 

We recommend using our `environment.yml` file to create a new conda environment to avoid issues with package incompatibility.

```
conda env create -f environment.yml
```
This will create a new conda environment with the name `instant` and has all dependancies installed. 

First we will initialise the Instant class object. This object will allow us to calculate the proximal pairs and find global colocalized genes. The primary argument is `threads` which controls the number of threads the program uses. The arguements `min_intensity` and `min_area` are used only for MERFISH data preprocessing and can be skipped otherwise.

```
obj = Instant(threads = threads, min_intensity = 10**0.75, min_area = 3)
```

To load MERFISH data, we use the function `preprocess_and_load_data()`. `preprocess_and_load_data()` is used to preprocess and load MERFISH data ***only***. The final data is stored as a pandas DataFrame and the following columns `['gene', 'uID', 'absX', 'absY']`. Below is an example table for a 2D data (if the data is 3D, `absZ` column is also expected) - 

| gene | uID | absX | absY |
| :---         |     :---:      |          ---: |           ---: |
| AKAP11   |  2    | -1401.666    | -2956.618     |
| SIPA1L3  |  3       | -1411.692      |   -2936.609     |
| THBS1  |  925       | -764.6989      |   -1604.828    |

```
obj.preprocess_and_load_data(expression_data = f'data/u2os/rep3/data.csv', barcode_data = 'data/u2os/codebook.csv')
```
If the data has been preprocessed, we can load it like below.
```
obj.load_preprocessed_data(data = f'data/u2os_new/data_processed.csv')
```

The dataframe is loaded in the object variable `df` and can be accessed through `obj.df`. After the data has been loaded, we can calculate each cell's proximal gene pairs using the `run_ProximalPairs()` function. The following arguments are used by `run_ProximalPairs()`
  - `distance_threshold`: *(Integer)* Distance threshold at which to consider 2 genes proximal.
  - `min_genecount`: *(Integer)* Minimum number of transcripts in each cell.
  - `pval_matrix_name`: *(String)* *(Optional)* if provided saves pvalue matrix using pickle at the input path.
  - `gene_count_name`: *(String)* *(Optional)* if provided saves gene expression count matrix using pickle at the input path.
```
obj.run_ProximalPairs(distance_threshold = 4, min_genecount = 20, 
    pval_matrix_name = f"data/u2os/rep3/rep3_pvals.pkl", 
    gene_count_name = f"data/u2os/rep3/rep3_gene_count.pkl")
```
Since calculating proximal pairs is the most time-consuming function, it is recommended to save the matrices for further analysis using our other features. Proximal Pairs calculation has 3 variants - 
  - `run_ProximalPairs()` - Designed for 2D subcellular spatial transcriptomics data.
  - `run_ProximalPairs3D()` - Designed for 3D subcellular spatial transcriptomics data with continuos z-axis.
  - `run_ProximalPairs3D_slice()` - Designed for 3D subcellular spatial transcriptomics data with discrete/sliced z-axis.

Next, we can find which gene pairs are significantly colocalized globally using the `run_GlobalColocalization()` function. The following arguments are used by `run_GlobalColocalization()`
  - `alpha_cellwise`: *(Float)* Pvalue signifcance threshold (>alpha_cellwise are converted to 1). Default = 0.05.
  - `min_transcript`: *(Float)* Gene expression lower threshold. Default = 0.
  - `high_precision`: *(Boolean*) High precision pvalue. Expect a longer compute time. Default = False.
  - `glob_coloc_name`: *(String)* *(Optional)* if provided, saves the global colocalization matrix as a CSV at the input path.
  - `exp_coloc_name`: *(String)* *(Optional)* if provided, saves the expected colocalization matrix as a CSV at the input path.
  - `unstacked_pvals_name`: *(String)* *(Optional)* if provided, saves a more interpretable global colocalization matrix as a CSV at the input path.
```
obj.run_GlobalColocalization(
    high_precision = False, 
    alpha_cellwise = 0.05,
    glob_coloc_name = f"data/u2os/rep3/global_colocalization.csv", 
    exp_coloc_name = f"data/u2os/rep3/expected_colocalization.csv", 
    unstacked_pvals_name = f"data/u2os/rep3/unstacked_global_pvals.csv")
```
Since performing the ProximalPairs is a time-consuming step but necessary for calculating Global Colocalization, we have added the support to load the outputs of ProximalPairs just so that you don't need to rerun them. Apart from loading the data using `load_preprocessed_data(),` you must use the functions below to load the p-value matrix and the gene counts. 
```
obj.load_pval_matrix(f"data/u2os/rep3/rep3_pvals.pkl")
obj.load_gene_count(f"data/u2os/rep3/rep3_gene_count.pkl")
```

The final outputs are 3 CSV files - 
  - Global colocalization:  contains pairwise significance value of d-colocalization
    | gene          | 5830417i10rik | Aatf     | Abcc1    | Abhd2    |
    | ------------- | ------------- | -------- | -------- | -------- |
    | 5830417i10rik | 1             | 0.317506 | 1        | 1        |
    | Aatf          | 0.317506      | 1        | 1        | 0.416612 |
    | Abcc1         | 1             | 1        | 1        | 0.055185 |
    | Abhd2         | 1             | 0.416612 | 0.055185 | 0.744798 |
  - Expected colocalization:  contains pairwise significance value of expected d-colocalization
    | gene          | 5830417i10rik | Aatf     | Abcc1    | Abhd2    |
    | ------------- | ------------- | -------- | -------- | -------- |
    | 5830417i10rik | 1.068986      | 1.978719 | 0.686994 | 0.999343 |
    | Aatf          | 1.978719      | 4.877825 | 1.69975  | 2.344674 |
    | Abcc1         | 0.686994      | 1.69975  | 0.795251 | 0.857116 |
    | Abhd2         | 0.999343      | 2.344674 | 0.857116 | 1.358516 |
  - Unstacked colocalization: contains pairwise significance value of d-colocalization in an interpretable format
    | g1g2          | gene_id1 | gene_id2 | p_val_cond | Expected coloc | Coloc. cells(Threshold<0.05) | Present cells | frac_cells |
    | ------------- | -------- | -------- | ---------- | -------------- | ---------------------------- | ------------- | ---------- |
    | Prpf8, Polr2a | Prpf8    | Polr2a   | \-2.5E-14  | 3.563221       | 112                          | 163           | 0.687117   |
    | Col1a1, Fn1   | Col1a1   | Fn1      | \-2.4E-14  | 5.350136       | 107                          | 179           | 0.597765   |
    | Fbln2, Fn1    | Fbln2    | Fn1      | \-1.7E-14  | 4.00996        | 74                           | 177           | 0.418079   |
    | Col1a1, Fbln2 | Col1a1   | Fbln2    | \-1.7E-14  | 4.629268       | 73                           | 177           | 0.412429   |
    | Col1a1, Bgn   | Col1a1   | Bgn      | \-1.6E-14  | 5.360849       | 71                           | 178           | 0.398876   |


Next, we will use InSTAnT's spatial modulation analyses to find spatially modulated gene pairs. We use the `run_spatial_modulation()` function for this. The following arguments are used by `run_spatial_modulation()`
  - `cell_locations`: *(String)* Path to a CSV file contains locations for each cell. It should be in sorted order.
  - `inter_cell_distance`: *(Float)* Maximum distance between cells at which they are considered proximal.
  - `spatial_modulation_name`: *(String)* Path and name of the output Excel file.
  - `alpha`: *(Float)* *(Optional)* p-value significance threshold (>alpha_cellwise are converted to 1). Default = 0.01.
  - `randomize`: *(Boolean)* *(Optional)* Shuffle cell locations. Default = False.
```
obj.run_spatial_modulation(f"data/u2os/rep3/cells_locations.csv", inter_cell_distance = 100, spatial_modulation_name = f"data/u2os/rep3/spatial_modulation.csv")
```
The `cell_locations` file should be a `CSV` file with the uID of the cells in sorted order and as the index of the file. The next 2 columns should be the x and y position of that cell respectively.

| uID  | x_centroid | y_centroid |
| ---- | ---------- | ---------- |
| 1054 | 5785.578   | 5699.618   |
| 1059 | 5757.25    | 5770.57    |
| 1067 | 5781.67    | 5238.706   |
| 1068 | 5784.023   | 5177.076   |
| 1069 | 5772.406   | 5173.247   |
| 1071 | 5757.652   | 5815.921   |

The output is an Excel file containing the gene pairs and the log-likelihood ratio of their spatial modulation. Below is an example -

| g1g2           | gene_id1 | gene_id2 | llr      | w_h1     | p_g_h1   | p_g_h0   |
| -------------- | -------- | -------- | -------- | -------- | -------- | -------- |
| MALAT1, MALAT1 | MALAT1   | MALAT1   | 113.2089 | 0.524039 | 0.803285 | 0.805375 |
| FASN, TLN1     | FASN     | TLN1     | 61.48302 | 0.404623 | 0.808779 | 0.808774 |
| COL5A1, THBS1  | COL5A1   | THBS1    | 58.72194 | 0.391671 | 0.729187 | 0.733704 |
| MALAT1, SRRM2  | MALAT1   | SRRM2    | 47.82571 | 0.353254 | 0.410644 | 0.409021 |
| COL5A1, FBN2   | COL5A1   | FBN2     | 47.62671 | 0.348113 | 0.574713 | 0.576769 |


Lastly, we will calculate the cell type specificity of InSTAnT categorized d-colocalized gene pair. We call it Differential Colocalization and the function `run_differentialcolocalization()` is used for it. There are 3 different modes in which it can be run - 
  - `1va`: Compares colocalization for genes in the input cell type vs all other cell types.
  - `1v1`: Compares colocalization for genes in the input cell type 1 vs input cell type 2.
  - `ava`: Compares colocalization for genes for all cell types vs all other cell types.
Requires `run_ProximalPairs()` to be run first to generate the p-value matrix for all cells. 
Arguments: 
            - `cell_type`: *(String)* Cell type to calculate differential colocalization for. Is ignored if mode == "ava".
            - `cell_labels`: *(String)* Path to the CSV file contains the cell type for each cell.
            - `file_location`: *(String)* Directory in which to store output files. if mode == "ava", creates a new directory in this path to store all results.
            - `cell_type_2`: *(String)* (Optional) Cell type 2 to calculate differential colocalization for. Required if mode == "1v1".
            - `mode`: *(String)* Either "1va" (One cell type vs All cell types), "1v1" (One cell type vs one cell type) or "ava" (All cell types vs All cell types).
            - `alpha`: *(Float)* *(Optional)* p-value significance threshold (>alpha_cellwise are converted to 1). Default = 0.01.
            - `folder_name`: *(String) *(Optional)* if mode == "ava", folder name inside the specified path in which to store results for each cell type. Default = "differential_colocalization".
```
obj.run_differentialcolocalization(cell_type = None, mode = "a2a", 
     cell_labels = f"data/u2os/rep3/cell_labels.csv", 
     file_location = f"data/u2os/rep3/",
     folder_name = "differential_colocalization")
```
The `cell_labels.csv` file must have 2 columns `uID` and `cell_type` which denote the cell number/ID and the type of the cell respectively.

| uID   | cell_type |
| ----- | --------- |
| 83434 | 27        |
| 83431 | 27        |
| 83447 | 27        |
| 7469  | 27        |
| 91749 | 27        |
| 83131 | 27        |

The function will create a folder named `folder_name` in the `file_location` location. The folder will contain a folder for each of the cell types selected to be analyzed and for each cell type contains 4 files.
  - `{cell_type}_conditional.csv`
    | gene          | 2010300C02Rik | Acsbg1 | Acta2 | Acvrl1 | Adamts2 |
    | ------------- | ------------- | ------ | ----- | ------ | ------- |
    | 2010300C02Rik | 1             | 1      | 1     | 1      | 1       |
    | Acsbg1        | 1             | 1      | 1     | 1      | 1       |
    | Acta2         | 1             | 1      | 1     | 1      | 1       |
    | Acvrl1        | 1             | 1      | 1     | 1      | 1       |
    | Adamts2       | 1             | 1      | 1     | 1      | 1       |
    | Adamtsl1      | 1             | 1      | 1     | 1      | 1       |
  - `{cell_type}_genemarkers.csv`
    | gene     | 0        |
    | -------- | -------- |
    | Tacr1    | 0        |
    | Tanc1    | 0.023006 |
    | Th       | 0        |
    | Thsd7a   | 0        |
    | Tle4     | 1.34E-19 |
    | Tmem132d | 7.33E-05 |
  - `{cell_type}_unconditional.csv`
    | gene          | 2010300C02Rik | Acsbg1   | Acta2    | Acvrl1   |
    |  ------------- | ------------- | -------- | -------- | -------- |
    | 2010300C02Rik | 0.768928      | 0.691713 | 0.957706 | 1        |
    | Acsbg1        | 0.691713      | 0.987682 | 0.139965 | 1        |
    | Acta2         | 0.957706      | 0.139965 | 0.536099 | 1        |
    | Acvrl1        | 1             | 1        | 1        | 0.139965 |
    | Adamts2       | 1             | 1        | 1        | 1        |
    | Adamtsl1      | 1             | 0.519267 | 1        | 1        |
  - `{cell_type}_unstacked.csv`
    | g1g2                 | gene_id1      | gene_id2 | ct_cond  | Unconditional | Conditional cells(Threshold<0.01) | Present cells | frac_cells |
    | -------------------- | ------------- | -------- | -------- | ------------- | --------------------------------- | ------------- | ---------- |
    | Epha4, Prox1         | Epha4         | Prox1    | 2.26E-14 | 7.81E-17      | 98                                | 4723          | 0.02075    |
    | 2010300C02Rik, Pdzd2 | 2010300C02Rik | Pdzd2    | 4.91E-13 | 2.68E-13      | 44                                | 4054          | 0.010853   |
    | Bhlhe22, Igfbp5      | Bhlhe22       | Igfbp5   | 1.03E-12 | 1.36E-15      | 52                                | 4623          | 0.011248   |
    | Nrn1, Plekha2        | Nrn1          | Plekha2  | 2.76E-12 | 3.07E-12      | 51                                | 4712          | 0.010823   |
    | 2010300C02Rik, Prox1 | 2010300C02Rik | Prox1    | 6.39E-12 | 5.03E-12      | 59                                | 4753          | 0.012413   |
