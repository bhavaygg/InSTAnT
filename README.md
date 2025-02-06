InSTAnT
=========

[![DOI](https://zenodo.org/badge/588719195.svg)](https://zenodo.org/doi/10.5281/zenodo.10994621)
[![License: MIT](https://img.shields.io/badge/license-MIT-C06524)](https://github.com/bhavaygg/InSTAnT/blob/main/LICENSE.txt)
[![PyPI - Version](https://img.shields.io/pypi/v/hatch-fancy-pypi-readme.svg)](https://pypi.org/project/sc-instant/)
[![Downloads](https://pepy.tech/badge/sc-instant)](https://www.pepy.tech/projects/sc-instant)


**InSTAnT** is a toolkit to identify gene pairs which are d-colocalized
from single molecule measurement data e.g. MERFISH or SeqFISH. A gene
pair is d-colocalized when their transcripts are within distance d
across many cells.

This repository contains implementation of PP Test and CPB test and demo
on a U2OS dataset. The dataset can be downloaded from here (Moffit et
al., 2016, PNAS ) -
http://zhuang.harvard.edu/MERFISHData/data_for_release.zip 

Paper Preprint Link - https://pubmed.ncbi.nlm.nih.gov/36747718/


### UPDATE: Added support for AnnData input/output. Added Frequent Subgraph Mining.

We recommend using our `environment.yml` file to create a new conda environment to avoid issues with package incompatibility.

```
conda env create -f environment.yml
```
This will create a new conda environment with the name `instant` and has all dependencies installed. 

Alternatively, the package can be installed using pip.

```
pip install sc-instant
```

First, we will initialize the Instant class object. This object will allow us to calculate the proximal pairs and find global colocalized genes. The primary argument is `threads` which controls the number of threads the program uses. If you run into memory issues with the default settings, we suggest setting the `precision_mode` to `low`. The arguments `min_intensity` and `min_area` are used only for MERFISH data preprocessing and can be skipped otherwise.

```
from InSTAnT import Instant
obj = Instant(threads = threads, precision_mode = 'high', min_intensity = 10**0.75, min_area = 3)
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
Since, subcellular spatial transcriptomic data is generally present as `.csv` file containing all the transcripts, the primary input format is `.csv`. However, we have included the ability to format the data into an AnnData object and save the results of the subsequent analysis in that object. We convert the input file to an AnnData object and save it with the same name into the same directory. All the subsequent analysis will be updated and saved into the same file. We also provide the functionality to save any of the results seperately as individual files.

**Note** - The function also supports loading an AnnData object. Currently, we only accept files in `.h5ad` format. Since, Anndata objects are not natively designed for subcellular datasets, we expect the `.csv` file containing the transcript information to be present in `adata.uns['transcripts']`. If you wish to run differential colocalization, cell type labels are required, which are expected in `adata.obs` while for spatial modulation, cell locations are required, which are expected in `adata.uns['cell_locations']`. If these files are not present in the AnnData object, they must be supplied to the specific funtions seperately. If an AnnData object is provided during this loading, all the subsequent outputs are saved and updated in the specified file. 

```
obj.load_preprocessed_data(data = f'data/u2os_new/data.h5ad')
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
Proximal Pairs calculation has 3 variants - 
  - `run_ProximalPairs()` - Designed for 2D subcellular spatial transcriptomics data.
  - `run_ProximalPairs3D()` - Designed for 3D subcellular spatial transcriptomics data with continous z-axis.
  - `run_ProximalPairs3D_slice()` - Designed for 3D subcellular spatial transcriptomics data with discrete/sliced z-axis.
(**Note** - For AnnData objects. These results are stored in `adata.uns['pp_test_d{distance_threshold}_pvalues']`.)

All subsequent analysis require `run_ProximalPairs()` to be run first to generate the p-value matrix for all cells.

Next, we can use the output to find which gene pairs are significantly colocalized globally using the `run_GlobalColocalization()` function. The following arguments are used by `run_GlobalColocalization()`
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
    unstacked_pvals_name = f"data/u2os/rep3/unstacked_global_pvals.xlsx")
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

(**Note** - For AnnData objects, only the unstacked file is saved in `adata.uns['cpb_results']`.)

Next, we will use InSTAnT's spatial modulation analyses to find spatially modulated gene pairs. We use the `run_spatial_modulation()` function for this. The following arguments are used by `run_spatial_modulation()`
  - `inter_cell_distance`: *(Float)* Maximum distance between cells at which they are considered proximal.
  - `cell_locations`: *(String)* *(Optional)* Path to file contains locations for each cell. Should be in sorted order. If not provided, cell locations are expected to be provided in `adata.uns['cell_locations']` in the AnnData file specified during initialization.
  - `spatial_modulation_name`: *(String)* *(Optional)* Path and name of the output Excel file.
  - `alpha`: *(Float)* *(Optional)* p-value significance threshold (>alpha_cellwise are converted to 1). Default = 0.01.
  - `randomize`: *(Boolean)* *(Optional)* Shuffle cell locations. Default = False.
```
obj.run_spatial_modulation(f"data/u2os/rep3/cells_locations.csv", inter_cell_distance = 100, spatial_modulation_name = f"data/u2os/rep3/spatial_modulation.csv")
```
The `cell_locations` file should be a `.csv` file with the uID of the cells in sorted order and as the index of the file. The next 2 columns should be the x and y position of that cell respectively. (if provided in the AnnData file during initialization, cell locations are expected in `adata.uns['cell_locations']`)

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

(**Note** - For AnnData objects. These results are stored in `adata.uns['spatial_modulation']`.)

Lastly, we will calculate the cell type specificity of InSTAnT categorized d-colocalized gene pair. We call it Differential Colocalization and the function `run_differentialcolocalization()` is used for it. There are 3 different modes in which it can be run - 
  - `1va`: Compares colocalization for genes in the input cell type vs all other cell types.
  - `1v1`: Compares colocalization for genes in the input cell type 1 vs input cell type 2.
  - `ava`: Compares colocalization for genes for all cell types vs all other cell types.

Arguments: 
 - `cell_type`: *(String)* Cell type to calculate differential colocalization for. Is ignored if mode == "ava".
 - `cell_labels`: *(String)* *(Optional)* Path to file contains cell type for each cell. If not provided, cell labels are expected to be provided in `adata.obs` in the AnnData file specified during initialization. 
 - `file_location`: *(String)* *(Optional)* Directory in which to store output files. if mode == "ava", creates a new directory in this path to store all results.
 - `cell_type_2`: *(String)* (Optional) Cell type 2 to calculate differential colocalization for. Required if mode == "1v1".
 - `mode`: *(String)* Either "1va" (One cell type vs All cell types), "1v1" (One cell type vs one cell type) or "ava" (All cell types vs All cell types).
 - `alpha`: *(Float)* *(Optional)* p-value significance threshold (>alpha_cellwise are converted to 1). Default = 0.01.
 - `alpha_dc`: (Float) (Optional) pvalue signifcance threshold for unconditional differential colocalization. Default = 5e-6.
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

The function will create a folder named `folder_name` in the `file_location` location. The folder will contain a folder for each of the cell types selected to be analyzed and for each cell type contains a single file will all relevant outputs.
  - `{cell_type}_unstacked.csv`
    | g1,g2           | p_uncond | p_cond_g1 | p_cond_g2 | p_g1_expression | p_g2_expression | min_g1g2 | g1_rank | g2_rank | min_rank |
    | --------------- | -------- | --------- | --------- | --------------- | --------------- | -------- | ------- | ------- | -------- |
    | Slc17a6, Syt4   | 2.27E-75 | 6.75E-08  | 6.67E-64  | 5.8E-304        | 1.48E-17        | 5.8E-304 | 1       | 22      | 1        |
    | Cbln1, Slc17a6  | 2.58E-58 | 6.29E-14  | 5.01E-12  | 2.1E-131        | 5.8E-304        | 5.8E-304 | 2       | 1       | 1        |
    | Gabra1, Slc17a6 | 1.01E-53 | 6.75E-45  | 5.89E-06  | 4.39E-13        | 5.8E-304        | 5.8E-304 | 27      | 1       | 1        |
    | Cbln1, Gpr165   | 6.64E-39 | 1.5E-08   | 6.27E-43  | 2.1E-131        | 0.999642        | 2.1E-131 | 2       | 95      | 2        |
    | Cbln1, Syt4     | 4.29E-36 | 1.55E-10  | 1.35E-32  | 2.1E-131        | 1.48E-17        | 2.1E-131 | 2       | 22      | 2        |

(**Note** - For AnnData objects. These results are stored in `adata.uns['differential_colocalization']`. `adata.uns['differential_colocalization']` is a dictionary with keys based on the analysis done. For `1va` and `ava`, the key is `{cell_type}` and for `1v1`, the key is `{cell_type)_vs_{cell_type_2}`).

### UPDATE: Frequent Subgraph Mining.

Finds networks of genes colocalized in many cells using [gSpan](https://github.com/betterenvi/gSpan). The networks found are potential candidates for groups of cells and genes exhibiting interesting subcellular patterns, in this case colocalization.
Such colocalization networks could potentially be underlying factors for specific biological processes.

```
obj.run_fsm(n_vertices = 4, alpha = 0.001)
```

Arguments: 
  - `n_vertices`: *(int)* Minimum number of vertices in the network.
  - `alpha`: *(Float)* *(Float)* pvalue signifcance threshold (>alpha_cellwise are converted to 1). Default = 0.001.
  - `clique`: *(Boolean)* Forces networks found to be fully connected. Number of edges is number of vertices choose 2. Default = False.
  - `n_edges`: *(int) *(Optional)* Minimum number of edges in the network. Works only when clique == False.
  - `fsm_name`: *(String)* *(Optional)* Name of adata.uns key to save output in. Be sure to specify if `run_fsm` is ran multiple times because it will overwrite. Default = "nV{x}_cliques" where x is the number of vertices.

THe function will create a dataframe in `adata.uns` with all the gene pair networks found. Ny default the key is `nV{x}_cliques` where x is the number of vertices but using the attribute `fsm_name` it can be changed. 

## Citation 

```
@article{kumar2024intracellular,
  title={Intracellular spatial transcriptomic analysis toolkit (InSTAnT)},
  author={Kumar, Anurendra and Schrader, Alex W and Aggarwal, Bhavay and Boroojeny, Ali Ebrahimpour and Asadian, Marisa and Lee, JuYeon and Song, You Jin and Zhao, Sihai Dave and Han, Hee-Sun and Sinha, Saurabh},
  journal={Nature communications},
  volume={15},
  number={1},
  pages={7794},
  year={2024},
  publisher={Nature Publishing Group UK London}
}
```
