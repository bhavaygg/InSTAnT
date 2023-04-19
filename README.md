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

First we will initialise the Instant class object. This object will allow us to calculate the proximal pairs and find global colocalized genes. The primary argument is `threads` which controls the number of threads the program uses. The arguements `min_intensity` and `min_area` are used only for MERFISH data preprocessing and can be skipped otherwise.

```
obj = Instant(threads = threads, min_intensity = 10**0.75, min_area = 3)
```

To load MERFISH data, we use the function `preprocess_and_load_data`. `preprocess_and_load_data` is used to preprocess and load MERFISH data only. The final data is stored as a pandas DataFrame and the following columns `['gene', 'uID', 'absX', 'absY']`. Below is an example table - 

| gene | uID | absX | absY |
| :---         |     :---:      |          ---: |           ---: |
| AKAP11   |  2    | -1401.666    | -2956.618     |
| SIPA1L3  |  3       | -1411.692      |   -2936.609     |
| THBS1  |  925       | -764.6989      |   -1604.828    |

```
obj.preprocess_and_load_data(expression_data = f'data/u2os/rep3/data.csv', barcode_data = 'data/u2os/codebook.csv')
```
If the data has been preprocessed, we can just load it like below
```
obj.load_preprocessed_data(data = f'data/u2os_new/data_processed.csv')
```
The datafrmae is loaded in the object variable `df` and can be accesses through `ovj.df`. After the data has been loaded, we can now calculate the proximal gene pairs for each cell using the `run_ProximalPairs()` function. The following arguments are used by `run_ProximalPairs()`
  - `distance_threshold`: *(Integer)* Distance threshold at which to consider 2 genes proximal.
  - `min_genecount`: *(Integer)* Minimum number of transcripts in each cell.
  - `pval_matrix_name`: *(String)* *(Optional)* if provided saves pvalue matrix using pickle at the input path.
  - `gene_count_name`: *(String)* *(Optional)* if provided saves gene expression count matrix using pickle at the input path.
```
obj.run_ProximalPairs(distance_threshold = 4, min_genecount = 20, 
    pval_matrix_name = f"data/u2os/rep3/rep3_pvals.pkl", 
    gene_count_name = f"data/u2os/rep3/rep3_gene_count.pkl")
```

Next, we can find which gene pairs that significantly colocalized globally using the `run_GlobalColocalization()` function. The following arguments are used by `run_GlobalColocalization()`
  - `alpha_cellwise`: *(Float)* Pvalue signifcance threshold (>alpha_cellwise are converted to 1). Default = 0.05.
  - `min_transcript`: *(Float)* Gene expression lower threshold. Default = 0.
  - `high_precision`: *(Boolean*) High precision pvalue. Expect longer computer. Default = False.
  - `glob_coloc_name`: *(String)* *(Optional)* if provided saves global colocalization matrix as a csv at the input path.
  - `exp_coloc_name`: *(String)* *(Optional)* if provided saves expected colocalization matrix as a csv at the input path.
  - `unstacked_pvals_name`: *(String)* *(Optional)* if provided saves interpretable global colocalization matrix as a csv at the input path.
```
obj.run_GlobalColocalization(
    high_precision = False, 
    alpha_cellwise = 0.05,
    glob_coloc_name = f"data/u2os/rep3/global_colocalization.csv", 
    exp_coloc_name = f"data/u2os/rep3/expected_colocalization.csv", 
    unstacked_pvals_name = f"data/u2os/rep3/unstacked_global_pvals.csv")
```
Since, performing the ProximalPairs is a time consuming step but necessary for calculating Global Colocalization, we have added the support to load the outputs of ProximalPairs just so that you don't need to run them again. Apart from loading the data using `load_preprocessed_data()`, you need to use the functions below to load the pvalue matrix and the gene counts. 
```
obj.load_pval_matrix(f"data/u2os/rep3/rep3_pvals.pkl")
obj.load_gene_count(f"data/u2os/rep3/rep3_gene_count.pkl")
```

The final outputs are 3 csv files - 
  - Global colocalization csv :  contains pairwise significance value of d-colocalization
    | gene          | 5830417i10rik | Aatf     | Abcc1    | Abhd2    |
    | ------------- | ------------- | -------- | -------- | -------- |
    | 5830417i10rik | 1             | 0.317506 | 1        | 1        |
    | Aatf          | 0.317506      | 1        | 1        | 0.416612 |
    | Abcc1         | 1             | 1        | 1        | 0.055185 |
    | Abhd2         | 1             | 0.416612 | 0.055185 | 0.744798 |
  - Expected colocalization csv :  contains pairwise significance value of expected d-colocalization
    | gene          | 5830417i10rik | Aatf     | Abcc1    | Abhd2    |
    | ------------- | ------------- | -------- | -------- | -------- |
    | 5830417i10rik | 1.068986      | 1.978719 | 0.686994 | 0.999343 |
    | Aatf          | 1.978719      | 4.877825 | 1.69975  | 2.344674 |
    | Abcc1         | 0.686994      | 1.69975  | 0.795251 | 0.857116 |
    | Abhd2         | 0.999343      | 2.344674 | 0.857116 | 1.358516 |
  - Unstacked colocalization csv : contains pairwise significance value of d-colocalization in a interpretable format
    | g1g2          | gene_id1 | gene_id2 | p_val_cond | Expected coloc | Coloc. cells(Threshold<0.05) | Present cells | frac_cells |
    | ------------- | -------- | -------- | ---------- | -------------- | ---------------------------- | ------------- | ---------- |
    | Prpf8, Polr2a | Prpf8    | Polr2a   | \-2.5E-14  | 3.563221       | 112                          | 163           | 0.687117   |
    | Col1a1, Fn1   | Col1a1   | Fn1      | \-2.4E-14  | 5.350136       | 107                          | 179           | 0.597765   |
    | Fbln2, Fn1    | Fbln2    | Fn1      | \-1.7E-14  | 4.00996        | 74                           | 177           | 0.418079   |
    | Col1a1, Fbln2 | Col1a1   | Fbln2    | \-1.7E-14  | 4.629268       | 73                           | 177           | 0.412429   |
    | Col1a1, Bgn   | Col1a1   | Bgn      | \-1.6E-14  | 5.360849       | 71                           | 178           | 0.398876   |
