
## Making CDS objects for Single Cell Data Analysis

Two weeks before session, distribute instructions on how to install Monocle version 3 (https://cole-trapnell-lab.github.io/monocle3/docs/installation/), common troubleshooting (https://cole-trapnell-lab.github.io/monocle3/docs/installation/#troubleshooting) and answers to common questions and issues (https://github.com/cole-trapnell-lab/monocle3/issues).

One week before session, hold an office hour to answer any remaining questions regarding Monocle installation. Students are ready to proceed with the course once they can run library(monocle3) without any errors.
_____________________________

Eliza's Notes/Overview:
- review big picture what single cell data looks like (with images)
- install and load v2 data
- import the elemnts of a CDS from a regular matrix, convert to sparse then make CDS ""
  - review difference between sparse and dense matrix (use image)
  - review the format for new Monocle 3 input
  - make a table comparing old and new input for new_cell_data_set
- also can load directly if already sparse with "load_mm_data"
- can load directly from cell ranger output "load_cellranger_data"
______________________________

Today we will be reviewing how to create the data structures that contain all of our single cell RNA-sequencing data. These structured are called cell data sets or CDSs. These CDS objects are multi-dimensional and are made up of 3 individual files:
  1. a matrix of counts of sequenced reads or unique molecular identifiers (UMIs) by genes (rows) and cell (columns)
  2. a list of gene IDs that represents all genes quantified (correspond to matrix rows)
  3. a list of cell IDs that represents all cells measured (correspond to matrix columns) and any additional single cell data (i.e. replicates, experimental conditions) associated with those cell IDs

We can start by checking that Monocle has been installed correctly. If library(monocle3) runs without any errors, Monocle 3 has been successfully installed!

```{r}
library(monocle3)
```

Monocle is a toolkit in R to analyze single-cell gene expression experiments. 

To get started, we want to download the data we will be working with today. This is single cell data from mouse that has been processed with 10x Genomics (v2) single cell platform.

1. Download the data as a zip from Github (https://github.com/fredhutchio/scRNAseq/tree/monocle/lung_data/Lung-10X_P7_8)
2. Unzip the data 
3. Move the data to a new folder titled "scRNAseq_class"

There are multiple methods to create a CDS object in Monocle that we will cover today- if you are starting with:
1. full matrix of counts + cell IDs + gene IDs
2. sparse matrix of counts + cell IDs + gene IDs
3. 10X Genomics CellRanger output
  a. file path to the 'outs' file
  b. individual CellRanger output components
4. CDS object

First we will assume you are starting with the least processed form of single cell gene expression data, which is a matrix of counts- either reads of UMIs. To load the data into a CDS, you must convert the matrix into a sparse matrix. Instead of a matrix of values, a sparse matrix has three columns: row, column and value. The position of each non-zero value in the matrix is stored in this format. If you have a lot of data and a lot of that data is made of zeros, this is a much more efficient way to store it AND it is the format neccessary to make a CDS object for Monocle.

![matrix_types](http://www.btechsmartclass.com/data_structures/ds_images/Triplet_Representation_of_Sparse_Matrix.png)

```{r}

data <- as(as.matrix(expression_matrix), 'sparseMatrix')
head(data)
```

```{r}

```

```{r}

```
