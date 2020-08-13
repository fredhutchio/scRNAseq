
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
Next, we want to download the data we will be working with today. This is single cell data from mouse that has been processed with 10x Genomics (v2) single cell platform.

1. Download the data as a zip from Github (https://github.com/fredhutchio/scRNAseq/tree/monocle/lung_data/Lung-10X_P7_8)
2. Unzip the data 
3. Move the data to a new folder titled "scRNAseq_class"

There are four methods to create a CDS object in Monocle:
1. load a full matrix of counts, list of 
2. load a sparse matrix of counts
3. 
4. 

![matrix_types](http://www.btechsmartclass.com/data_structures/ds_images/Triplet_Representation_of_Sparse_Matrix.png)

1) starting with a regular matrix, 2) starting with a sparse matrix, 3) with cell ranger output directory and 4) with load_mm_data

```{r}

```

```{r}

```

```{r}

```
