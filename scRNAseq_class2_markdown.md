# Analyzing Single-Cell RNA Sequencing Data in R

# Class 2: Preprocessing Cell Dataset Objects and Generating UMAP Plots

This class is based on materials provided by [Monocle 3](https://cole-trapnell-lab.github.io/monocle3/) and the Brotman
Baty Institute (via attended [tutorials](https://brotmanbaty.org/calendar/)).

## Objectives

In the last class we learned how to generate a cell dataset (CDS) object
and examined the structure of a CDS and how it can be modified
(i.e. adding columns). Today, we will be learning how to filter and do
an initial analysis on sc-RNA-seq data contained in a CDS.

By the end of this lesson, you should be able to:

  - determine appropriate filters for data contained in a CDS object
    (i.e. UMI and mitochondrial read cutoffs)
  - reduce the dimensions of a filtered CDS object using principle
    component analysis
  - generate cell clusters from data stored in a CDS object
  - create a UMAP plot to visualize your cell clusters

## Loading Required Packages

Load the packages that will be required for today’s class as follows.

``` r
# Load the monocle 3 package from last time so we have the set
# of tools designed for working with data stored in a CDS
# object
library(monocle3)

# Load the ggplot2 functions that we will use to generate 
# today's plots
library(ggplot2)

# Load the Matrix package from last time so that we have the 
# tools to more efficiently access the CDS object elements
library(Matrix)

# Load the dplyr package so that we can use some of its 
# functions for later dataframe manipulation
library(dplyr)
```

## Load the Data

Go to [this page](https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_microfluidic_emulsion_v2_/5968960)
and download the attached files by clicking the “Download All” button.
This should result in the download of a zip file called “5968960.zip”.

Navigate to the zip file and unzip it (can usually be done by double
clicking on the .zip file).

Open the new folder and unzip the “droplet.zip” file. Make sure to note
where this folder is stored on your computer.

The “genes.tsv” file in each of the datasets does not contain a column
called “gene\_short\_name”. However, we can still create a
cell\_data\_set object and rename the column that includes the gene
short names after creation.

FIXME (needs to link to github documents)

``` r
# Create a cell_data_set from the Lung-10X_P7_8 dataset

cds <- load_mm_data(mat_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/matrix.mtx", 
      feature_anno_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/genes.tsv", 
      cell_anno_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/barcodes.tsv")
```

Now, let's take a look at what's in the cds object.

``` r
# Take a peek at the cds

cds
```
![cds_view1](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/cds_view1.png)

Note that the "rowData names" indicates that there is a row called "V2." Note that all of the actual "rownames" are gene short names. Therefore, we want to rename "V2" as "gene_short_names" so that downstream Monocle 3 functions will work properly.

``` r
# Rename "V2" row as "gene_short_name"

names(rowData(cds))[names(rowData(cds))=="V2"] <- "gene_short_name"
```
Now, let's check that the rowData name "V2" was changed to "gene_short_name"
``` r
# Take another peek at the cds

cds
```
![cds_view2](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/cds_view2.png)

## Filtering a CDS Object

### Determining a Good UMI Cutoff

We can distinguish between “real cells” and detection/sequencing errors
based on the amount of RNAs (UMIs) detected in each cell. In order to
assess a good cutoff for the UMIs in our sample, we can generate a plot
showing the distribution of the UMIs in our sample.

Note that under colData names there’s a listing called “n.umi”. This
will be a column that lists the number of unique RNA molecules (or UMIs)
found in each individual cell. Let’s look at a snapshot of the column
data.

``` r
# Look at a snapshot of the column data found in the cds

colData(cds)
```
![coldata](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/colData_view.png)

Now let’s look at the distrubtion of our UMIs by generating a plot
showing how many cells contain a specific amount of UMIs (i.e. a density
plot).

NOTE: the “load\_mm\_data” function from Monocle automatically sets a
minimum cutoff for umis of 100. You can turn this off by using
“umi\_cutoff = 0” as an input to the function if needed.

FIXME (might need to introduce quick plot? also, need to add in plot
example images once this is run on data)

``` r
# Plot the distribution of UMIs 

qplot(colData(cds)$n.umi, geom="density")
```

![umi_plot1](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/UMI_standard.png)

The distribution of UMIs is not always very interpretable with a plot of
the direct counts. However, log transforming the count data can often
clear up the trend and make it easier to determine where to set a
cutoff.

``` r
# Plot the log10 distribution of UMIs

qplot(log10(colData(cds)$n.umi), geom="density")
```

![umi_log](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/UMI_log.png)

From both of these plots, we can conservatively say that most cells have
somewhere between 100 and 10,000 UMIs on the standard density plot or between
2 and 4 on the log10 scale.

Note that the minimum detected number of UMIs seems to be well above our
minimum cutoff, but this works for this data. You should never use less
than 100 UMIs as a cutoff (<100 indicates a generally bad sample), but
you can choose to increase it if necessary/desired.

Now, let’s apply the determined UMI cutoff.

``` r
# Filter out cells with more than 100 and less than 10000 UMIs

cds_goodumi <- cds[,pData(cds)$n.umi > 100 & pData(cds)$n.umi < 10000]
```
``` r
# Check that cells have been removed

cds_goodumi
```
![cds_goodumi](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/cds_umigood.png)

Note that we have now filtered out 93 cells based on the UMI cutoffs.

### Determining a Good Mitochondrial Gene Expression Cutoff

FIXME

(genes found in the mitochondrial dna don’t seem to be included in this
lung dataset, only mitochondrial activity associated genes)

One of the ways that we can tell a cell was alive and of good quality
going in is the number of mitochondrial genes it exhibits. To calculate
this, we first need a list of mitochondrial genes found in our organism,
then we need to calculate what percentage of the overall RNAs/UMIs in
the cell correspond to mitochondrial genes.

In the mouse genome, mitochondrial genes are prefaced with “mt” which
makes them easier to identify in our dataset.

### Normalizing the Data Based on ERCC Spike-In Controls

Published datasets often include a set of gene transcript identifiers at
the end of the gene list that have names beginning with “ERCC.” These
are a set of spike-in controls that can be used to normalize cell
RNA/UMI contents. This is done by calculating the amount of intrinsic
RNAs/UMIs measured over the number of ERCC transcripts/UMIs measured.

This is not a requirement and its usefulness is somewhat debated in the
context of single-cell transcriptomics (see [this
paper](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5549838/)), but is
often a useful tool.

FIXME (contemplating simply not using this as it is not in most current
monocle packages; rather size factor is utilized instead in a lot of the
calculations, but it would be useful for people to be aware of this if
they are trying to analyze other people’s data)

``` r
# Create a new cds with normalized transcript levels

# FIXME
```

## Preprocessing the Data Prior to Clustering

Given that single-cell RNA sequencing datasets are often very large, we
need to first determine how many principle components can adequately
explain most of the variance in our dataset. We can then use this
knowledge to reduce our matrix of data to only the critical components
which will significantly speed up downstream steps.

Let's start with calculating the first 100 principle components.

``` r
# Preprocess cds with first 100 principle components 

cds_preprocess100 <- preprocess_cds(cds_goodumi, num_dim = 100)
```
Now, let's evaluate how much variance is explained by each of these princple components.

``` r
# Look at the amount of variance explained by each parinciple component

plot_pc_variance_explained(cds_preprocess100)
```
![plot_pc_variance](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/plot_pc_variance.png)

## Clustering scRNAseq Data

In order to cluster the cells, we first have to reduce the
dimensionality of our data down to 2 dimensions (X,Y) in order to better
graph it (either with a tSNE or UMAP plot).

Both t-SNE and UMAP are ways of taking high dimensional data and visually representing them in 2D space. In general, UMAP is considered to be better because it preservers more spatial relationships. We will be using UMAP today.

Note that the reduce_dimension function automatically assumes principle component anlaysis (PCA) is the desired dimension reduction technique. There are others and you can set a different one by designating "preprocess_method" internally to the "reduce_dimension" function (i.e. cds_reddim <- reduce_dimension(cds_preprocess, preprocess_method="PCA") .

``` r
# Reduce dimensions of our data to 2D

cds_reddim <- reduce_dimension(cds_preprocess50)
```

Now, let's look at the spatial distribution of our cells that have been reduced in dimensions based on the first 50 principle components. We will use UMAP for visualization (alternatives to UMAP include t-SNE). 

``` r
# Look at the spatial distribution of our cells on a UMAP plot

plot_cells(cds_reddim)
```
![umap_plot1](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/umap_nocluster.png)

Note that cells will visually appear to be in clusters because they have similar RNA/UMI contents and therefore similar placement in UMAP space. However, clusters are not recognized as an entity by your computer yet even if they are visible by eye. Now, we want to tell our computer to recognize cells that are in similar locations in UMAP space as clusters as they are likely the same/similar cell types.

``` r
# Cluster cells that are spatially related in UMAP space

cds_clustered <- cluster_cells(cds_reddim)
```
``` r
# Plot the cells newly recognized as clusters

plot_cells(cds_clustered)
```
![umap_clustered](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/umap_cluster_noresindication.png)

Sometimes clusters that are very close together won’t be recognized as
separate with the automatic parameters in each of these functions. If
this happens you can run something like follows (note that our cluster
identification is good, so we should see the same thing).

``` r
# Cluster cells and change the resolution parameter
# Resolution ranges from 0 to 1
# Best to keep resolution between 0 and 1e-2 if possible

cds_clustered_res <- cluster_cells(cds_reddim, resolution=1e-5)
```
``` r
# Plot clustered cells with defined resolution

plot_cells(cds_clustered_res)
```
![umap_cluster_withres](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/umap_cluster_res.png)

## Wrapping Up

Today we learned:

  - how to evaluate which cells are good enough for analysis
  - how to preprocess cells for analysis
  - how to identify clusters in UMAP space
