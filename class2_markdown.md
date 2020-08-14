# Analyzing Single-Cell RNA Sequencing Data in R

# Class 2: Preprocessing Cell Dataset Objects and Generating UMAP Plots

This class is based on materials provided by Monocle 3 and the Brotman
Baty Institute.

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

Go to
<https://figshare.com/articles/dataset/Single-cell_RNA-seq_data_from_microfluidic_emulsion_v2_/5968960>
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

cds <- load_mm_data(mat_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/matrix.mtx", feature_anno_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/genes.tsv", cell_anno_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/barcodes.tsv")
```

``` r
# Now, take a look at the cds
cds
```
![cds_view1](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/cds_view1.png)

``` r
# Note that the "rowData names" includes a indication of
# something called "V2". We need to rename this to 
# gene_short_names since the rownames are all gene short 
# names.

names(rowData(cds))[names(rowData(cds))=="V2"] <- "gene_short_name"

#Now lets check that "V2" was changed to "gene_short_name"

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

![umi_plot1](https://github.com/fredhutchio/scRNAseq/blob/monocle/class2_figures/UMI_standard.png | width = 200)
