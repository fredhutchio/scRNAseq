# Analyzing Single-Cell RNA Sequencing Data in R 

# Class 4: Pseuodtime and Differential Expression Analysis

This class is based on materials provided by [Monocle 3](https://cole-trapnell-lab.github.io/monocle3/) and the Brotman
Baty Institute (via attended [tutorials](https://brotmanbaty.org/calendar/)).

## Objectives

In the last class we learned how to generate a classifier and use it to classify cells with Garnett. Today, we will learn how to determine the pseudotime associated with all of the clusters and how to look at the differential expression of genes between clusters. 

By the end of this lesson, you should be able to:

* determine the order of clusters in pseudotime
* identify genes differentially expressed in clusters of interest
* generate gene modules associated with clusters of interest

## Load Required Packages

Load the packages that will be required for today's class as follows.

```{r, message=FALSE}
# Load the monocle 3 package from last time so we have the set
# of tools designed for working with data stored in a CDS
# object
library(monocle3)

# Load the ggplot2 functions that we will use to generate 
# today's plots
library(ggplot2)

# Load the dplyr functions so that we will need downstream
# for filtering
library(dplyr)
```

## Reload CDS from Previous Lessons

We need to reload the cds from last time, which can be done as follows.

```{r, include=TRUE, message=FALSE}
# Load matrix file into cds
cds <- load_mm_data(mat_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/matrix.mtx",
                    feature_anno_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/genes.tsv",
                    cell_anno_path = "~/Downloads/5968960/droplet/Lung-10X_P7_8/barcodes.tsv")

# Rename V2 row to gene_short_name 
names(rowData(cds))[names(rowData(cds))=="V2"] <- "gene_short_name"

# Filter out cells based on umi levels
cds_goodumi <- cds[,pData(cds)$n.umi > 100 & pData(cds)$n.umi < 10000]

# Preprocess using 50 PCs
cds_preprocess50 <- preprocess_cds(cds_goodumi, num_dim = 50)

# Reduce dimensions to 2D of the preprocessed cds
cds_reddim <- reduce_dimension(cds_preprocess50)

# Cluster cells
cds_clustered <- cluster_cells(cds_reddim)

# Plot the cells to make sure everything looks like it did
# last time
plot_cells(cds_clustered)
```

![reloaded_cds](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/reloaded_UMAP.png)

## Differential Gene Expression

There are several ways in which we can evaluate differentially expressed genes between clusters of interest. 

* We can look at the top genes that are expressed in each cluster. 
* We can look at the expression of genes of interest across all clusters. These may or may not be differentially expressed, but can help identify cell types if they are.

We will explore both of these methods here.

### Looking at the Top Differentially Expressed Genes in Each Cluster

Let's determine which genes are differentially expressed in each cluster relative to other clusters.

```{r, include=TRUE, message=FALSE}
# Determine which genes are differentially expressed in clusters
diff_genes <- top_markers(cds_clustered,
                          group_cells_by="cluster",
                          reference_cells = 1000,
                          cores = 8)
```

The top_markers function generates various metrics that describe the expression of genes in each cluster. Let's take a look at what's included in diff_genes.

```{r, message=FALSE}
head(diff_genes)
```
![head_diff_genes](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/head_diff_genes.png)

There are multiple measures that we can assess. For example, we can look at the pseudo_R2 values which describes how well a given genes describe a cluster based on an underlying logistic regression model. 

Let's look at the top 5 genes that define each cluster based on pseudo_R2.

```{r, include=TRUE, message=FALSE}
# Identify the top 5 different genes found in each cluster
top5_specific_markers <- diff_genes %>%
                            filter(fraction_expressing >= 0.10) %>%
                            group_by(cell_group) %>%
                            top_n(5, pseudo_R2)
  
# Identify the gene names that are associated with each of 
# the top differentially expressed genes in each cluster
top5_specific_marker_ids <- unique(top5_specific_markers %>%
                                     pull(gene_id))

# Look at the expression level of each of the top genes
plot_genes_by_group(cds_clustered, 
                    top5_specific_marker_ids, 
                    group_cells_by="cluster", 
                    ordering_type="cluster_row_col", 
                    max.size=3)
```

![top_marker_genes](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/top5_specific_markers.png)

### Looking at the Expression of Genes of Interest

Let's look at the expression of some genes that are usually associated with lung tissues, since we know that our data includes lung cells.

```{r, message=FALSE}
# Look at lung-specific genes: Sftpc, sftpb, Scgb1a1, Ager, Slc34a2, Cldn18

plot_cells(cds_clustered, genes=c("Sftpc", "Sftpb", "Scgb1a1","Ager", "Slc34a2", "Cldn18"))
```

![specific_genes](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/genes_of_interest.png)

You can look at any gene in this way. Note that if a gene is not found within a dataset, the expression map will not be depicted. If this happens, make sure to double check the gene short name is correct as sometimes the colloquial name is not the recognized name (i.e. CD43 would actually be called Spn). 

Additionally, note that the representation of the gene name must match exactly to what is found in the CDS. For example, if you were to use all capital letters for the above genes, you would get an error stating that none of the genes are found in the dataset.

## Pseudotime Analysis

It is often useful to know how your cells link to each other across time. Is one cell cluster derived from another, for example? While we can't definitively say that one cell type gives rise to another, we can make a very, very good guess. The relationships between cells in "time" is measured therefore by inferring "pseudotime," where pseudotime is established by indicating a cell type as the precursor, and then finding the next most similar cell type, followed by the next most similar cell type, and so on. 

The first step in this process is learning a graph that represents the links between cell types/clusters.

```{r}
# Find the relationships of clusters to each other
cds_learn_graph <- learn_graph(cds_clustered)

# Plot the new version with a transposed graph
plot_cells(cds_learn_graph)
```
![learn_graph](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/learn_graph_plot.png)

Now, we have to tell the graph where the most likely starting point is. Pseudotime is best used when there is some sort of termporal component to your data. For example, you might be combining data collected from different embryo time points or maybe you're looking at a time course of samples with a drug added. With whichever temporal element you choose, it is typically a good idea to look at how cells are distributed in the clusters based on time points. 

Based on the distribution of your cells from those time points, you can then say which cluster contains the most cells from the earliest expected time point and tell Monocle to start there to determine lineages of cells that arise from that cell population. 

However, the dataset that we are using doesn't really have a temporal component. That being said, much like you plotted by cell_type in the Garnett class, you can always graph by an "embryo_time" or "time_point" component as long as you add that column to your cds as previously shown.

For now, let's move on to the next step and use a random cluster (we'll go with cluster 1) as our "start" point for pseudotime so that we can demonstrate how to order cells in pseudotime.


```{r}
# Use cluster 1 as a proxy "start" point for pseudotime
cds_order <- order_cells(cds_learn_graph)
```
![cds_order](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/order_cells_selected_node.png)

Note that in the above image, the "start" point has already been chosen and is highlighted in red. 

When the above window pops up, put your mouse over the "2" dot that is overlayed on cluster 1, then click it. Next, click the "Choose/unchoose" button on the left of the screen. Now, click "Done". You have now chosen a node over cluster 1 as the "start" for pseudotime.

We can now color our plot by pseudotime so we can see how the trajectory reflects time. 

```{r}
# Plot clusters with coloring indicating pseudotime
plot_cells(cds_order,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_branch_points = FALSE,
           graph_label_size=1.5)
```

![pseudotime_coloring](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/color_cells_by_pseudotime.png)

## Generating Gene Modules

One thing we might be interested in doing with our data is developing a sense of gene groups that have similar differential expression patterns across our data, either from cluster to cluster or across pseudotime. Groups of genes that have similar differential expression patterns are also called "gene modules." 

This can be done using a method called "graph-autocorrelation analysis" that was introduced as part of Monocle 3 (you can find more on the statistic that's used in this method on the Monocle 3 website linked at the beginning of this lesson; note that most of this code is only a slightly different version of that found on the Monocle 3 website). This method involves using a function called graph_test to determine genes that are differentially expressed between clusters or across pseudotime. Because clusters are the more relevant comparison in our dataset given the lack of temporal data, we will focus on genes differentially expressed between clusters.


```{r}
# Use graph_test to find differentially expressed genes between clusters
graph_test_result <- graph_test(cds_clustered, 
                                neighbor_graph="knn", 
                                cores=8)

# Determine the gene ids associated with those genes that are differentially
# expresed between clusters based on the statistic calculated by graph_test
# Note that rownames are all gene_short_names and the q_value cutoff is the 
# standard 0.05
deg_ids <- row.names(subset(graph_test_result, q_value < 0.05))
```

Now we can find gene modules containing genes that have similar expression patterns.

```{r}
# Find genes that have similar differential expression patterns across 
# clusters and store the result in a dataframe
# Note the resolution is the same as used to define the clusters
gene_modules <- find_gene_modules(cds_clustered[deg_ids,], resolution = 1e-2)

# Look at what's included in the new gene_modules dataframe
head(gene_modules)
```

![head_gene_modules](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/head_gene_modules.png)

Note that the new gene_modules dataframe includes a column called "module" containing a number. This is the number indicating which group of genes it most resembles in differential expression across the clusters. However, we likely want to know how these gene groups are expressed across clusters and what makes each group of genes different. For example, is module 1 expressed predominantly in endothelial cells whereas module 2 is predominantly downregulated in endothelial cells? 

Let's take a look at the expression of these groups of genes, or gene modules, in each cluster.

```{r}
# Create a tibble (something sort of like a list) containing a list of which clusters each cell belongs to
# Note that this particular step is most relevant if you've created a cds with only a subset of your clusters
cell_groups <- tibble::tibble(cell=row.names(colData(cds_clustered)), cell_group=clusters(cds_clustered)[colnames(cds_clustered)])

#Take a look at what this tibble looks like
head(cell_groups)
```
![head_cell_groups](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/head_cell_group.png)

We now have a list of barcodes/cells and which clusters they belong to. Now we can create a matrix containing the aggregated expression values associated with each gene module across the cells in each cluster.

```{r}
# Create matrix with the expression level that each group of genes has in each cluster
aggregate_exp <- aggregate_gene_expression(cds_clustered, gene_modules, cell_groups)

# Rename each of the rows in aggregate_exp with the proper Module name/number
row.names(aggregate_exp) <- stringr::str_c("Module", row.names(aggregate_exp))

# Renameeach of the columns in aggregrate_exp with the proper cluster name/number
colnames(aggregate_exp) <- stringr::str_c("Cluster", colnames(aggregate_exp))

# Take a peek at aggregate_exp and what it contains
head(aggregate_exp)
```
![head_agg_exp](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/head_agg_exp.png)

Now, let's finally look at the expression of each gene module in each cluster so that we can determine gene modules/groups of interest based on the cell type(s) we are most interested in.

```{r}
# Create heatmap showing the expression level of each gene module in each cluster
pheatmap::pheatmap(aggregate_exp, cluster_rows=TRUE, cluster_cols=TRUE, scale="column", clustering_method="ward.D2", fontsize=6)
```
![heatmap](https://github.com/fredhutchio/scRNAseq/blob/monocle/class4_figures/gene_module_heatmap.png)

Note that modules with more similar expression patterns across clusters are grouped closer together and clusters that have more similar expression patterns are closer together. 

If a module has an expression pattern of interest (i.e. maybe it's highly expressed in a cell type of interest) you can pull out the genes found in that module by looking back at the gene_modules dataframe.

FIX ME 
(maybe add in a brief example of how to do this; probably not wholly necessary since we've already shown them how to view that dataframe)

## Wrapping Up
Today we learned:

* how to look at the most differentially expressed genes in each cluster
* how to look at the expression of a gene across clusters
* how to determine pseudotime as inferred from a known starting cell population
* how to identify groups of genes with similar differential expression patterns across clusters (can also be extended to pseudotime)
