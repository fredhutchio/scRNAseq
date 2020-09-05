# Single Cell Gene Expression Pipeline: Cell Ranger

FIXME: intro

## Getting set up

Here is the short version of instructions found in the [remote computing](remote_compute.md) instructions.

- Log on to snail (if necessary while off campus): `ssh username@snail.fhcrc.org`
- Log on to rhino: `ssh username@rhino`
- Grab a node: `grabnode`
- Load software: `ml cellranger`

This material is modified from an original tutorial from 10 Genomics available
[here](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/tutorial_ov).

### Project organization

FIXME

    mkdir ~/yard/run_cellranger_mkfastq
    cd ~/yard/run_cellranger_mkfastq

### Downloading data

other data upload options [here](https://raw.githubusercontent.com/fredhutchio/tfcb_2019/master/lectures/lecture18/lecture18.md)

FIXME

    wget http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-1.2.0.tar.gz
    wget http://cf.10xgenomics.com/supp/cell-exp/cellranger-tiny-bcl-simple-1.2.0.csv
    tar -zxvf cellranger-tiny-bcl-1.2.0.tar.gz

    cat cellranger-tiny-bcl-simple-1.2.0.csv

    tree -L 2 cellranger-tiny-bcl-1.2.0

## Running Cell Ranger

FIXME: about Cell Ranger

### Demultiplexing Illumina base call files (BCL) with `cellranger mkfastq`

FIXME

    cellranger mkfastq --help

FIXME

    cellranger mkfastq --id=tutorial_walk_through \
    --run=cellranger-tiny-bcl-1.2.0 \
    --csv=cellranger-tiny-bcl-simple-1.2.0.csv

### Align sequence reads to reference using `cellranger count`

    mkdir ~/yard/run_cellranger_count
    cd ~/yard/run_cellranger_count

    wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar
    tar -xvf pbmc_1k_v3_fastqs.tar

    wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-3.0.0.tar.gz
    tar -zxvf refdata-cellranger-GRCh38-3.0.0.tar.gz

    cellranger count --help

    cellranger count --id=run_count_1kpbmcs \
    --fastqs=run_cellranger_count/pbmc_1k_v3_fastqs \
    --sample=pbmc_1k_v3 \
    --transcriptome=run_cellranger_count/refdata-cellranger-GRCh38-3.0.0

### Aggregate samples using `cellranger aggr`

    cellranger aggr --help

    mkdir ~/yard/run_cellranger_aggr
    cd ~/yard/run_cellranger_aggr

    wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_molecule_info.h5
    wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_molecule_info.h5

    pwd

Copy the output of the command above, use to create csv below

    nano pbmc_aggr.csv


    cellranger aggr --id=1k_10K_pbmc_aggr --csv=pbmc_aggr.csv

    ls -1 1k_10K_pbmc_aggr/outs/

### Rerun an analysis with different parameters using `cellranger reanalyze`

    mkdir ~/yard/run_cellranger_reanalyze
    cd ~/yard/run_cellranger_reanalyze


    cellranger reanalyze --help

    wget http://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5

create parameters file:

    nano reanalyze_10k_pbmcs.csv

    num_principal_comps,14
    max_clusters,15

    cellranger reanalyze --id=10k_pbmc_reanalyze_pc_clust --matrix=pbmc_10k_v3_filtered_feature_bc_matrix.h5 --params=reanalyze_10k_pbmcs.csv

    ls -1 10k_pbmc_reanalyze_pc_clust/outs/
