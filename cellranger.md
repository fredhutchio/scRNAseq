# Single Cell Gene Expression Pipeline: Cell Ranger

FIXME: intro

## Logging on to rhino

FIXME: overview rhino

### On Fred Hutch's campus

Check to make sure you're connected to the Marconi wireless network. Open the program you'll be using to access Unix on your computer (Terminal on Mac or WSL/Git bash on Windows) and type the following (where `username` is your HutchNetID):

    ssh username@rhino

Hit enter to execute the command. The first time you execute this on a computer, you'll receive a response similar to:

> The authenticity of host 'rhino.fhcrc.org' can't be established.
> Are you sure you want to continue connecting (yes/no)?

Type `y` and press enter.

You'll be prompted to enter your username and then password (note that your cursor will not move as you type your password; this is normal!).

After you have successfully logged in, you'll see the information about the cluster printed to your screen. You'll be ready to start entering commands when you see a prompt like the following appear:

> username@rhino2:~$


### Off campus log-in

For more information about remote login, please see [this entry](https://sciwiki.fredhutch.org/scicomputing/access_methods/#access-via-a-remote-location) in the Fred Hutch Biomedical Data Science Wiki.

The short version: logging in off campus requires an additional step to connect to the campus network (where username is your HutchNetID):

    ssh username@snail.fhcrc.org

You'll see a message printed to the screen that starts with:

>  Welcome to the Fred Hutchinson Cancer Research Center

Then you can enter:

    ssh rhino

You'll then be able to interact with the cluster as if you were on campus.


## Loading software

Your compute environment is the collection of machinery, software, and networks on which you perform tasks. Managing environments is challenging, but there are tools available to help you stay on the right track.

Rhino has hundreds of pieces of software installed, so it uses a module system to make software available for use. To see what software you currently have loaded:

    module list

If you're looking for a specific package, you can search using terms related to the name of the software:

    module spider cellranger

If you see the package you want, you can load it with:

    module load cellranger

...or the abbreviated version, which is the same as above:

    ml cellranger

Be careful with multiple versions of the same software! You can load specific versions too (where XXX is the version number):

    ml cellranger/XXX

FIXME: Link to software available on rhino

FIXME: loading modules every time you log on


## Getting set up

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
