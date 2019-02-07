# magic-docker

[![Travis CI Build](https://api.travis-ci.com/KrishnaswamyLab/magic-docker.svg?branch=master)](https://travis-ci.com/KrishnaswamyLab/magic-docker)
[![Docker Build](https://img.shields.io/docker/pulls/scottgigante/magic.svg?style=flat)](https://cloud.docker.com/repository/docker/scottgigante/magic)

### Quick Start

```
docker pull scottgigante/magic
docker run scottgigante/magic --help
```

### Usage Example

```
docker run scottgigante/magic --filename https://github.com/KrishnaswamyLab/MAGIC/raw/master/data/HMLE_TGFb_day_8_10.csv.gz --gene-labels both --min-library-size 2000 --normalize --transform sqrt --knn 5 --decay 15 --output magic_output.csv
```

### Validation

You can check that PHATE is installed and running correctly with our [test dataset](https://raw.githubusercontent.com/KrishnaswamyLab/magic-docker/master/magic-validate.csv) as follows:

```
docker run scottgigante/magic --validate
```

which will compare the results with expected output stored in [magic-validate.csv](https://github.com/KrishnaswamyLab/magic-docker/blob/master/magic-validate.csv).

### Help

If you have any questions or require assistance using MAGIC, please visit <https://github.com/KrishnaswamyLab/MAGIC> or contact us at <https://krishnaswamylab.org/get-help>.

```
Run MAGIC for imputation of high-dimensional data.

optional arguments:
  -h, --help            show this help message and exit

Data IO:
  --filename FILENAME   Input data. Allowed types: csv, tsv, mtx, hdf5/h5 (10X format), directory/zip (10X format)
  --validate            Run MAGIC on a test dataset to ensure output is correct.
  --sparse              Use sparse data format
  --dense               Use dense data format
  --gene-names          Use gene name headers in data file (csv, tsv, fcs)
  --no-gene-names       Do not use gene names (csv, tsv, fcs, mtx)
  --gene-name-file FILE
                        Use gene name headers in FILE (csv, tsv, fcs, mtx)
  --cell-names          Use cell name headers in data file (csv, tsv, fcs)
  --no-cell-names       Do not use cell names (csv, tsv, fcs, mtx)
  --cell-name-file FILE
                        Use cell name headers in FILE (csv, tsv, fcs, mtx)
  --cell-axis {row,column}
                        States whether cells are on rows or columns (csv, tsv, mtx)
  --gene-labels {symbol,id,both}
                        Choice of gene labels for 10X data (dir, zip, hdf5)
  --genome GENOME       Genome name for 10X HDF5 data (hdf5)
  --metadata-channels CHANNEL [CHANNEL ...]
                        Names of channels to remove from fcs data (fcs)

Preprocessing:
  --min-library-size COUNTS
                        Filter cells with less than COUNTS counts
  --no-cell-filter      Do not filter cells
  --min-cells-per-gene CELLS
                        Filter genes with less than CELLS non-zero cells
  --no-gene-filter      Do not filter genes
  --normalize           Normalize cells by total UMI count (library size)
  --no-normalize        Do not normalize cells
  --transform {sqrt,log,arcsinh}
                        Sublinear data transformation function
  --no-transform        Do not transform data
  --pseudocount PCOUNT  Pseudocount to add to genes prior to log transform
  --cofactor COFACTOR   Factor by which to divide genes prior to arcsinh transform

Kernel Computation:
  -k KNN, --knn KNN     Number of nearest neighbors on which to build kernel
  -a DECAY, --decay DECAY
                        Sets decay rate of kernel tails
  --no-decay            Do not use alpha decay
  --pca N_PCA           Number of principal components to use for
                        neighborhoods
  --no-pca              Do not use PCA
  --knn-dist DISTANCE   Distance metric to use for calculating neighborhoods.
                        Recommended values are "euclidean" and "cosine"
  -t THREADS, --threads THREADS
                        Use THREADS threads. If -1 all CPUs are used
  --seed SEED           Integer random seed
  -v, --verbose         Print verbose output
  -q, --quiet           Do not print verbose output
  -vv, --debug          Print debugging output

MAGIC:
  --t-magic T           Level of diffusion for MAGIC
  --pca-only            Return PCA on the smoothed matrix
  --all-genes           Return the entire smoothed matrix
  --gene-list GENE [GENE ...]
                        List of genes to return from MAGIC, either as integer indices or column names.
  --output FILE         Output CSV file to save smoothed data matrix

For help, visit magic.readthedocs.io or krishnaswamylab.org/get-help
```
