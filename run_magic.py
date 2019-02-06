import pandas as pd
import scprep
import magic

from utils import check_filetype, check_load_args, check_transform_args


def run_magic_from_file(
        filename,
        # data loading params
        sparse=True,
        gene_names=None,
        cell_names=None,
        cell_axis=None,
        gene_labels=None,
        allow_duplicates=None,
        genome=None,
        metadata_channels=None,
        # filtering params
        min_library_size=2000,
        min_cells_per_gene=10,
        # normalization params
        library_size_normalize=True,
        transform='sqrt',
        pseudocount=None,
        cofactor=None,
        # kernel params
        knn=5,
        decay=15,
        n_pca=100,
        knn_dist='euclidean',
        n_jobs=1,
        random_state=42,
        verbose=1,
        # magic params
        t_magic='auto',
        genes=None,
        # output params
        output='magic.csv'):
    """Run MAGIC on a file

    Parameters
    ----------
    filename : str
        Allowed types: csv, tsv, mtx, hdf5/h5 (10X format),
        directory/zip (10X format)
    sparse : bool (recommended: True for scRNAseq, False for CyTOF)
        Force data sparsity. If `None`, sparsity is determined by data type.
    gene_names : str, list or bool
        Allowed values:
        - if filetype is csv or fcs, `True` says gene names are data
        headers, `str` gives a path to a separate csv or tsv file containing
        gene names, list gives an array of gene names, `False` means
        no gene names are given
        - if filetype is mtx, `str` gives a path to a separate csv or tsv file
        containing gene names, list gives an array of gene names, or `False`
        means no gene names are given
        - if filetype is hdf5, h5, directory or zip, must be `None`.
    cell_names : str, list or bool
        Allowed values:
        - if filetype is csv or fcs, `True` says cell names are data
        headers, `str` gives a path to a separate csv or tsv file containing
        cell names, list gives an array of cell names, `False` means
        no cell names are given
        - if filetype is mtx, `str` gives a path to a separate csv or tsv file
        containing cell names, list gives an array of cell names, or `False`
        means no gene names are given
        - if filetype is hdf5, h5, directory or zip, must be `None`.
    cell_axis : {'row', 'column'}
        States whether cells are on rows or columns. If cell_axis=='row',
        data is of shape [n_cells, n_genes]. If cell_axis=='column', data is of
        shape [n_genes, n_cells]. Only valid for filetype mtx and csv
    gene_labels : {'symbol', 'id', 'both'}
        Choice of gene labels for 10X data. Recommended: 'both'
        Only valid for directory, zip, hdf5, h5
    allow_duplicates : bool
        Allow duplicate gene names in 10X data. Recommended: True
        Only valid for directory, zip, hdf5, h5
    genome : str
        Genome name. Only valid for hdf5, h5
    metadata_channels : list of str (recommended: ['Time', 'Event_length', 'DNA1', 'DNA2', 'Cisplatin', 'beadDist', 'bead1'])
        Names of channels in fcs data which are not real measurements.
        Only valid if datatype is fcs.
    min_library_size : int or `None`, optional (default: 2000)
        Cutoff for library size normalization. If `None`,
        library size filtering is not used
    min_cells_per_gene : int or `None`, optional (default: 10)
        Minimum non-zero cells for a gene to be used. If `None`,
        genes are not removed
    library_size_normalize : `bool`, optional (default: True)
        Use library size normalization
    transform : {'sqrt', 'log', 'arcsinh', None}
        How to transform the data. If `None`, no transformation is done
    pseudocount : float (recommended: 1)
        Number of pseudocounts to add to genes prior to log transformation
    cofactor : float (recommended: 5)
        Factor by which to divide genes prior to arcsinh transformation
    knn : int, optional, default: 10
        number of nearest neighbors on which to build kernel
    decay : int, optional, default: 15
        sets decay rate of kernel tails.
        If None, alpha decaying kernel is not used
    n_pca : int, optional, default: 100
        Number of principal components to use for calculating
        neighborhoods. For extremely large datasets, using
        n_pca < 20 allows neighborhoods to be calculated in
        roughly log(n_samples) time.
    knn_dist : string, optional, default: 'euclidean'
        recommended values: 'euclidean', 'cosine'
        Any metric from `scipy.spatial.distance` can be used
        distance metric for building kNN graph.
    n_jobs : integer, optional, default: 1
        The number of jobs to use for the computation.
        If -1 all CPUs are used. If 1 is given, no parallel computing code is
        used at all, which is useful for debugging.
        For n_jobs below -1, (n_cpus + 1 + n_jobs) are used. Thus for
        n_jobs = -2, all CPUs but one are used
    random_state : integer or numpy.RandomState, optional, default: None
        The generator used to initialize random PCA
        If an integer is given, it fixes the seed
        Defaults to the global `numpy` random number generator
    verbose : `int` or `boolean`, optional (default: 1)
        If `True` or `> 0`, print status messages
    t_magic : int, optional, default: 'auto'
        power to which the diffusion operator is powered for MAGIC.
        This sets the level of diffusion. If 'auto', t is selected
        according to the Procrustes disparity of the diffused data
    genes : list or {"all_genes", "pca_only"}, optional (default: None)
        List of genes to return from MAGIC,
        either as integer indices or column names
        if input data is a pandas DataFrame. If "all_genes", the entire
        smoothed matrix is returned. If "pca_only", PCA on the smoothed
        data is returned. If None, the entire matrix is also
        returned, but a warning may be raised if the resultant matrix
        is very large.
    output : str, optional (default: 'magic.csv')
        Output CSV file to save smoothed data matrix
    """
    # check arguments
    filetype = check_filetype(filename)
    load_fn, load_kws = check_load_args(filetype,
                                        sparse=sparse,
                                        gene_names=gene_names,
                                        cell_names=cell_names,
                                        cell_axis=cell_axis,
                                        gene_labels=gene_labels,
                                        allow_duplicates=allow_duplicates,
                                        genome=genome,
                                        metadata_channels=metadata_channels)
    transform_fn, transform_kws = check_transform_args(transform=transform,
                                                       pseudocount=pseudocount,
                                                       cofactor=cofactor)

    # load data
    # example: scprep.io.load_csv("data.csv")
    # https://scprep.readthedocs.io/en/stable/reference.html#module-scprep.io
    data = load_fn(filename, **load_kws)
    data = scprep.sanitize.check_numeric(data, copy=True)

    # filter data
    # https://scprep.readthedocs.io/en/stable/reference.html#module-scprep.filter
    if min_library_size is not None:
        data = scprep.filter.filter_library_size(data,
                                                 cutoff=min_library_size)
    if min_cells_per_gene is not None:
        data = scprep.filter.remove_rare_genes(data,
                                               min_cells=min_cells_per_gene)

    # normalize data
    # https://scprep.readthedocs.io/en/stable/reference.html#module-scprep.normalize
    if library_size_normalize:
        data = scprep.normalize.library_size_normalize(data)

    # transform data
    # example: data = scprep.transform.sqrt(data)
    # https://scprep.readthedocs.io/en/stable/reference.html#module-scprep.transform
    if transform is not None:
        data = transform_fn(data, **transform_kws)

    # run MAGIC
    # https://magic.readthedocs.io/
    magic_op = magic.MAGIC(knn=knn, decay=decay, t=t_magic, n_pca=n_pca,
                           knn_dist=knn_dist,
                           n_jobs=n_jobs, random_state=random_state,
                           verbose=verbose)
    magic_data = magic_op.fit_transform(data, genes=genes)

    # save as csv
    magic_data = pd.DataFrame(magic_data)
    if cell_axis in ['col', 'column']:
        magic_data = magic_data.T
    magic_data.to_csv(output)
