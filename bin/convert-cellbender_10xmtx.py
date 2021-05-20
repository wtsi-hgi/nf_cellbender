#!/usr/bin/env python


__author__ = 'Henry Taylor'
__date__ = '2020-07-09'
__version__ = '0.0.1'

import argparse
import os
from distutils.version import LooseVersion
from typing import Dict
import tables
import scipy
import scipy.io
import scipy.sparse
import gzip
import pandas as pd
import numpy as np
import anndata
import scanpy as sc


def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(
    file: str,
    analyzed_barcodes_only: bool = True
) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count
            matrix. True to load a limited set of barcodes: only those
            analyzed by the algorithm. This allows relevant latent
            variables to be loaded properly into adata.obs and adata.obsm,
            rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """
    d = dict_from_h5(file)
    X = scipy.sparse.csc_matrix(
        (d.pop('data'), d.pop('indices'), d.pop('indptr')),
        shape=d.pop('shape')
    ).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print(
                'Warning: analyzed_barcodes_only=True, but the key ',
                '"barcodes_analyzed_inds" or "barcode_indices_for_latents" ',
                'is missing from the h5 file. ',
                'Will output all barcodes, and proceed as if ',
                'analyzed_barcodes_only=False'
            )

    print(d.keys())

    # Construct the count matrix.
    if 'gene_names' in d.keys():
        gene_symbols = d.pop('gene_names').astype(str)
    else:
        gene_symbols = d.pop('name').astype(str)
    adata = anndata.AnnData(
        X=X,
        obs={'barcode': d.pop('barcodes').astype(str)},
        var={
            'gene_ids': d.pop('id').astype(str),
            'gene_symbols': gene_symbols
        }
    )
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_ids', inplace=True)

    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print(
                'Unable to load data into AnnData: ', key, value, type(value)
            )

    if analyzed_barcodes_only:
        cols = adata.obs.columns[
            adata.obs.columns.str.startswith('barcodes_analyzed')
            | adata.obs.columns.str.startswith('barcode_indices')
        ]
        for col in cols:
            try:
                del adata.obs[col]
            except Exception:
                pass

    return adata


def cellbender_to_tenxmatrix(
    adata,
    out_file='',
    out_dir='tenx_from_adata',
    verbose=True
):
    """Write 10x like data from 10x H5.

    Parameters
    ----------
    adata : pandas.DataFrame
        Description of parameter `adata`.
    output_file : string
        Description of parameter `output_file`.
    out_dir : string
        Description of parameter `out_dir`.
    verbose : boolean
        Description of parameter `verbose`.

    Returns
    -------
    execution_code : int
    """
    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir, exist_ok=True)

    # Set up out_file
    if out_file != '':
        out_file = '{}-'.format(out_file)

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # First filter out any cells that have 0 total counts
    zero_count_cells = adata.obs_names[np.where(adata.X.sum(axis=1) == 0)[0]]
    if verbose:
        print("Filtering {}/{} cells with 0 counts.".format(
            len(zero_count_cells),
            adata.n_obs
        ))
    adata = adata[adata.obs_names.difference(zero_count_cells, sort=False)]

    # Save the barcodes.
    out_f = os.path.join(
        out_dir,
        '{}barcodes.tsv.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
    pd.DataFrame(adata.obs.index).to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the features.
    # CellBender output: gene_symbols as index, enseble as `gene_ids`
    # Type is stored in adata.var.type
    out_f = os.path.join(
        out_dir,
        '{}features.tsv.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))

    # Check if read in by gene_symbol or gene_id
    if 'gene_ids' in adata.var.columns:
        gene_var = 'gene_ids'
    elif 'gene_symbols' in adata.var.columns:
        gene_var = 'gene_symbols'
    # elif 'id' in adata.var.columns:
    #     gene_var = 'id'
    # elif 'name' in adata.var.columns:
    #     gene_var = 'name'
    else:
        raise Exception(
            'Could not find "gene_symbols" or "gene_ids" in adata.var'
        )
    df_features = pd.DataFrame(
        data=[
            adata.var.index.values,
            adata.var.loc[:, gene_var].values,
            ["Gene Expression"] * adata.n_vars
        ]
    ).T

    df_features.to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the count-adjusted matrix
    out_mtx = adata.X.transpose()
    if not isinstance(out_mtx, scipy.sparse.csr.csr_matrix):
        out_mtx = scipy.sparse.csr_matrix(out_mtx)
    out_f = os.path.join(
        out_dir,
        '{}matrix.mtx.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
        # print(out_mtx.dtype)  # it looks like even counts stored as float
    with gzip.open(out_f, 'wb', compresslevel=9) as f:
        scipy.io.mmwrite(
            f,
            out_mtx,
            comment='metadata_json: {"format_version": 2}'
            # field='real'  # can be integer if counts otherwise real
        )

    if verbose:
        print(
            'Done.'
        )

    return 0


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read an H5 10x-formatted object and write matrix files similar to
            10x output. This script requires pandas >1.0.1
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_10x',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '-g', '--genome',
        action='store',
        dest='g',
        default="background_removed",
        help='Filter expression to genes within this genome. For CellBender,\
        this will be "background_removed"'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files.\
            (default: no tag basename)'
    )

    parser.add_argument(
        '-od', '--output_dir',
        action='store',
        dest='od',
        default='tenx_from_adata',
        help='Basename of output directory.\
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the AnnData file
    # adata = sc.read_10x_h5(filename=options.h5, genome=options.g)
    # Load remove-background output data.
    # We need to do this because of bug here:
    # https://github.com/broadinstitute/CellBender/issues/57
    adata = anndata_from_h5(
        options.h5,
        analyzed_barcodes_only=True
    )
    print(adata)

    # Run the conversion function.
    _ = cellbender_to_tenxmatrix(
        adata,
        out_file=options.of,
        out_dir=options.od
    )


if __name__ == '__main__':
    main()
