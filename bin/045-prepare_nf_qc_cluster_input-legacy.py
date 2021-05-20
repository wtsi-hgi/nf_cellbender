#!/usr/bin/env python


__author__ = 'Henry Taylor'
__date__ = '2020-07-09'
__version__ = '0.0.1'

import argparse
import os
import pandas as pd


def prepare_qc_cluster_input(
    experiment_ids,
    matrix_dirs,
    ncells_expected,
    out_file='',
    verbose=True
):
    """Write a .tsv file that contains the information needed for qc_cluster
    pipeline.

    Parameters
    ----------
    experiment_ids : array
        Description of parameter `adata`.
    matrix_dirs : array
        Description of parameter `matrix_dirs`.
    ncells_expected : array
        Description of parameter `ncells_expected`.
    output_file : string
        Description of parameter `output_file`.
    verbose : boolean
        Description of parameter `verbose`.

    Returns
    -------
    execution_code : int
    """

    # Set up out_file
    if out_file != '':
        out_file = '{}'.format(out_file)


    # First, check and make sure files are actually there
    for dir in matrix_dirs:
        if os.path.isdir(dir):
            exp_files = ['barcodes.tsv.gz', 'features.tsv.gz', 'matrix.mtx.gz']
            for file in exp_files:
                file_path = '{}/{}'.format(dir, file)
                if not os.path.isfile(file_path):
                    raise ValueError('The file "{}" does not exist.'.format(file_path))
        else:
            raise ValueError('The directory "{}" does not exist.'.format(dir))

    # Create dataframe -- we don't really care about short_experiment_id,
    # so just make it the same as experiment_ids
    df = pd.DataFrame({
        "experiment_id": experiment_ids,
        "data_path_10x_format": matrix_dirs,
        "short_experiment_id": experiment_ids,
        "ncells_expected": ncells_expected
    })

    df.to_csv(
        out_file,
        sep='\t',
        index=False,
        header=True
    )

    return 0


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Format a file to
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-id', '--experiment_ids',
        action='store',
        dest='ids',
        required=True,
        help='List of comma-separated experiment ids.'
    )

    parser.add_argument(
        '-dir', '--matrix_dirs',
        action='store',
        dest='dirs',
        required=True,
        help='List of comma-separated 10x formatted matrix dirs. Indices should\
            correspond with the experiment ids.'
    )

    parser.add_argument(
        '-n', '--ncells_expected',
        action='store',
        dest='n_cells',
        required=True,
        help='List of comma-separated number of cells expected. Indices should\
            correspond with the experiment ids.'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        required=True,
        help='Output file.'
    )

    options = parser.parse_args()

    # Split all inputs
    exp_ids = options.ids.split(",")
    dirs = options.dirs.split(",")
    ncells = options.n_cells.split(",")

    # Run the conversion function.
    _ = prepare_qc_cluster_input(
        exp_ids,
        dirs,
        ncells,
        out_file=options.of
    )


if __name__ == '__main__':
    main()
