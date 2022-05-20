#!/usr/bin/env bash

"""A small helper script to select, rename and reorder columns 
from the priorized catalogue table into a form described in the 
GLEAM-X survey description paper. 
"""

from argparse import ArgumentParser 
from pathlib import Path
from typing import List, Union
import logging 
import sys 


from astropy.table import Table
import numpy as np

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

# The frequency label mapping between the columns in the DR1 catalogue
# paper and those in the priorized fitting code
PAPER_TO_PRIOR_FREQ = {
    # The reference 60MHz image
    'wide': 'wide',
    # The narrow band
    '076' : '072-080',
    '084' : '080-088',
    '092' : '088-095',
    '099' : '095-103',
    '107' : '103-111',
    '115' : '111-118',
    '122' : '118-126',
    '130' : '126-134',
    '143' : '139-147',
    '151' : '147-154',
    '158' : '154-162',
    '166' : '162-170',
    '174' : '170-177',
    '181' : '177-185',
    '189' : '185-193',
    '197' : '193-200',
    '204' : '200-208',
    '212' : '208-216',
    '220' : '216-223',
    '227' : '223-231',
    # The 30MHz wide band
    'W_087' : '072-103',
    'W_118' : '103-134',
    'W_154' : '139-170',
    'W_185' : '170-200',
    'W_215' : '200-231'
}

# Split them out
PAPER_FREQ_NAMES = PAPER_TO_PRIOR_FREQ.keys()
PRIOR_FREQ_NAMES = PAPER_TO_PRIOR_FREQ.values()

# The prefix for each column in the paper catalogue and the corresponding
# item in the priorized fitted catalogue. Note that the priorized catalogue
# does not have ['err_abs_flux_pct', 'err_fit_flux_pct'] columns. 
PAPER_TO_PRIOR_NAMES = {
    'background' : 'background',
    'local_rms' : 'local_rms',
    'peak_flux' : 'peak_flux',
    'err_peak_flux' : 'err_peak_flux',
    'int_flux' : 'int_flux',
    'err_int_flux' : 'err_int_flux',
    'a' : 'a',
    'err_a' : 'err_a',
    'b' : 'b',
    'err_b' : 'err_b',
    'pa' : 'pa',
    'err_pa' : 'err_pa',
    'residual_mean' : 'residual_mean',
    'residual_std' : 'residual_std',
    'psf_a' : 'psf_a',
    'psf_b' : 'psf_b',
    'psf_pa' : 'psf_pa',
    'Name' : 'src_name',
    'sp_int_flux_fit_200' : 'pl_norm', 
    'err_sp_int_flux_fit_200' : 'pl_norm_err',
    'sp_alpha' : 'pl_alpha',
    'err_sp_alpha' : 'pl_alpha_err',
    'sp_reduced_chi2' : 'pl_rchi2',
    'csp_int_flux_fit_200' : 'cpl_norm',
    'err_csp_int_flux_fit_200' : 'cpl_norm_err',
    'csp_alpha' : 'cpl_alpha',
    'err_csp_alpha' : 'cpl_alpha_err',
    'csp_beta' : 'cpl_q',
    'err_csp_beta' : 'cpl_q_err',
    'csp_reduced_chi2' : 'cpl_rchi2',
    'ra_str' : 'ra_str',
    'dec_str' : 'dec_str',
    'RAJ2000' : 'ref_ra',
    'DEJ2000' : 'ref_dec',
    'err_RAJ2000' : 'err_ra',
    'err_DEJ2000' : 'err_dec',
}

# These are the names of the columns required for each frequency
PAPER_COL_FREQ = [ 
      'background', 'local_rms', 'peak_flux','err_peak_flux','int_flux','err_int_flux','a','b','pa','residual_mean','residual_std','psf_a','psf_b','psf_pa'
]

# For my own sanity and to cover myself
assert all([i in PAPER_TO_PRIOR_NAMES for i in PAPER_COL_FREQ]), "A column from does not have an expected mapping"

def create_paper_columns() -> List[str]:
    """Generate the column names as described in the GLEAM-X survey paper

    Returns:
        List[str]: Column names of the final gleam-x table
    """
    paper_freq_names = PAPER_FREQ_NAMES
    paper_col_freq = PAPER_COL_FREQ

    # Creates the initial mapping
    paper_columns = [f"{q}_{f}" for f in paper_freq_names for q in paper_col_freq]
    
    # Now insert extra non-freq dependent items
    paper_columns.insert(0, 'Name')
    paper_columns.insert(3, 'ra_str')
    paper_columns.insert(4, 'dec_str')
    paper_columns.insert(5, 'RAJ2000')
    paper_columns.insert(6, 'err_RAJ2000')
    paper_columns.insert(7, 'DEJ2000')
    paper_columns.insert(8, 'err_DEJ2000')
    paper_columns.insert(14, 'err_a_wide')
    paper_columns.insert(16, 'err_b_wide')
    paper_columns.insert(18, 'err_pa_wide')
    paper_columns.insert(22, 'err_abs_flux_pct')
    paper_columns.insert(23, 'err_fit_flux_pct')

    # Now insert the items at the tail end
    items = [
        'sp_int_flux_fit_200', 
        'err_sp_int_flux_fit_200',
        'sp_alpha',
        'err_sp_alpha',
        'sp_reduced_chi2',
        'csp_int_flux_fit_200',
        'err_csp_int_flux_fit_200',
        'csp_alpha',
        'err_csp_alpha',
        'csp_beta',
        'err_csp_beta',
        'csp_reduced_chi2'
    ]

    paper_columns += items

    for i, item in enumerate(paper_columns):
        logger.debug(f"{i+1:3d} {item}")

    logger.debug(f"Total number of columns: {len(paper_columns)}")
    logger.debug(f"Number of freq_names: {len(paper_freq_names)}")
    logger.debug(f"Number of col_freq names: {len(paper_col_freq)}")

    return paper_columns


def priorized_catalogue_lookup(col: str) -> str:
    """Given a column name from the paper table, generate the corresponding
    column name from the set of priorized table column names

    Args:
        col (str): The name from the paper table that needs to be looked up

    Raises:
        ValueError: Raised when no mapping value is known

    Returns:
        str: Corresponding column name in the priorized table
    """

    mapped_col = None

    if col in PAPER_TO_PRIOR_NAMES:
        return PAPER_TO_PRIOR_NAMES[col]

    # Currently these columns have no equivilant in the 
    # prior catalogue. Putting them here explicitly. 
    if col in ['err_abs_flux_pct', 'err_fit_flux_pct']:
        return col

    if '_' in col:
        splits = col.split('_')
        
        # Horrid hack to avoid the requested W for wideband lookup
        split_to = -2 if splits[-2] == 'W' else -1
        paper_freq = '_'.join(splits[split_to:])
        paper_prefix = '_'.join(splits[:split_to])
        
        if paper_freq in PAPER_TO_PRIOR_FREQ:
            prior_freq = PAPER_TO_PRIOR_FREQ[paper_freq]
            prior_prefix = PAPER_TO_PRIOR_NAMES[paper_prefix]
            
            if paper_freq in ['W_087', 'W_118', 'W_154', 'W_185', 'W_215']:
                mapped_col = f"{prior_prefix}_W_{prior_freq.replace('-','_')}MHz"
            elif paper_freq != 'wide':
                mapped_col = f"{prior_prefix}_N_{prior_freq.replace('-','_')}MHz"
            else:
                mapped_col = f"{prior_prefix}"

    if mapped_col is None:
        raise ValueError(f"Mapped value for {col} does not exist. Current values {paper_freq=} {paper_prefix=}")

    return mapped_col


def create_mapping() -> tuple:
    """Generates the set of columns for the paper and priorized tables

    Returns:
        tuple: Tuple containing the two tuples of column names 
    """
    paper_columns = tuple(create_paper_columns())
    prior_columns = tuple([priorized_catalogue_lookup(pc) for pc in paper_columns])

    return paper_columns, prior_columns

def print_mapping() -> None:
    """Print out the column mapping from the paper to the priorized table
    """
    logger.info("Returning the column name mapping")
    paper_columns, prior_columns = create_mapping()

    for prior_col, paper_col in zip(prior_columns, paper_columns):
        logger.info(f"{paper_col} {prior_col}")


def clean_table(cata_tab: Table) -> Table:
    """When writing an astropy table fill values are used in place of values that have been masked. 
    When importing with from_pandas() a fill value replaces NaN values. This is undesirable. 

    Args:
        cata_tab (astropy.table.Table): The newly created table with bad fill values
    
    Returns:
        astropy.table.Table: The fixed table that will not insert 1e20 values when written out as fits. 
    """
    logger.info(f"Table Columns: {cata_tab.columns=}")
    logger.info(f"Table Columns: {len(cata_tab.columns)=}")

    for c in cata_tab.columns:
        try:
            logger.info(f"Cleaning column {c}, {np.sum(cata_tab[c].mask)} nans to fill")
            cata_tab[c][cata_tab[c].data.mask] = np.nan
        except Exception as e:
            logger.debug(f"Error cleaning {c} - {e}")
            pass

    for c in ('err_RAJ2000', 'err_DEJ2000', 'err_peak_flux_wide', 'err_a_wide', 'err_b_wide', 'residual_mean_wide','residual_std_wide'):
        logger.info(f"Column {c=} {cata_tab[c].dtype=}")
        try:
            # This 1e20 is a magic astropy Table fill value for nan / non-finite numbers
            max_mask = np.isclose(cata_tab[c], 1e20)
            if np.any(max_mask):
                # Using -1, as errors on these quantities are aegean issues quantifying the unceraintity on the wide band
                # It appears the sub-band quantities are OK
                cata_tab[c][max_mask] = -1 

                logger.info(f"Detected astropy table {np.sum(max_mask)} fill values in {c}, reseting to nan")
            else:
                logger.info(f"Confused Tim, {np.max(cata_tab[c])}")
        except TypeError as err:
            logger.debug(f"TypeError reported for {c}")

    return cata_tab

def apply_mapping(cata_path: Union[str,Path]) -> None:
    """Takes the priorized catalogue table and selected the subset of columns, renames them, and 
    performs some basic cleaning up of masked values to appropriately write to a FITS table

    Args:
        cata_path (str): Path to the priorized catalogue table
    """

    logger.info("Creating the column mapping")
    paper_columns, prior_columns = create_mapping()

    cata_path = Path(cata_path)

    logger.info(f"Loading {cata_path}")
    cata_df = Table.read(cata_path).to_pandas()

    # Insert temp columns
    cata_df['err_abs_flux_pct'] = 8. # Section 4.2.1 of GLEAM-X survey paper
    cata_df['err_fit_flux_pct'] = 2. # Section 4.4 of GLEAM-X survey paper

    unknown_cols = [pc for pc in prior_columns if pc not in cata_df.columns]
    for uc in unknown_cols:
        logger.info(f"Unknown column: {uc}")
    
    if bool(unknown_cols):
        logger.info(f"Unkown columns: {len(unknown_cols)}")

    assert all([i in cata_df.columns for i in prior_columns])

    # Make a subset in the order we expect
    logger.info(f"Extracting and renaming {len(prior_columns)} from {cata_path}")
    sub_cata_df = cata_df[list(prior_columns)]
    sub_cata_df.columns = paper_columns
    
    # now clean and save
    out_name = cata_path.with_name(f"{cata_path.stem}_paper.fits")
    logger.info(f"Saving to {out_name}")
    
    paper_tab = clean_table(Table.from_pandas(sub_cata_df))
    paper_tab.write(out_name, overwrite=True)


    # bad_col = paper_tab['err_DEJ2000']
    # logger.debug(
    #     f"{np.min(bad_col)=} {np.max(bad_col)=} {np.sum(np.isfinite(bad_col))} {len(bad_col)=}"
    # )

    # import ipdb
    # ipdb.set_trace()


if __name__ == '__main__':
    parser = ArgumentParser(description="Rename and subset columns into the final gleam-x catalogue")

    parser.add_argument('catalogue', type=str, help='Path to catalogue to correct')
    parser.add_argument('-v','--verbose', default=False, action='store_true', help='Extra logging')
    parser.add_argument('-m','--mapping', default=False, action='store_true', help='Print the column mapping that will be applied and exit')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    if args.mapping:
        print_mapping()
        sys.exit(0)
    else:
        apply_mapping(
            args.catalogue
        )
    