import logging
from pathlib import Path
from argparse import ArgumentParser
from typing import Tuple

from tqdm import tqdm
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(logging.StreamHandler())

# limits denote the lower left and upper right corner of a selection region
RA_LIMITS = (60, 195)
DEC_LIMITS = (-32.7, -20.7)

def check_row(
    img_shape: Tuple[int], 
    row_count: int, 
    img_wcs: WCS, 
    ra_limits: Tuple[float, float]=RA_LIMITS,
    dec_limits: Tuple[float, float]=DEC_LIMITS
    ) -> np.ndarray:
    """Convert pixel-coordinates into sky positions and see if they fall in an allowed region 

    Args:
        img_shape (Tuple[int]): Image dimension of data being considered
        row_count (int): Row to consider in of the image data
        img_wcs (WCS): Coordinate transform for the data
        ra_limits (Tuple[float, float], optional): Lower and upper RA bounds for a lower-left and upper-right corner. Defaults to RA_LIMITS.
        dec_limits (Tuple[float, float], optional): Lower and upper Dec bounds for a lower-left and upper-right corner. Defaults to DEC_LIMITS.

    Returns:
        np.ndarray: _description_
    """

    row_val = np.zeros(img_shape[1]) + row_count
    col_val = np.arange(img_shape[1])

    sky_coords = img_wcs.all_pix2world(
        np.array((col_val, row_val)).T,
        0,
        ra_dec_order=True
    )

    ra_vals = sky_coords[:,0]
    dec_vals = sky_coords[:,1]

    dec_mask = (dec_limits[0] <= dec_vals) & (dec_vals < dec_limits[1])
    
    logger.debug(sky_coords)

    if ra_limits[0] < ra_limits[1]:
        ra_mask = (ra_limits[0] <= ra_vals) & (ra_vals < ra_limits[1])
    else:
        # Wraps around 360 -> 0
        ra_mask = (ra_limits[0] <= ra_vals) | (ra_vals < ra_limits[1])

    mask = ra_mask & dec_mask
    logger.debug(f"{row_count=}, {img_shape=}, Pixels in region: {np.sum(mask)} - {ra_vals.shape} - ra={np.sum(ra_mask)} - dec={np.sum(dec_mask)}")

    return mask


def save_fits_img(data: np.ndarray, hdr: dict, path: Path=Path('test.fits')) -> None:
    """Simple wrapper to save fits files

    Args:
        data (np.ndarray): Data to write to a fits file
        hdr (dict): The header of the corresponding data
        path (Path, optional): Where to save data to. Defaults to 'test.fits'.
    """
    logger.info(f"Writing mask to {path} with shape {data.shape=}")
    fits.writeto(
        path,
        data.astype(float),
        hdr,
        overwrite=True
    )


def construct_bounding_region(mask: np.ndarray, header: dict=None):
    assert len(mask.shape) == 2, "Only two dimensions are supported"

    logger.info("Finding row min")
    row_min = 0
    while np.all(~mask[row_min, :]):
        row_min += 1
    
    logger.info("Finding row max")
    row_max = 0
    while np.all(~mask[-row_max, :]):
        row_max += 1
    
    logger.info("Finding col min")
    col_min = 0
    while np.all(~mask[:, col_min]):
        col_min += 1
    
    logger.info("Finding col max")
    col_max = 0
    while np.all(~mask[:, -col_max]):
        col_max += 1

    logger.info(f"{row_min=} {row_max=} {col_min=} {col_max}")

    if header is not None:
        header['CRPIX1'] -= col_min 
        header['CRPIX2'] -= row_min 
    
    row_slice = slice(row_min, -row_max)
    col_slice = slice(col_min, -col_max)
    hdr_slice = header

    logger.debug(f"{row_slice=} {col_slice}")

    return row_slice, col_slice, hdr_slice


def nan_image(img: Path, trim: bool=False, save_mask: bool=False):
    logger.info(f"Opening {img}")

    with fits.open(img) as img_fits:
        img_hdr = img_fits[0].header
        img_data = img_fits[0].data

    img_wcs = WCS(img_hdr)

    logger.info(f"Loaded data shape: {img_data.shape}")

    logger.debug(img_wcs)
    logger.debug(f"{img_data.shape=}")

    row_masks = [check_row(img_data.shape, i, img_wcs) for i in tqdm(range(img_data.shape[0]))]
    img_mask = np.vstack(row_masks)

    # Apply mask to the image here
    img_data[~img_mask] = np.nan

    if save_mask:
        save_fits_img(
            img_mask,
            img_hdr,
            path=img.with_name(f"{img.stem}_mask.fits")
        )

    if not trim:
        return 

    logger.info("Applying trim operation")
    slices = construct_bounding_region(img_mask, header=img_hdr)
    row_slice, col_slice, hdr_slice = slices

    save_fits_img(
        img_data[row_slice, col_slice],
        hdr_slice,
        path=img.with_name(f"{img.stem}_nan.fits")
    )

    if save_mask:
        save_fits_img(
            img_mask[row_slice, col_slice],
            hdr_slice,
            path=img.with_name(f"{img.stem}_mask_sliced.fits")
        )



if __name__ == '__main__':
    parser = ArgumentParser(description='Apply a NaN mask to the region outside the DR1 prime area')

    parser.add_argument('image', type=Path, help='Path to image to apply nan masking to')
    parser.add_argument('-v', '--verbose', action='store_true', help='More logging')
    parser.add_argument('-t','--no-trim', action='store_true', help='Disable the trim operaiton applied to the extracted region')
    parser.add_argument('-m', '--save-mask', action='store_true', help='Save all the masks for QA')

    args = parser.parse_args()

    if args.verbose:
        logger.setLevel(logging.DEBUG)

    nan_image(
        args.image,
        trim=not args.no_trim,
        save_mask=args.save_mask
    )